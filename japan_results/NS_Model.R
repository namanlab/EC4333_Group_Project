library(tidyverse)
library(minpack.lm)  # For non-linear least squares fitting
library(forecast)
library(YieldCurve)
library(gridExtra)

df_pr = read.csv("df_japan_spot_rates.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD))  

f2 <- function(beta1,theta,tau){
  beta1*(1-exp(-theta/tau))/(theta/tau)
}

f3 <- function(beta2,theta,tau){
  beta2*((1-exp(-theta/tau))*tau/theta - exp(-theta/tau))
}

f4 <- function(beta3,theta,tau2){
  beta3*((1-exp(-theta/tau2))*tau2/theta - exp(-theta/tau2))
}

ns <- function(theta, beta0,beta1,beta2,tau){
  beta0 + f2(beta1,theta,tau) + f3(beta2,theta,tau)
}

extendns <- function(theta, beta0,beta1,beta2,beta3,tau,tau2){
  beta0 + f2(beta1,theta,tau) + f3(beta2,theta,tau) + f4(beta3,theta,tau2) 
}

ns_solve_fn = function(x, params){
  ns(x, params[1], params[2],params[3],params[4])
}


exns_solve_fn = function(x, params){
  extendns(x, params[1], params[2],params[3],params[4],params[5],params[6])
}


rem_outl <- function(data, time_column = "TIME_PERIOD", value_column = "OBS_VALUE", frequency = 12, threshold = 2) {
  ts_data <- ts(data[[value_column]], frequency = frequency)
  stl_fit <- stl(ts_data, s.window = "periodic")
  residuals <- stl_fit$time.series[, "remainder"]
  residual_sd <- sd(residuals, na.rm = TRUE)
  outlier_threshold <- threshold * residual_sd
  outliers <- abs(residuals) > outlier_threshold
  data <- data %>%
    mutate(Residual = residuals,
           Outlier = outliers)
  outliers_data <- data %>% filter(Outlier == F)
  return(outliers_data)
}


detect_outliers_stl <- function(data, time_column = "TIME_PERIOD", value_column = "OBS_VALUE", frequency = 12, threshold = 2) {
  ts_data <- ts(data[[value_column]], frequency = frequency)
  stl_fit <- stl(ts_data, s.window = "periodic")
  residuals <- stl_fit$time.series[, "remainder"]
  residual_sd <- sd(residuals, na.rm = TRUE)
  outlier_threshold <- threshold * residual_sd
  outliers <- abs(residuals) > outlier_threshold
  data <- data %>%
    mutate(Residual = residuals,
           Outlier = outliers)
  outliers_data <- data %>% filter(Outlier == TRUE)
  return(outliers_data)
}

fit_model_fn <- function(date_val = "2007-11-06", fn_main = exns_solve_fn,
                         init_params = c(1, -1, 1, 1, 1, 100), viz = FALSE,
                         outlier_handle = TRUE, frequency = 12, stl_threshold = 10){
  
  residuals_fn <- function(params, MAT, VAL) {
    res <- VAL - sapply(MAT, fn_main, params = params)
    return(res)
  }
  
  # Filter data for the specified date and type
  data_rel <- df_pr %>% filter(TIME_PERIOD == date_val)
  
  # Exclude STL-detected outliers from initial model fitting if outlier_handle is TRUE
  if (outlier_handle) {
    stl_outliers <- detect_outliers_stl(data_rel, time_column = "TIME_PERIOD", 
                                        value_column = "OBS_VALUE", frequency = frequency, 
                                        threshold = stl_threshold)
    data_clean <- anti_join(data_rel, stl_outliers, by = c("TIME_PERIOD", "MAT"))
    outlier_ct <- nrow(stl_outliers)
  } else {
    data_clean <- data_rel
    outlier_ct <- 0
  }
  
  train_index <- sample(1:nrow(data_clean), size = round(nrow(data_clean)*0.6), replace = F, )
  train_data <- data_clean[train_index, ]
  test_data <- data_clean[-train_index, ]
  
  # Initial model fitting on training data
  fit <- nls.lm(par = init_params, 
                fn = function(x) residuals_fn(x, MAT = train_data$MAT, VAL = train_data$OBS_VALUE))
  params_est <- fit$par
  
  # Calculate predictions and residuals for both train and test data
  train_data <- train_data %>% mutate(Predicted_Yield = sapply(MAT, fn_main, params = params_est),
                                      Residual = OBS_VALUE - Predicted_Yield)
  test_data <- test_data %>% mutate(Predicted_Yield = sapply(MAT, fn_main, params = params_est),
                                    Residual = OBS_VALUE - Predicted_Yield)
  
  # Calculate error metrics for training data
  train_mse <- mean(train_data$Residual^2)
  train_mae <- mean(abs(train_data$Residual))
  
  # Calculate error metrics for testing data
  test_mse <- mean(test_data$Residual^2)
  test_mae <- mean(abs(test_data$Residual))
  
  # Calculate R-squared and Adjusted R-squared on the training data
  rss <- sum(train_data$Residual^2)
  tss <- sum((train_data$OBS_VALUE - mean(train_data$OBS_VALUE))^2)
  r2 <- 1 - rss / tss
  n <- nrow(train_data)
  k <- length(init_params)
  adj_r2 <- 1 - (1 - r2) * ((n - 1) / (n - k - 1))
  
  # Display results
  lst_vals <- list(train_mse = train_mse, train_mae = train_mae, 
                   test_mse = test_mse, test_mae = test_mae, 
                   r2 = r2, adj_r2 = adj_r2, outlier_ct = outlier_ct, 
                   params_est = params_est)
  cat("Processed date:", date_val, "\nTrain MSE:", train_mse, "Train MAE:", train_mae, 
      "\nTest MSE:", test_mse, "Test MAE:", test_mae, 
      "\nR2:", r2, "Adjusted R2:", adj_r2, "STL Outliers:", outlier_ct, "\n")
  
  # Optional visualization
  if (viz){
    data_long <- bind_rows(train_data, test_data, .id = "Dataset") %>%
      pivot_longer(cols = c(OBS_VALUE, Predicted_Yield), 
                   names_to = "Category", 
                   values_to = "Value")
    
    # Create the plot using ggplot2
    p1 <- ggplot(data_long, aes(x = MAT, y = Value, color = Category, linetype = Category)) +
      geom_point(data = filter(data_long, Category == "OBS_VALUE"), size = 1) +
      geom_line(size = 1) +
      labs(
        title = paste("Fitted vs Observed Values on", date_val),
        subtitle = paste("Train MSE:", round(train_mse, 5), "Test MSE:", round(test_mse, 5)),
        x = "Maturity (MAT)",
        y = "Spot Rates",
        color = "Legend",
        linetype = "Legend"
      ) +
      facet_wrap(~ Dataset, ncol = 1, scales = "free_y") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      ) +
      scale_color_manual(values = c("OBS_VALUE" = "blue", "Predicted_Yield" = "red"))
    print(p1)
  }
  
  return(lst_vals)
}


res = fit_model_fn(viz = T, outlier_handle = F)
fit_model_fn(fn_main = ns_solve_fn, viz = T, outlier_handle = F)


# Function to iterate over time periods for a given data frame and store results
process_time_periods <- function(TIME_PERIOD, fn) {
  # Initialize lists to store separate results
  all_results <- list()
  params_est_results <- list()
  
  # Get unique time periods to iterate over
  unique_dates <- unique(TIME_PERIOD)
  
  for (date in unique_dates) {
    # Fit the model and get results for each date
    res <- fit_model_fn(date_val = date, fn = fn, outlier_handle = F)
    params_est <- as.list(res$params_est)
    names(params_est) <- paste0("param_", seq_along(params_est))  # Name the params_est items
    params_est$TIME_PERIOD <- date
    other_results <- res
    other_results$TIME_PERIOD <- date
    other_results$params_est <- NULL  # Remove params_est from other results
    
    # Store params_est and other results in separate lists
    params_est_results[[as.character(date)]] <- params_est
    all_results[[as.character(date)]] <- other_results
  }
  
  # Combine lists into data frames
  params_est_df <- bind_rows(params_est_results)
  all_results_df <- bind_rows(all_results)
  
  # Return both data frames as a list
  return(list(params_est_df = params_est_df, all_results_df = all_results_df))
}

df_results_exns <- process_time_periods(df_pr$TIME_PERIOD, exns_solve_fn)
df_params_exns <- df_results_exns$params_est_df
df_all_exns <- df_results_exns$all_results_df

df_results_ns <- process_time_periods(df_pr$TIME_PERIOD, ns_solve_fn)
df_params_ns <- df_results_ns$params_est_df
df_all_ns <- df_results_ns$all_results_df

df_results_exns <- df_all_exns
df_results_ns <- df_all_ns

df_results_exns_params <- df_params_exns
df_results_ns_params <- df_params_ns
# theta, beta0,beta1,beta2,beta3,tau,tau2
# theta, beta0,beta1,beta2,tau
colnames(df_results_exns_params)[1:6] <- c("beta0","beta1","beta2","beta3","tau","tau2")
colnames(df_results_ns_params)[1:4] <- c("beta0","beta1","beta2","tau")
df_results_ns_params <- df_results_ns_params[,c(1:4, 7)]

write.csv(df_results_exns, "../data/japan_df_results_exns.csv", row.names = F)
write.csv(df_results_ns, "../data/japan_df_results_ns.csv", row.names = F)
write.csv(df_results_exns_params, "../data/japan_df_results_exns_params.csv", row.names = F)
write.csv(df_results_ns_params, "../data/japan_df_results_ns_params.csv", row.names = F)



################################################################################
################################# RESULTS: ERR #################################
################################################################################

colnames(df_results_exns)
df_results_exns <- df_results_exns %>%
  filter(
    train_mse >= (quantile(train_mse, 0.25) - 1.5 * IQR(train_mse)) & train_mse <= (quantile(train_mse, 0.75) + 1.5 * IQR(train_mse)),
    train_mae >= (quantile(train_mae, 0.25) - 1.5 * IQR(train_mae)) & train_mae <= (quantile(train_mae, 0.75) + 1.5 * IQR(train_mae)),
    test_mse >= (quantile(test_mse, 0.25) - 1.5 * IQR(test_mse)) & test_mse <= (quantile(test_mse, 0.75) + 1.5 * IQR(test_mse)),
    test_mae >= (quantile(test_mae, 0.25) - 1.5 * IQR(test_mae)) & test_mae <= (quantile(test_mae, 0.75) + 1.5 * IQR(test_mae))
  )

df_results_ns <- df_results_ns %>%
  filter(
    train_mse >= (quantile(train_mse, 0.25) - 1.5 * IQR(train_mse)) & train_mse <= (quantile(train_mse, 0.75) + 1.5 * IQR(train_mse)),
    train_mae >= (quantile(train_mae, 0.25) - 1.5 * IQR(train_mae)) & train_mae <= (quantile(train_mae, 0.75) + 1.5 * IQR(train_mae)),
    test_mse >= (quantile(test_mse, 0.25) - 1.5 * IQR(test_mse)) & test_mse <= (quantile(test_mse, 0.75) + 1.5 * IQR(test_mse)),
    test_mae >= (quantile(test_mae, 0.25) - 1.5 * IQR(test_mae)) & test_mae <= (quantile(test_mae, 0.75) + 1.5 * IQR(test_mae))
  )

plot_line_fn_errs <- function(col_sel){
  df_plot <- df_results_exns %>%
    select(all_of(col_sel), TIME_PERIOD) %>%
    rem_outl(time_column = "TIME_PERIOD", value_column = col_sel) %>%
    mutate(TIME_PERIOD = as.Date(TIME_PERIOD))
  p1 = df_plot %>%
    ggplot() +
    annotate("rect", xmin = ym("2019-05"), xmax = ym("2020-02"), ymin = -Inf, ymax = Inf,
             fill = "red", alpha = 0.6) +
    annotate("rect", xmin = ym("2016-03"), xmax = ym("2016-10"), ymin = -Inf, ymax = Inf,
             fill = "red", alpha = 0.6) +
    geom_line(aes(x = TIME_PERIOD, y = .data[[col_sel]])) +
    theme_bw() +
    labs(x = "", y = col_sel, title = "Extended NSS") +
    theme(legend.position = "bottom") +
    scale_color_viridis_d()
  
  df_plot <- df_results_ns %>%
    select(all_of(col_sel), TIME_PERIOD) %>%
    rem_outl(time_column = "TIME_PERIOD", value_column = col_sel) %>%
    mutate(TIME_PERIOD = as.Date(TIME_PERIOD))
  p2 <- df_plot %>%
    ggplot() +
    annotate("rect", xmin = ym("2019-05"), xmax = ym("2020-02"), ymin = -Inf, ymax = Inf,
             fill = "red", alpha = 0.6) +
    annotate("rect", xmin = ym("2016-03"), xmax = ym("2016-10"), ymin = -Inf, ymax = Inf,
             fill = "red", alpha = 0.6) +
    geom_line(aes(x = TIME_PERIOD, y = .data[[col_sel]])) +
    theme_bw() +
    labs(x = "", y = col_sel, title = "NS")  +
    theme(legend.position = "bottom") +
    scale_color_viridis_d()
  
  pf1 = grid.arrange(p2, p1, nrow = 1)
  
  dfn1 = df_results_exns %>%
    filter(TIME_PERIOD < ym("2020-02"),
           TIME_PERIOD > ym("2019-05"))
  dfn2 = df_results_exns %>%
    filter(TIME_PERIOD < ym("2016-10"),
           TIME_PERIOD > ym("2016-03"))
  dfnexns = rbind(dfn1, dfn2)
  dfpexns = setdiff(df_results_exns, dfnexns)
  dfnexns$period = "Negative"
  dfpexns$period = "Positive"
  p3 <- rbind(dfnexns, dfpexns) %>%
    ggplot() +
    geom_boxplot(aes(x = period, y = .data[[col_sel]], fill = period), outliers = F) +
    theme_bw() + guides(fill = "none") + labs(title = "Extended NSS")
  
  dfn1 = df_results_ns %>%
    filter(TIME_PERIOD < ym("2020-02"),
           TIME_PERIOD > ym("2019-05"))
  dfn2 = df_results_ns %>%
    filter(TIME_PERIOD < ym("2016-10"),
           TIME_PERIOD > ym("2016-03"))
  dfnns = rbind(dfn1, dfn2)
  dfpns = setdiff(df_results_ns, dfnns)
  dfnns$period = "Negative"
  dfpns$period = "Positive"
  p4 <- rbind(dfnns, dfpns) %>%
    ggplot() +
    geom_boxplot(aes(x = period, y = .data[[col_sel]], fill = period), outliers = F) +
    theme_bw() + guides(fill = "none") + labs(title = "NS")
  
  pf2 = grid.arrange(p4, p3, nrow = 1)
  return(list(pf1, pf2))
}

results_train_mse <- plot_line_fn_errs("train_mse")
results_test_mse <- plot_line_fn_errs("test_mse")
results_train_mae <- plot_line_fn_errs("train_mae")
results_test_mae <- plot_line_fn_errs("test_mae")
results_r2 <- plot_line_fn_errs("r2")
results_adj_r2 <- plot_line_fn_errs("adj_r2")


################################################################################
################################# RESULTS: PAR #################################
################################################################################


colnames(df_results_exns_params)
df_results_exns_params <- df_results_exns_params %>%
  mutate(`beta0 + beta1` = beta1 + beta0) %>%
  filter(
    beta0 > quantile(beta0, 0.25, na.rm = TRUE) - 1.5 * IQR(beta0, na.rm = TRUE),
      beta0 < quantile(beta0, 0.75, na.rm = TRUE) + 1.5 * IQR(beta0, na.rm = TRUE),
      beta1 > quantile(beta1, 0.25, na.rm = TRUE) - 1.5 * IQR(beta1, na.rm = TRUE),
      beta1 < quantile(beta1, 0.75, na.rm = TRUE) + 1.5 * IQR(beta1, na.rm = TRUE),
      beta2 > quantile(beta2, 0.25, na.rm = TRUE) - 1.5 * IQR(beta2, na.rm = TRUE),
      beta2 < quantile(beta2, 0.75, na.rm = TRUE) + 1.5 * IQR(beta2, na.rm = TRUE),
      beta3 > quantile(beta3, 0.25, na.rm = TRUE) - 1.5 * IQR(beta3, na.rm = TRUE),
      beta3 < quantile(beta3, 0.75, na.rm = TRUE) + 1.5 * IQR(beta3, na.rm = TRUE),
      tau > quantile(tau, 0.25, na.rm = TRUE) - 1.5 * IQR(tau, na.rm = TRUE),
      tau < quantile(tau, 0.75, na.rm = TRUE) + 1.5 * IQR(tau, na.rm = TRUE),
      tau2 > quantile(tau2, 0.25, na.rm = TRUE) - 1.5 * IQR(tau2, na.rm = TRUE),
      tau2 < quantile(tau2, 0.75, na.rm = TRUE) + 1.5 * IQR(tau2, na.rm = TRUE)
  )
df_results_ns_params <- df_results_ns_params %>%
  mutate(`beta0 + beta1` = beta1 + beta0) %>%
  filter(
    beta0 > quantile(beta0, 0.25, na.rm = TRUE) - 1.5 * IQR(beta0, na.rm = TRUE),
    beta0 < quantile(beta0, 0.75, na.rm = TRUE) + 1.5 * IQR(beta0, na.rm = TRUE),
    beta1 > quantile(beta1, 0.25, na.rm = TRUE) - 1.5 * IQR(beta1, na.rm = TRUE),
    beta1 < quantile(beta1, 0.75, na.rm = TRUE) + 1.5 * IQR(beta1, na.rm = TRUE),
    beta2 > quantile(beta2, 0.25, na.rm = TRUE) - 1.5 * IQR(beta2, na.rm = TRUE),
    beta2 < quantile(beta2, 0.75, na.rm = TRUE) + 1.5 * IQR(beta2, na.rm = TRUE),
    tau > quantile(tau, 0.25, na.rm = TRUE) - 1.5 * IQR(tau, na.rm = TRUE),
    tau < quantile(tau, 0.75, na.rm = TRUE) + 1.5 * IQR(tau, na.rm = TRUE)
  )


plot_line_fn_params <- function(col_sel){
  
  df_plot <- df_results_exns_params %>%
    select(all_of(col_sel), TIME_PERIOD,) %>%
    rem_outl(time_column = "TIME_PERIOD", value_column = col_sel) %>%
    mutate(TIME_PERIOD = as.Date(TIME_PERIOD))
  p1 = df_plot %>%
    ggplot() +
    annotate("rect", xmin = ym("2019-05"), xmax = ym("2020-02"), ymin = -Inf, ymax = Inf,
             fill = "red", alpha = 0.6) +
    annotate("rect", xmin = ym("2016-03"), xmax = ym("2016-10"), ymin = -Inf, ymax = Inf,
             fill = "red", alpha = 0.6) +
    geom_line(aes(x = TIME_PERIOD, y = .data[[col_sel]])) +
    theme_bw() +
    labs(x = "", y = col_sel, title = "Extended NSS") +
    theme(legend.position = "bottom") +
    scale_color_viridis_d()
  if (col_sel %in% c("beta3", "tau2")){
    return(list(p1))
  }
  
  df_plot <- df_results_ns_params %>%
    select(all_of(col_sel), TIME_PERIOD) %>%
    rem_outl(time_column = "TIME_PERIOD", value_column = col_sel) %>%
    mutate(TIME_PERIOD = as.Date(TIME_PERIOD))
  p2 <- df_plot %>%
    ggplot() +
    annotate("rect", xmin = ym("2019-05"), xmax = ym("2020-02"), ymin = -Inf, ymax = Inf,
             fill = "red", alpha = 0.6) +
    annotate("rect", xmin = ym("2016-03"), xmax = ym("2016-10"), ymin = -Inf, ymax = Inf,
             fill = "red", alpha = 0.6) +
    geom_line(aes(x = TIME_PERIOD, y = .data[[col_sel]])) +
    theme_bw() +
    labs(x = "", y = col_sel, title = "NS")  +
    theme(legend.position = "bottom") +
    scale_color_viridis_d()
  
  pf1 = grid.arrange(p2, p1, nrow = 1)
  
  
  return(list(pf1))
}

results_beta0 <- plot_line_fn_params("beta0")
results_beta1 <- plot_line_fn_params("beta1")
results_beta10 <- plot_line_fn_params("beta0 + beta1")
results_beta2 <- plot_line_fn_params("beta2")
results_tau <- plot_line_fn_params("tau")
results_beta3 <- plot_line_fn_params("beta3")
results_tau2 <- plot_line_fn_params("tau2")


