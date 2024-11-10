library(tidyverse)
library(minpack.lm)  # For non-linear least squares fitting
library(forecast)
library(YieldCurve)
library(gridExtra)

df_pr_aaa = read.csv("df_spot_rates_aaa.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD))  %>%
  rename("AAA" = OBS_VALUE)
df_pr_all = read.csv("df_spot_rates_all.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD)) %>%
  rename("All" = OBS_VALUE)
df_pr = df_pr_aaa %>% inner_join(df_pr_all) %>%
  pivot_longer(3:4, names_to = "Type", values_to = "OBS_VALUE")


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

fit_vasicek_fn <- function(date_val = "2005-09-09", type = "AAA", split_ratio = 0.5,
                           init_params = c(0.1, 0.03, 0.02), viz = FALSE,
                           outlier_handle = TRUE, nmax = 50, frequency = 12, stl_threshold = 10) {
  
  # Residual function for Vasicek model
  residuals_vasicek <- function(params, MAT, VAL) {
    kappa <- params[1]
    theta <- params[2]
    sigma <- params[3]
    predicted_yield <- function(T) theta + (MAT - theta) * exp(-kappa * T) + (sigma^2 / (2 * kappa^2)) * (1 - exp(-2 * kappa * T))
    res <- VAL - sapply(MAT, predicted_yield)
    return(res)
  }
  
  # Filter data for the specified date and type
  data_rel <- df_pr %>% filter(TIME_PERIOD == date_val, Type == type)
  
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
  
  # Split data into training and testing sets
  set.seed(123)
  train_index <- sample(seq_len(nrow(data_clean)), size = split_ratio * nrow(data_clean))
  train_data <- data_clean[train_index, ]
  test_data <- data_clean[-train_index, ]
  
  # Initial model fitting on training data
  fit <- nls.lm(par = init_params, 
                fn = function(x) residuals_vasicek(x, MAT = train_data$MAT, VAL = train_data$OBS_VALUE),
                niter(nmax))
  params_est <- fit$par
  kappa_est <- params_est[1]
  theta_est <- params_est[2]
  sigma_est <- params_est[3]
  
  # Calculate predictions and residuals for both train and test data
  train_data <- train_data %>% mutate(
    Predicted_Yield = sapply(MAT, function(T) theta_est + (T - theta_est) * exp(-kappa_est * T) + (sigma_est^2 / (2 * kappa_est^2)) * (1 - exp(-2 * kappa_est * T))),
    Residual = OBS_VALUE - Predicted_Yield
  )
  test_data <- test_data %>% mutate(
    Predicted_Yield = sapply(MAT, function(T) theta_est + (T - theta_est) * exp(-kappa_est * T) + (sigma_est^2 / (2 * kappa_est^2)) * (1 - exp(-2 * kappa_est * T))),
    Residual = OBS_VALUE - Predicted_Yield
  )
  
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
                   r2 = r2, adj_r2 = adj_r2, outlier_ct = outlier_ct)
  cat("Processed date:", date_val, "\nTrain MSE:", train_mse, "Train MAE:", train_mae, 
      "\nTest MSE:", test_mse, "Test MAE:", test_mae, 
      "\nR2:", r2, "Adjusted R2:", adj_r2, "STL Outliers:", outlier_ct, "\n")
  
  # Optional visualization
  if (viz) {
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


fit_vasicek_fn(viz = T, split_ratio = 0.2, type = "All")
