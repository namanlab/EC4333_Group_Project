library(tidyverse)
library(minpack.lm)  # For non-linear least squares fitting
library(forecast)
library(YieldCurve)
library(gridExtra)
rm("ns")
library(splines)  # For natural cubic splines


# Load data
df_pr_aaa = read.csv("data/df_spot_rates_aaa.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD)) %>%
  rename("AAA" = OBS_VALUE)
df_pr_all = read.csv("data/df_spot_rates_all.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD)) %>%
  rename("All" = OBS_VALUE)
df_pr = df_pr_aaa %>% inner_join(df_pr_all) %>%
  pivot_longer(3:4, names_to = "Type", values_to = "OBS_VALUE")

# Outlier removal functions remain the same
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

# Outlier detection using STL decomposition
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

# Modified fit_model_fn to use natural cubic spline
fit_model_fn <- function(date_val = "2005-09-09",
                         type = "AAA",
                         num_knots = 3,  # Number of knots for spline
                         viz = FALSE,
                         outlier_handle = TRUE,
                         frequency = 12,
                         stl_threshold = 10) {
  
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
  
  # Select training and testing data based on specific maturities
  train_index <- which(data_rel$MAT %in% c(3, 6, 9, c(1, 5, 10, 20, 30) * 12))
  train_data <- data_clean[train_index, ]
  test_data <- data_clean[-train_index, ]
  
  # Fit natural cubic spline with specified number of knots
  spline_model <- lm(OBS_VALUE ~ ns(MAT, df = num_knots), data = train_data)
  
  # Predictions and residuals for both train and test data
  train_data <- train_data %>%
    mutate(Predicted_Yield = predict(spline_model, newdata = train_data),
           Residual = OBS_VALUE - Predicted_Yield)
  test_data <- test_data %>%
    mutate(Predicted_Yield = predict(spline_model, newdata = test_data),
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
  k <- num_knots
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



fit_model_fn(viz = T, type = "All")
fit_model_fn(date_val = "2019-05-02", num_knots = 4,
             outlier_handle = T, viz = TRUE)
fit_model_fn(date_val = "2019-05-02", fn_main = ns_solve_fn, 
             outlier_handle = T, viz = TRUE)


# Function to iterate over time periods and different knot configurations, recording results
process_time_periods <- function(TIME_PERIOD, type = "AAA", knot_range = 1:4) {
  # Initialize list to store results for each combination of time period and knots
  all_results <- list()
  
  # Get unique time periods to iterate over
  unique_dates <- unique(TIME_PERIOD)
  
  for (date in unique_dates) {
    for (knots in knot_range) {
      # Fit the model for each date and knot configuration
      res <- fit_model_fn(date_val = date, type = type, num_knots = knots, 
                          outlier_handle = FALSE, frequency = 12, stl_threshold = 10)
      
      # Record the results for this date and knot configuration
      result_entry <- res
      result_entry$TIME_PERIOD <- date
      result_entry$Knots <- knots  # Add knots information
      
      # Append to the results list
      all_results[[paste0(date, "_", knots)]] <- result_entry
    }
  }
  
  # Combine list into a data frame
  all_results_df <- bind_rows(all_results)
  
  # Return the data frame with all results
  return(all_results_df)
}


# For Spline model with type = "AAA"
df_all_spline_aaa <- process_time_periods(df_pr$TIME_PERIOD, type = "AAA")
df_all_spline_aaa$Type <- "AAA"

# For Spline model with type = "All"
df_all_spline_all <- process_time_periods(df_pr$TIME_PERIOD, type = "All")
df_all_spline_all$Type <- "All"

df_results_spline <- rbind(df_all_spline_aaa, df_all_spline_all)

write.csv(df_results_spline, "data/df_results_spline.csv", row.names = F)
df_results_spline = read.csv("data/df_results_spline.csv")


################################################################################
################################# RESULTS: ERR #################################
################################################################################

colnames(df_results_spline)
df_results_spline <- df_results_spline %>%
  filter(train_mse <= 0.1, train_mae <= 0.1, test_mse <= 0.1, test_mae <= 0.1)
df_results_spline <- df_results_spline %>%
  filter(train_mse <= 0.1, train_mae <= 0.1, test_mse <= 0.1, test_mae <= 0.1)

plot_line_fn_errs <- function(col_sel){
  df_plot <- df_results_spline %>%
    select(all_of(col_sel), TIME_PERIOD, Type, Knots) %>%
    rem_outl(time_column = "TIME_PERIOD", value_column = col_sel) %>%
    mutate(TIME_PERIOD = as.Date(TIME_PERIOD))
  p1 = df_plot %>%
    ggplot() +
    annotate("rect", xmin = ym("2019-05"), xmax = ym("2022-02"), ymin = 0, ymax = Inf,
             fill = "red", alpha = 0.6) +
    annotate("rect", xmin = ym("2016-06"), xmax = ym("2016-11"), ymin = 0, ymax = Inf,
             fill = "red", alpha = 0.6) +
    geom_line(aes(x = TIME_PERIOD, y = .data[[col_sel]], color = Type)) +
    theme_bw() +
    labs(x = "", y = col_sel, title = "Natural Cubic Splines") +
    theme(legend.position = "bottom") +
    scale_color_viridis_d() +
    facet_wrap(~Knots, scales = "free_y")
  print(p1)
  
  
  dfnexns = df_results_spline %>%
    filter(TIME_PERIOD < ym("2022-02"), TIME_PERIOD > ym("2019-05"))
  df0exns = df_results_spline %>%
    filter(TIME_PERIOD < ym("2016-11"), TIME_PERIOD > ym("2016-06"))
  dfpexns = df_results_spline %>% 
    filter(TIME_PERIOD > ym("2016-11"),  TIME_PERIOD < ym("2019-05"))
  df0exns$period = "Zero"
  dfnexns$period = "Negative"
  dfpexns$period = "Positive"
  p2 <- rbind(dfnexns, df0exns, dfpexns) %>%
    ggplot() +
    geom_boxplot(aes(x = period, y = .data[[col_sel]], fill = period), outliers = F) +
    theme_bw() + guides(fill = "none") + labs(title = "Natural Cubic Splines") +
    facet_wrap(~Knots, scales = "free_y")
  print(p2)
  
  p3 = df_results_spline %>%
    select(all_of(col_sel), TIME_PERIOD, Type, Knots) %>%
    rem_outl(time_column = "TIME_PERIOD", value_column = col_sel) %>%
    mutate(TIME_PERIOD = as.Date(TIME_PERIOD)) %>%
    mutate(year_val = year(TIME_PERIOD)) %>%
    ggplot() +
    geom_tile(aes(x = year_val, y = Knots, fill = .data[[col_sel]])) +
    theme_bw() +
    labs(x = "", y = "Knots", title = "Natural Cubic Splines") +
    scale_fill_gradient(trans = "log10", low = "green", high = "red") +
    scale_x_continuous(breaks = 2004:2023) +
    geom_vline(xintercept = 2019.42, linetype = "dashed") +
    geom_vline(xintercept = 2022.16, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90))
  print(p3)
  
  return(list(p1, p2, p3))
}

results_train_mse <- plot_line_fn_errs("train_mse")
results_test_mse <- plot_line_fn_errs("test_mse")
results_train_mae <- plot_line_fn_errs("train_mae")
results_test_mae <- plot_line_fn_errs("test_mae")
results_r2 <- plot_line_fn_errs("r2")
results_adj_r2 <- plot_line_fn_errs("adj_r2")




