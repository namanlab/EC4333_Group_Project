library(tidyverse)
library(minpack.lm)
library(forecast)
library(YieldCurve)
library(gridExtra)


df1 = read.csv("data/df_results_exns.csv")
df1$Model = "NSS"
df2 = read.csv("data/df_results_ns.csv")
df2$Model = "NS"
df3 = read.csv("data/df_results_spline.csv") %>% filter(Knots == 3) %>%
  select(-Knots)
df3$Model = "3-Spline"

df <- rbind(df1, df2, df3) %>%
  filter(train_mse <= 0.1, train_mae <= 0.1, test_mse <= 0.1, test_mae <= 0.1)


colnames(df)

df %>%
  mutate(TIME_PERIOD = as.Date(TIME_PERIOD)) %>%
  ggplot() +
  annotate("rect", xmin = ym("2019-05"), xmax = ym("2022-02"), ymin = 0, ymax = Inf,
           fill = "red", alpha = 0.6) +
  annotate("rect", xmin = ym("2016-06"), xmax = ym("2016-11"), ymin = 0, ymax = Inf,
           fill = "red", alpha = 0.6) +
  geom_area(aes(x = TIME_PERIOD, y = test_mse, fill = Model)) +
  theme_bw() +
  labs(x = "", title = "Extended NSS") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d() +
  facet_wrap(~Type)

df %>%
  mutate(TIME_PERIOD = as.Date(TIME_PERIOD)) %>%
  ggplot() +
  annotate("rect", xmin = ym("2019-05"), xmax = ym("2022-02"), ymin = 0, ymax = Inf,
           fill = "red", alpha = 0.4) +
  annotate("rect", xmin = ym("2016-06"), xmax = ym("2016-11"), ymin = 0, ymax = Inf,
           fill = "red", alpha = 0.4) +
  geom_area(aes(x = TIME_PERIOD, y = test_mse, fill = Model), position = "fill") +
  geom_vline(xintercept = ym("2019-05"), linetype = "dashed", color = "white", size = 1) +
  geom_vline(xintercept = ym("2022-02"), linetype = "dashed", color = "white", size = 1) +
  theme_bw() +
  labs(x = "", title = "Comparison of Model Errors") +
  theme(legend.position = "bottom") +
  scale_fill_viridis_d() +
  facet_wrap(~Type)


df %>%
  mutate(TIME_PERIOD = as.Date(TIME_PERIOD)) %>%
  ggplot() +
  annotate("rect", xmin = ym("2019-05"), xmax = ym("2022-02"), ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.6) +
  annotate("rect", xmin = ym("2016-06"), xmax = ym("2016-11"), ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.6) +
  annotate("rect", xmin = ym("2006-03"), xmax = ym("2008-03"), ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = 0.3) +
  annotate("rect", xmin = ym("2021-08"), xmax = ym("2023-11"), ymin = -Inf, ymax = Inf,
           fill = "blue", alpha = 0.3) +
  geom_line(aes(x = TIME_PERIOD, y = adj_r2, color = Model)) +
  theme_bw() +
  labs(x = "", title = "Extended NSS") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d() +
  facet_wrap(~Type)

ggplotly(p)


# qualitative analysis




df_pr_aaa = read.csv("data/df_spot_rates_aaa.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD)) %>%
  rename("AAA" = OBS_VALUE)
df_pr_all = read.csv("data/df_spot_rates_all.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD)) %>%
  rename("All" = OBS_VALUE)
df_pr = df_pr_aaa %>% inner_join(df_pr_all) %>%
  pivot_longer(3:4, names_to = "Type", values_to = "OBS_VALUE")


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
                   values_to = "Value") %>%
      mutate(Dataset = ifelse(Dataset == 1, "Train", "Test"))
    
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
      facet_wrap(~ Dataset, ncol = 2, scales = "free_y") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom"
      ) +
      scale_color_manual(values = c("OBS_VALUE" = "blue", "Predicted_Yield" = "red"))
    print(p1)
  }
  
  return(lst_vals)
}

"2006-11-09"
"2007-03-20"
"2023-07-07"
# abnormal_eg_1.png

fit_model_fn(date_val = "2023-07-07", viz = T, type = "All", num_knots = 3,
             outlier_handle = F)
