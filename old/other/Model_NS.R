library(tidyverse)
library(minpack.lm)  # For non-linear least squares fitting

df_pr = read.csv("df_spot_rates.csv")


# MODEL

f2 <- function(beta1,theta,tau){
  beta1*(1-exp(-theta/tau))/(theta/tau)
}

f3 <- function(beta2,theta,tau){
  beta2*((1-exp(-theta/tau))/(theta/tau) - exp(-theta/tau))
}

f4 <- function(beta3,theta,tau2){
  beta3*((1-exp(-theta/tau2))/(theta/tau2) - exp(-theta/tau2))
}

ns <- function(theta, beta0,beta1,beta2,tau){
  beta0 + f2(beta1,theta,tau) + f3(beta2,theta,tau)
}

extendns <- function(theta, beta0,beta1,beta2,beta3,tau,tau2){
  beta0 + f2(beta1,theta,tau) + f3(beta2,theta,tau) + f4(beta3,theta,tau2) 
}

ns_solve_fn = function(x, params){
  extendns(x, params[1], params[2],params[3],params[4],params[5],params[6])
}

str(df_pr)



# Filter data for the specified date
date_val <- "2019-05-02"
data_rel <- df_pr %>% filter(TIME_PERIOD == date_val)

# Define the residual function for optimization
residuals_ns <- function(params, MAT, VAL) {
  res <- VAL - sapply(MAT, ns_solve_fn, params = params)
  return(res)
}

# Initial guesses for the parameters (can be adjusted)
init_params <- c(1, -1, 1, 1, 1, 1)

# Perform non-linear least squares fitting
fit <- nls.lm(par = init_params, 
              fn = function(x) residuals_ns(x, MAT = data_rel$MAT, VAL = data_rel$OBS_VALUE))

# Extract the fitted parameters
params_est <- fit$par
print(params_est)

# Predict values using the fitted parameters
data_rel <- data_rel %>%
  mutate(
    Predicted_Yield = sapply(MAT, ns_solve_fn, params = params_est),
    Residual = OBS_VALUE - Predicted_Yield
  )

# Identify outliers (absolute residual > 2 * standard deviation of residuals)
residual_threshold <- 2 * sd(data_rel$Residual)
data_clean <- data_rel %>% filter(abs(Residual) <= residual_threshold)

dplyr::setdiff(data_rel, data_clean)

# Calculate RSS for the non-outlier data
rss <- sum((data_clean$OBS_VALUE - data_clean$Predicted_Yield)^2)

# Pivot non-outlier data to long format for ggplot
data_long <- data_clean %>%
  pivot_longer(cols = c(OBS_VALUE, Predicted_Yield), 
               names_to = "Type", 
               values_to = "Value")

# Create the plot using ggplot2
ggplot(data_long, aes(x = MAT, y = Value, color = Type, linetype = Type)) +
  geom_point(data = filter(data_long, Type == "OBS_VALUE"), size = 1) +
  geom_line(size = 1) +
  labs(
    title = paste("Fitted vs Observed Values on", date_val),
    subtitle = paste("Residual Sum of Squares (RSS):", round(rss, 4)),
    x = "Maturity (MAT)",
    y = "Spot Rates",
    color = "Legend",  # Customize legend title
    linetype = "Legend"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  scale_color_manual(values = c("OBS_VALUE" = "blue", "Predicted_Yield" = "red"))




