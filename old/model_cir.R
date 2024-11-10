library(tidyverse)

df <- read.csv("data/df_spot_rates_aaa.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD))

df_date <- df %>% filter(TIME_PERIOD == "2011-09-09")
obs_value <- df_date$OBS_VALUE

# OBTAINING CIR PARAMETER ESTIMATES FOR CIR USING THE OBSERVED VALUES
ols_cir <- function(obs_value, dt) {
  n_steps <- length(obs_value)
  rs <- obs_value[1:n_steps - 1]
  rt <- obs_value[2:n_steps]
  
  y <- (rt - rs) / sqrt(rs)
  z1 <- dt / sqrt(rs)
  z2 <- dt * sqrt(rs)
  data <- data.frame(y, z1, z2)

  model <- lm(y ~ 0 + z1 + z2, data = data)
  resid <- model$residuals
  
  ssresidual <- (t(resid) %*% resid)[1]
  sstotal <- (t(y) %*% y)[1] - (sum(y)^2) / length(y)
  adj_r_squared <-  1 - ssresidual/(length(y)-3) / sstotal/(length(y)-1)
  
  coefficients <- model$coefficients
  k <- -coefficients["z2"]
  theta <- coefficients["z1"] / k
  sigma <- sqrt(var(resid) / dt)
  cir_params <- c(k, theta, sigma) # alpha, mu
  
  return(list(y, model$fitted.values, cir_params, adj_r_squared))
}

model <- ols_cir(obs_value, 1/length(obs_value))
y <- model[[1]]
fitted <- model[[2]]
cir_params <- model[[3]]
adj_r_squared <- model[[4]]

df_actual_fitted <- data.frame(tail(df_date$MAT, -1), y, fitted)
colnames(df_actual_fitted) <- c("maturity", "y", "fitted")

## PLOTTING THE ACTUAL AND FITTED VALUES FOR THE OLS MODEL
ggplot(df_actual_fitted) +
  geom_line(mapping = aes(x=maturity, y=y), color="red") +
  geom_line(mapping = aes(x=maturity, y=fitted), color="blue") +
  annotate("text", x=280, y=0.06, label=sprintf("adj_r_squared: %s", adj_r_squared))


# SIMULATING THE YIELD CURVE USING THE FITTED CIR PARAMETERS FROM ABOVE
sim_cir <- function(k, theta, sigma, r0, N) {
  dt <- 1 / N
  yield_paths <- c(r0)
  
  for (t in 2:N) {
    z <- rnorm(1) # generate 1 sample from standard normal
    r <- yield_paths[[t-1]]
    print(k * (theta - r) * dt)
    yield <- r + k * (theta - r) * dt + sigma * sqrt(dt) * sqrt(max(0.01, r)) * z
    yield_paths <- append(yield_paths, yield)
  }
  
  return(yield_paths)
}

yield_paths <- sim_cir(cir_params[1], cir_params[2], cir_params[3], obs_value[1], length(obs_value))
yield_paths_df <- data.frame(df_date$MAT, yield_paths)

## PLOTTING THE SIMULATED YIELD CURVE USING THE FITTED COEFFICIENTS
ggplot() + 
  geom_line(data=df_date, aes(x=MAT, y=OBS_VALUE), color="red") +
  geom_line(data=yield_paths_df, aes(x=df_date.MAT, y=yield_paths), color="blue")

