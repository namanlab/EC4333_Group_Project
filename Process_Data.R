library(tidyverse)
library(minpack.lm)  # For non-linear least squares fitting

# All years (AAA bonds data) from https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/index.en.html
df = read_csv("data/data_aaa.csv")

df = df %>% dplyr::select(DATA_TYPE_FM, TIME_PERIOD, OBS_VALUE) %>% na.omit()
df_pr <- df %>% filter(str_detect(DATA_TYPE_FM, "SR_")) %>%
  mutate(DATA_TYPE_FM = str_split_fixed(DATA_TYPE_FM, "_", 2)[,2]) %>%
  dplyr::select(TIME_PERIOD, MAT = DATA_TYPE_FM, OBS_VALUE)
C = unique(df_pr$MAT)

convert_to_months <- function(value) {
  # Handle the "Y_Y" format
  if (grepl("Y_\\d+Y", value)) {
    # Split by the underscore, sum the years, and convert to months
    parts <- strsplit(value, "_")[[1]]
    total_years <- sum(as.numeric(gsub("Y", "", parts)))
    return(total_years * 12)
  }
  # Handle the "Y_M" format
  if (grepl("Y_\\d+M", value)) {
    # Split by the underscore, sum the years, and convert to months
    parts <- strsplit(value, "_")[[1]]
    total_years <- sum(as.numeric(gsub("Y", "", parts[1])))
    total_mnths <- sum(as.numeric(gsub("M", "", parts[2])))
    return(total_years * 12 + total_mnths)
  }
  # Handle the "YM" format
  if (grepl("Y\\d+M", value)) {
    years_months <- as.numeric(gsub("Y|M", "", unlist(strsplit(value, "Y|M"))))
    return(years_months[1] * 12 + years_months[2])
  }
  # Handle the "Y" format (years only)
  if (grepl("Y$", value)) {
    years <- as.numeric(gsub("Y", "", value))
    return(years * 12)
  }
  # Handle the "M" format (months only)
  if (grepl("M$", value)) {
    months <- as.numeric(gsub("M", "", value))
    return(months)
  }
  # Return NA if format doesn't match
  return(NA)
}

# Apply the conversion to the MAT column in your data frame
df_pr <- df_pr %>%
  mutate(MAT = sapply(MAT, convert_to_months))
df_pr <- na.omit(df_pr) %>% arrange(TIME_PERIOD, MAT)

write.csv(df_pr, "data/df_spot_rates_aaa.csv", row.names = F)



################################################################################
################################### ALL BONDS ##################################
################################################################################

library(tidyverse)
library(minpack.lm)  # For non-linear least squares fitting

# All years (All bonds data) from https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/index.en.html
df = read_csv("data/data_all.csv")

df = df %>% dplyr::select(DATA_TYPE_FM, TIME_PERIOD, OBS_VALUE) %>% na.omit()
df_pr2 <- df %>% filter(str_detect(DATA_TYPE_FM, "SR_")) %>%
  mutate(DATA_TYPE_FM = str_split_fixed(DATA_TYPE_FM, "_", 2)[,2]) %>%
  dplyr::select(TIME_PERIOD, MAT = DATA_TYPE_FM, OBS_VALUE)
C = unique(df_pr2$MAT)

convert_to_months <- function(value) {
  # Handle the "Y_Y" format
  if (grepl("Y_\\d+Y", value)) {
    # Split by the underscore, sum the years, and convert to months
    parts <- strsplit(value, "_")[[1]]
    total_years <- sum(as.numeric(gsub("Y", "", parts)))
    return(total_years * 12)
  }
  # Handle the "Y_M" format
  if (grepl("Y_\\d+M", value)) {
    # Split by the underscore, sum the years, and convert to months
    parts <- strsplit(value, "_")[[1]]
    total_years <- sum(as.numeric(gsub("Y", "", parts[1])))
    total_mnths <- sum(as.numeric(gsub("M", "", parts[2])))
    return(total_years * 12 + total_mnths)
  }
  # Handle the "YM" format
  if (grepl("Y\\d+M", value)) {
    years_months <- as.numeric(gsub("Y|M", "", unlist(strsplit(value, "Y|M"))))
    return(years_months[1] * 12 + years_months[2])
  }
  # Handle the "Y" format (years only)
  if (grepl("Y$", value)) {
    years <- as.numeric(gsub("Y", "", value))
    return(years * 12)
  }
  # Handle the "M" format (months only)
  if (grepl("M$", value)) {
    months <- as.numeric(gsub("M", "", value))
    return(months)
  }
  # Return NA if format doesn't match
  return(NA)
}

# Apply the conversion to the MAT column in your data frame
df_pr2 <- df_pr2 %>%
  mutate(MAT = sapply(MAT, convert_to_months))
df_pr2 <- na.omit(df_pr2) %>% arrange(TIME_PERIOD, MAT)

write.csv(df_pr2, "data/df_spot_rates_all.csv", row.names = F)







