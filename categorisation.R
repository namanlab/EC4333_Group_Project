library(tidyverse)
library(minpack.lm)  # For non-linear least squares fitting
library(forecast)
library(YieldCurve)
library(gridExtra)

df_pr_aaa = read.csv("data/df_spot_rates_aaa.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD))  %>%
  rename("AAA" = OBS_VALUE)
df_pr_all = read.csv("data/df_spot_rates_all.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD)) %>%
  rename("All" = OBS_VALUE)
df_pr = df_pr_aaa %>% inner_join(df_pr_all) %>%
  pivot_longer(3:4, names_to = "Type", values_to = "OBS_VALUE")

df_pr %>% 
  group_by(TIME_PERIOD, Type) %>%
  summarise(val = mean(OBS_VALUE, na.rm = T)) %>%
  ungroup() %>%
  ggplot() +
  annotate("rect", xmin = ym("2019-05"), xmax = ym("2022-02"), ymin = -0.9, ymax = 5,
           fill = "red", alpha = 0.6) +
  annotate("rect", xmin = ym("2016-06"), xmax = ym("2016-11"), ymin = -0.9, ymax = 5,
           fill = "red", alpha = 0.6) +
  geom_line(aes(x = TIME_PERIOD, y = val, color = Type)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Average Yield") +
  scale_color_viridis_d()

# Time Periods: 
# Close to 0: 2016-06 to 2016-11
# Negative: 2019-05 to 2022-02
# Positive: all other time periods

dfn = df_pr %>%
  filter(TIME_PERIOD < ym("2022-02"),
         TIME_PERIOD > ym("2019-05"))

df0 = df_pr %>%
  filter(TIME_PERIOD < ym("2016-11"),
         TIME_PERIOD > ym("2016-06"))

dfp = df_pr %>%
  filter(TIME_PERIOD > ym("2016-11"),
         TIME_PERIOD < ym("2019-05"))

