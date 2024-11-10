library(tidyverse)
library(minpack.lm)  # For non-linear least squares fitting
library(forecast)
library(YieldCurve)
library(gridExtra)

df_pr = read.csv("df_japan_spot_rates.csv") %>%
  mutate(TIME_PERIOD = ymd(TIME_PERIOD)) 

df_pr %>% 
  group_by(TIME_PERIOD) %>%
  summarise(val = mean(OBS_VALUE, na.rm = T)) %>%
  ungroup() %>%
  ggplot() +
  annotate("rect", xmin = ym("2019-05"), xmax = ym("2020-02"), ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.6) +
  annotate("rect", xmin = ym("2016-03"), xmax = ym("2016-10"), ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.6) +
  geom_line(aes(x = TIME_PERIOD, y = val)) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Average Yield") +
  scale_color_viridis_d()


dfn1 = df_pr %>%
  filter(TIME_PERIOD < ym("2020-02"),
         TIME_PERIOD > ym("2019-05"))
dfn2 = df_pr %>%
  filter(TIME_PERIOD < ym("2016-10"),
         TIME_PERIOD > ym("2016-03"))
dfn = rbind(dfn1, dfn2)

dfp = df_pr %>%
  filter(TIME_PERIOD > ym("2016-11"),
         TIME_PERIOD < ym("2019-05"))

