
library(tidyverse)
library(fullhouse)
library(parallel)
library(Pareto)
library(mvtnorm)
library(Lahman)

paradox_prob_results_2022 = read_csv("paradox_prob_results_2022.csv")
paradox_prob_results_2023 = read_csv("paradox_prob_results_2023.csv")
bat_season_para = read_csv("bat_season_para.csv")
pops_full_time = read_csv("BA_normal.csv")

## compute n and N for all years
sample_sizes = pops_full_time %>% group_by(yearID) %>% 
  summarise(n = sum(full_time == "Y"), 
            pops = max(pops), 
            thres = max(thres))

## compute n and N for all years
pops_full_time = read_csv("BA_normal.csv")
sample_sizes = pops_full_time %>% group_by(yearID) %>% 
  summarise(n = sum(full_time == "Y"), 
            pops = max(pops), 
            thres = max(thres))

## get estimated seasonal means and sds for batting average 
## for full time players
bat_season_para = read_csv("bat_season_para.csv")
foo = pops_full_time %>% 
  rename(year = yearID) %>% 
  select(playerID, year, thres)
mean_sd = bat_season_para %>% 
  left_join(foo) %>% 
  filter(AB >= thres) %>% 
  summarise(mean_BA = mean(BA), sd_BA = sd(BA), .by = year) %>% 
  arrange(year)

## Normal distribution means and standard deviations 
## over time
batting_averages_raw_full_time = Batting %>% 
  rename(year = yearID) %>% 
  left_join(sample_sizes %>% 
              rename(year = yearID) %>% 
              select(year, thres)) %>% 
  filter(AB >= thres) %>% 
  mutate(BA = H/AB) 

batting_averages_raw_full_time %>% 
  summarise(Mean = mean(BA), SD = sd(BA), .by = year) %>% 
  rename("Standard Deviation" = SD) %>% 
  pivot_longer(cols = Mean:`Standard Deviation`, names_to = "quantity") %>% 
  ggplot() + 
  aes(x = year, y = value) + 
  #geom_line() + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 20)) +
  facet_wrap(.~ quantity, scales = "free_y") + 
  labs(title = "Batting average means and standard deviations", 
       subtitle = "Full-time players from 1871-2023") + 
  theme_minimal()

## t-test difference in 2022 and 2023 batting average means
foo = batting_averages_raw_full_time %>% 
  filter(year %in% c(2022,2023)) 
bar = split(foo, f = as.factor(foo$year))
baz = lapply(bar, function(xx) xx$BA)
t.test(baz[[1]], baz[[2]])
### 2022
shapiro.test(baz[[1]])
### 2023
shapiro.test(baz[[2]])

## hits per game
Teams %>% 
  filter(yearID %in% c(2022,2023)) %>% 
  summarise(sum_H = sum(H), .by = yearID) %>% 
  mutate(HpG = sum_H / (162*30)) %>% 
  pull(HpG) %>% 
  diff()

## BABIP difference
Teams %>% 
  filter(yearID %in% c(2022,2023)) %>% 
  summarise(sum_H = sum(H), 
            sum_AB = sum(AB),
            sum_HR = sum(HR), 
            sum_SO = sum(SO), 
            sum_SF = sum(SF), .by = yearID) %>% 
  mutate(BABIP = (sum_H - sum_HR)/(sum_AB - sum_SO - sum_HR + sum_SF)) %>% 
  pull(BABIP) %>% 
  diff()



## Extreme right tail of batting average distributions over time
batting_averages_raw_full_time %>% 
  summarise(mean = mean(BA), sd = sd(BA), .by = year) %>% 
  mutate(three_sigma = mean + 3*sd) %>% 
  ggplot() + 
  aes(x = year, y = three_sigma) + 
  geom_line() + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 20)) +
  labs(title = "Three standard deviations from mean batting averages", 
       subtitle = "1871-2023") + 
  ylab("") + 
  theme_minimal()



## Paradox Probability plot
#paradox_prob_results = bind_rows(paradox_prob_results_2022, 
#                                 paradox_prob_results_2023, 
#                                 .id = "year_target") %>% 
#  mutate(year_target = ifelse(year_target == 1, 2022, 2023))
ggplot(data = paradox_prob_results_2023, aes(x = year, y = norm_prob)) + 
  geom_line() +
  geom_smooth() + 
  ylim(0,1) + 
  labs(title = "Normal Paradox Probabilities for Max Batting Averages", 
       subtitle = "Comparing year 2023 season to all others") + 
  xlab("Year") + 
  ylab("Paradox Probability") +
  theme_minimal()


## Era-adjusted batting leaders comparison assuming 
## batting averages are normally distributed
qux = bat_season_para %>% 
  left_join(sample_sizes %>% 
              rename(year = yearID) %>% 
              select(year, thres)) %>% 
  filter(AB >= thres) %>% 
  summarise(max_BA = max(BA), .by = year) %>% 
  arrange(year)
Arraez_2023 = qux %>% filter(year == 2023) %>% pull(max_BA)
McNeil_2022 = qux %>% filter(year == 2022) %>% pull(max_BA)

qux %>% 
  mutate(era_paradox = as.numeric(max_BA > Arraez_2023)) %>% 
  ggplot() + 
  aes(x = year, y = era_paradox) + 
  geom_point() +
  #geom_smooth() + 
  geom_smooth(method = "glm", formula = y ~ poly(x, 3)) +   
  ylim(0,1) + 
  labs(title = "Era-adjusted Max Batting Averages Exceeding 2023 Max") + 
  xlab("Year") + 
  ylab("") +
  theme_minimal()  



## Batting leaders comparison 
foobar = batting_averages_raw_full_time %>% 
  summarise(max_BA = max(BA), .by = year) %>% 
  arrange(year)
Arraez_2023_raw = foobar %>% filter(year == 2023) %>% pull(max_BA)
McNeil_2022_raw = foobar %>% filter(year == 2022) %>% pull(max_BA)

mu2022 = mean(baz[[1]]); sd2022 = sd(baz[[1]]); n2022 = length(baz[[1]])
pnorm(McNeil_2022_raw, mean = mu2022, sd = sd2022, lower = TRUE)^n2022

mu2023 = mean(baz[[1]]); sd2023 = sd(baz[[1]]); n2023 = length(baz[[1]])
pnorm(Arraez_2023_raw, mean = mu2023, sd = sd2023, lower = TRUE)^n2023

## Max order stat percentiles
### and shapiro wilk test results
extreme_prob = batting_averages_raw_full_time %>% 
  summarise(max_BA = max(BA), 
            n = n(), 
            mean_BA = mean(BA), 
            sd_BA = sd(BA),
            shap_test = shapiro.test(BA)$p.value,
            .by = year) %>% 
  mutate(extreme_prob = pnorm(max_BA, mean = mean_BA, sd = sd_BA)^n) %>% 
  arrange(year)

more_extreme_Arraez = 
  extreme_prob[extreme_prob["extreme_prob"] >= 
               extreme_prob[extreme_prob$year == 2023, "extreme_prob"],  ] %>% 
  select(year, max_BA, n, shap_test, extreme_prob)
nrow(more_extreme_Arraez)
nrow(more_extreme_Arraez[more_extreme_Arraez$year > 1947, ])

more_extreme_Arraez_names = more_extreme_Arraez %>% left_join(
batting_averages_raw_full_time %>% 
  select(playerID, year, BA) %>% 
  left_join(People %>% 
              mutate(name = paste(nameFirst, nameLast, sep = " ")) %>% 
              select(playerID, name)) %>% 
  select(-playerID) %>%
  group_by(year) %>%
  filter(BA == max(BA, na.rm = TRUE)) %>%
  ungroup() %>%
  select(year, name, BA) %>% 
  rename(max_BA = BA)) %>% 
  select(year, name, max_BA, n, shap_test, extreme_prob) %>% 
  mutate(extreme_prob = round(extreme_prob, 4), 
         shap_test = round(shap_test, 4))

library(xtable)  
xtable(more_extreme_Arraez_names, 
       digits = c(0, 0, 0, 3, 0, 4, 4)) %>% 
  print(include.rownames = FALSE)


### Uniform Q-Q plot
frac_seq = (1:nrow(extreme_prob) / nrow(extreme_prob))
plot(frac_seq,  sort(extreme_prob$extreme_prob), 
     main = "Uniform Q-Q Plot for Extreme Batting Averages", 
     xlab = "Theoretical Quantiles", 
     ylab = "Observed Quantiles")
lines(frac_seq, frac_seq)



foobar %>% 
  mutate(era_paradox = as.numeric(max_BA > Arraez_2023_raw)) %>% 
  ggplot() + 
  aes(x = year, y = era_paradox) + 
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ poly(x, 3)) + 
  ylim(0,1) + 
  labs(title = "Max Batting Averages Exceeding 2023 Max") + 
  xlab("Year") + 
  ylab("") +
  theme_minimal()  

blah = qux %>% 
  mutate(era_paradox = as.numeric(max_BA > Arraez_2023_raw)) 

lapply(1:10, FUN = function(j){ 
  AIC(glm(era_paradox ~ poly(year, j), data = blah, family = "binomial"))
})  

lapply(1:10, FUN = function(j){ 
  AIC(glm(era_paradox ~ poly(year, j), data = baz %>% 
            mutate(era_paradox = as.numeric(max_BA > Arraez_2023), 
                   family = "binomial")))
})  


## compare paradox probabilities to  
## probability of max BA exceeding 2023 Luis Arraez
foobar %>% 
  mutate(era_paradox = as.numeric(max_BA > Arraez_2023_raw)) %>% 
  ggplot(aes(x = year)) + 
  # 6th order polynomial from AIC (raw BA)
  geom_smooth( # formula = y ~ poly(x, 6)
    aes(y = era_paradox, color = "Smoothed Prob (Raw BA)", linetype = "Smoothed Prob (Raw BA)"),
    se = FALSE
  ) + 
  # Paradox probability line
  geom_line(
    data = paradox_prob_results_2023, 
    aes(y = norm_prob, color = "Paradox Prob", linetype = "Paradox Prob"), 
    lwd = 1
  ) +
  ylim(0,1) + 
  labs(
    title = "Comparison of Paradox Probabilities and Results", 
    subtitle = "Comparing 2023 season to all others (1871-2022)",
    x = "Year",
    y = "",
    color = "Legend",
    linetype = "Legend"
  ) +
  scale_color_manual(values = c("Smoothed Prob (Raw BA)" = "black", "Paradox Prob" = "blue")) +
  scale_linetype_manual(values = c("Smoothed Prob (Raw BA)" = "dashed", "Paradox Prob" = "solid")) +
  theme_minimal()


foobar %>% 
  mutate(era_paradox = as.numeric(max_BA > McNeil_2022_raw)) %>% 
  ggplot() + 
  aes(x = year, y = era_paradox) + 
  # 6th order polynomial from AIC (raw BA)
  geom_smooth(#method = "glm", formula = y ~ poly(x, 6), 
    color = "black", lty = 2) + 
  geom_line(data = paradox_prob_results_2022, 
            aes(y = norm_prob), color = "blue", lwd = 1, lty = 1) + 
  # 5th order polynomial from AIC (era-adjusted BA)
  geom_smooth(data = qux %>% 
                mutate(era_paradox = as.numeric(max_BA > McNeil_2022)), 
              #method = "glm", formula = y ~ poly(x, 5), 
              aes(y = era_paradox), color = "red", lty = 2) +
  ylim(-0.2,1.2) + 
  labs(title = "Comparison of Paradox Probabilities to Results in 2022") + 
  xlab("Year") + 
  ylab("") +
  theme_minimal()  

