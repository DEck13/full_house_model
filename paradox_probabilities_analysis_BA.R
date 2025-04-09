
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
  summarise(mean = mean(BA), sd = sd(BA), .by = year) %>% 
  pivot_longer(cols = mean:sd, names_to = "quantity") %>% 
  ggplot() + 
  aes(x = year, y = value) + 
  #geom_line() + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 20)) +
  facet_wrap(.~ quantity, scales = "free_y") + 
  labs(title = "Batting average means and standard deviations", 
       subtitle = "1871-2023") +
  theme_minimal()

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
  labs(title = "Normal Paradox Probabilities Over Time", 
       subtitle = "Comparing year 2023 season to all others") + 
  xlab("Year") + 
  ylab("Paradox Probability") +
  theme_minimal()


## Era-adjusted batting leaders comparison assuming 
## batting averages are normally distributed
baz = bat_season_para %>% 
  left_join(sample_sizes %>% 
              rename(year = yearID) %>% 
              select(year, thres)) %>% 
  filter(AB >= thres) %>% 
  summarise(max_BA = max(BA), .by = year) %>% 
  arrange(year)
Arraez_2023 = baz %>% filter(year == 2023) %>% pull(max_BA)

baz %>% 
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
qux = batting_averages_raw_full_time %>% 
  summarise(max_BA = max(BA), .by = year) %>% 
  arrange(year)
Arraez_2023_raw = qux %>% filter(year == 2023) %>% pull(max_BA)

qux %>% 
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


