
library(fullhouse)
library(tidyverse)
data("batters_adjusted")
load("research/full_house_model/baseball_simu.RData")

foo = batters_adjusted %>%
  filter(AB >= 320) %>% 
  group_by(year) %>% 
  summarise(maxBA = max(H/AB)) %>% 
  arrange(desc(maxBA)) %>% 
  as.data.frame()

bar = foo %>% 
  mutate(resp = ifelse(maxBA <= 0.3704348, 1,0)) %>% 
  filter(year < 2023)

#m1 = glm(resp ~ poly(year,1), data = bar, family = "binomial")
#summary(m1)
#baz = cbind(bar, preds = predict(m1, type = "response"))
#baz

ggplot(bar) + 
  aes(x = year, y = resp) + 
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              formula = y ~ poly(x,1),
              se = FALSE) + 
  geom_point() +
  theme_minimal()

bat_raw = read_csv("~/website/era-adjustment-website/src/data/raw_batter.csv")
#View(bat_dat %>% filter(name == "Jeff McNeil"))
#View(bat_raw %>% filter(name == "Jeff McNeil"))

#bat_raw = bat_dat
foo_raw = bat_raw %>%
  filter(AB >= 320) %>% 
  group_by(yearID) %>% 
  summarise(maxBA = max(H/AB)) %>% 
  arrange(desc(maxBA)) %>% 
  as.data.frame()

maxBA2023 = foo_raw %>% filter(yearID == 2023) %>% pull(maxBA)




bar_raw = foo_raw %>% 
  mutate(resp = ifelse(maxBA >= maxBA2023, 1,0)) %>% 
  filter(yearID < 2023)

#m1_raw = glm(resp ~ poly(yearID,2), data = bar_raw, family = "binomial")
#summary(m1_raw)
#baz_raw = cbind(bar_raw, preds = predict(m1_raw, type = "response"))
#baz_raw

ggplot(bar_raw) + 
  aes(x = yearID, y = resp) + 
  geom_smooth(method = "glm", 
              formula = y ~ poly(x,2),
              method.args = list(family = "binomial"), 
              se = TRUE) + 
  geom_point() +
  theme_minimal()



