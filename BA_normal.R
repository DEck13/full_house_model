rm(list=ls())
library(tidyverse)
library(orderstats)
library(Pareto)
library(doParallel)
library(splines)
library(retrosheet)
library(kableExtra)
library(Lahman)
ncores <- detectCores() - 1
years <- 1871:2023
setwd("~/Desktop")
load("pop_data.RData")
load("bat_dat.RData")
load("pit_dat.RData")
source('source.R')

batters <- bat_dat

cutoff <- do.call(rbind, mclapply(years, mc.cores = ncores, FUN = function(xx){
  m <- batters %>% filter(yearID == xx) %>% filter(PA >= 75)
  data.frame(thres = median(m$PA), yearID = xx)
}))

batters <- merge(batters, cutoff, by = 'yearID')

batters <- batters %>% select(playerID, yearID, lgID, name, age, lgID, AB, 
                              PA, obs_hits, H, thres, bWAR, pops) %>%
  mutate(comp = ifelse(AB != 0, H / AB, 0), full_time = ifelse(PA >= thres, 'Y', 'N'))

batters_schell <- do.call(rbind, mclapply(years, mc.cores = ncores, FUN = function(xx){
  int<- batters %>% filter(yearID == xx)
  lg_avg <- sum(int$H)/sum(int$AB)
  int %>% mutate(comp = (H + lg_avg * 25)/(AB + 25))
}))

batters_talent_AVG <- do.call(rbind, mclapply(years, mc.cores = ncores, function(yy){
  talent_computing_para(dataset = batters_schell, year = yy, alpha = 1.16)
})) %>% arrange(-talent)

no_strike <- batters_talent_AVG %>%
  filter(yearID < 1990, yearID >= 1977, yearID != 1981,
         full_time == 'Y', lgID == 'NL') %>%
  arrange(talent)

min_refAVG <- min(no_strike$comp)

t <- isoreg(no_strike$talent, no_strike$comp)

talent_new <- quantile(no_strike$talent, probs = seq(0,250)/250)
comp_new <- as.stepfun(t)(talent_new)
pop_new <- round(mean(no_strike$pops))

mean_new <- mean(comp_new)
sd_new <- sd(comp_new)

career_kAB_AVG_1st <- do.call(rbind, mclapply(1:30000, function(zz){
  int <- career_talent_para(dataset = batters_talent_AVG, 
                            snippet = batters_talent_AVG[zz,], alpha = 1.16)
  int
}, mc.cores = ncores))

career_kAB_AVG_2nd <- do.call(rbind, mclapply(30001:60000, function(zz){
  int <- career_talent_para(dataset = batters_talent_AVG, 
                            snippet = batters_talent_AVG[zz,], alpha = 1.16)
  int
}, mc.cores = ncores))

career_kAB_AVG_3rd <- do.call(rbind, mclapply(60001:nrow(batters_talent_AVG), function(zz){
  int <- career_talent_para(dataset = batters_talent_AVG, 
                            snippet = batters_talent_AVG[zz,], alpha = 1.16)
  int
}, mc.cores = ncores))


career_kAB <- rbind(career_kAB_AVG_1st, 
                    career_kAB_AVG_2nd, 
                    career_kAB_AVG_3rd)

mapped_batters_AVG_para <- career_kAB %>% 
  mutate(adj_AVG = ifelse(adj_comp < min_refAVG, min_refAVG, adj_comp))
