
library(fullhouse)
library(tidyverse)
library(parallel)
library(Pareto)
library(mvtnorm)
ncores = detectCores() - 1

data("batters_adjusted")
load("~/research/full_house_model/baseball_simu.RData")

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







years = 1871:2023
batters = bat_dat
cutoff = do.call(rbind, mclapply(years, mc.cores = ncores, FUN = function(xx){
  m = batters %>% filter(yearID == xx) %>% filter(PA >= 75)
  data.frame(thres = median(m$PA), yearID = xx)
}))
batters = merge(batters, cutoff, by = 'yearID')
batters = batters %>% select(playerID, yearID, lgID, name, age, lgID, AB, 
                              PA, obs_hits, H, thres, bWAR, pops) %>%
  mutate(comp = ifelse(AB != 0, H / AB, 0), full_time = ifelse(PA >= thres, 'Y', 'N'))
batters_schell = do.call(rbind, mclapply(years, mc.cores = ncores, FUN = function(xx){
  int= batters %>% filter(yearID == xx)
  lg_avg = sum(int$H)/sum(int$AB)
  int %>% mutate(comp = (H + lg_avg * 25)/(AB + 25))
}))

para_dat = batters_schell %>% filter(full_time == 'Y') %>% group_by(yearID) %>%
  summarise(mu = mean(comp), sd = sd(comp), n = n(), N = round(mean(pops)))

set.seed(13)

simu = function(n1, n2, N1, N2, n1j, n2k, mu1, mu2, sigma1, sigma2){
  simu_paradox = function(n1, n2, N1, N2, n1j, n2k, mu1, mu2, sigma1, sigma2){
    
    Y1 = sort(rnorm(n1, mu1, sigma1))
    Y2 = sort(rnorm(n2, mu2, sigma2))
    
    order_pnorm = function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4){
      p = pnorm(q, mean = mean, sd = sd)
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    order_pnorm_vec = function(q, mean = 0, sd = 1){
      q = sort(q) 
      n = length(q)
      unlist(lapply(1:n, function(j){
        order_pnorm(q[j], k = j, n = n, mean = mean, sd = sd)
      }))
    }
    
    u1 = order_pnorm_vec(Y1, mean = mu1, sd = sigma1)
    u2 = order_pnorm_vec(Y2, mean = mu2, sd = sigma2)
    
    order_Pareto_vec = function(u, t = 1, alpha = 1.16, npop = 1e4){
      n = length(u)
      if(length(npop) == 1) npop = rep(npop, n)
      unlist(lapply(1:n, function(j){
        qPareto(qbeta(u[j], j + npop[j]-n, n + 1 - j), t = t, alpha = alpha)
      }))
    }
    
    # try it out on the percentiles obtained from order_norm_vec
    X1 = order_Pareto_vec(u1, npop = N1)
    X2 = order_Pareto_vec(u2, npop = N2)
    return(c(Y1[n1-n1j], Y2[n2-n2k], X1[n1-n1j], X2[n2-n2k]))
  }
  
  res = do.call(rbind, mclapply(1:1000, function(xx){
    simu_paradox(n1 = n1, n2 = n2, N1 = N1, N2 = N2, n1j = n1j, n2k = n2k, 
                 mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2)
  }, mc.cores = 7)) 
  res = as.data.frame(res)
  colnames(res) = c('Y1', 'Y2', 'X1', 'X2')
  return(nrow(res[res$Y1 < res$Y2 & res$X1 > res$X2,]))
}

prob_paradox = function(n1, n2, N1, N2, mu1, mu2, sigma1, sigma2){
  x = c(rep((mu2-mu1)/sigma1, n1), rep(0, (n2-1)))
  mu = rep(0, (n1+n2-1))
  sigma = diag((n1+n2-1))
  sigma[1:n1, 1:n1] = sigma[1:n1, 1:n1] + (sigma2^2/sigma1^2)
  sigma[(n1+1):(n1+n2-1), 1:n1] = 
    sigma[(n1+1):(n1+n2-1), 1:n1] + (sigma2/sigma1)
  sigma[1:n1, (n1+1):(n1+n2-1)] = 
    sigma[1:n1, (n1+1):(n1+n2-1)] + (sigma2/sigma1)
  sigma[(n1+1):(n1+n2-1), (n1+1):(n1+n2-1)] = 
    sigma[(n1+1):(n1+n2-1), (n1+1):(n1+n2-1)] + 1
  
  return(n2*pmvnorm(lower = -Inf, upper = x, mean = mu, sigma = sigma) - 
           (N2/(N1+N2)))
}

para_AVG = do.call(rbind, mclapply(1:152, mc.cores = ncores, FUN = function(xx){
  # res_test = replicate(100, simu(n1 = para_dat$n[xx], 
  #                               n2 = para_dat$n[153], 
  #                               N1 = para_dat$N[xx], 
  #                               N2 = para_dat$N[153], 
  #                               n1j = 0, 
  #                               n2k = 0, 
  #                               mu1 = para_dat$mu[xx], 
  #                               mu2 = para_dat$mu[153], 
  #                               sigma1 = para_dat$sd[xx], 
  #                               sigma2 = para_dat$sd[153]))
  
  num_test = replicate(100, prob_paradox(n1 = para_dat$n[153],
                                               n2 = para_dat$n[xx],
                                               N1 = para_dat$N[153],
                                               N2 = para_dat$N[xx],
                                               mu1 = para_dat$mu[153],
                                               mu2 = para_dat$mu[xx],
                                               sigma1 = para_dat$sd[153],
                                               sigma2 = para_dat$sd[xx]))
  num_test_mean = mean(num_test)
  num_test_se = sd(num_test) / 10
  
  # data.frame(yearID = para_dat$yearID[xx], 
  #            simu_prob = mean(res_test)/1000, 
  #            norm_prob = ifelse(num_test < 0, 0, round(num_test, 4)))
  data.frame(yearID = para_dat$yearID[xx], 
             norm_prob = ifelse(num_test_mean < 0, 0, round(num_test_mean, 4)),
             num_test_se = num_test_se)
}))

qknitr::kable(para_AVG)



ggplot(data = para_AVG, aes(x = yearID, y = norm_prob)) + geom_line() +
  geom_smooth()


