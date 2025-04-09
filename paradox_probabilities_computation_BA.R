library(tidyverse)
library(fullhouse)
library(parallel)
library(Pareto)
library(mvtnorm)
library(Lahman)

# set seed for parallel job
set.seed(13)
RNGkind(kind = "L'Ecuyer-CMRG")

# load in annual information about full-time 
# player status
pops_full_time = read_csv("BA_normal.csv")

# compute n and N for all years
sample_sizes = pops_full_time %>% group_by(yearID) %>% 
  summarise(n = sum(full_time == "Y"), 
            pops = max(pops), 
            thres = max(thres))


# get estimated seasonal means and sds for batting average 
# for full-time players
## calculate batting average for all players
BA_year_player = Batting %>% 
  rename(year = yearID) %>% 
  summarise(BA = H/AB, AB, .by = c(playerID, year))

## get annual thresholds for full-time players
foo = pops_full_time %>% 
  rename(year = yearID) %>% 
  select(playerID, year, thres)

## restrict attention to BAs for full-time players
mean_sd = BA_year_player %>% 
  left_join(foo) %>% 
  filter(AB >= thres) %>% 
  summarise(mean_BA = mean(BA), sd_BA = sd(BA), .by = year) %>% 
  arrange(year)

# paradox probability calculator function for normal distributions
# n2 P_W(w) - N2/(N1 + N2)
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


# compute paradox probabilities for all baseball seasons 
# with respect to a target season
## Monte Carlo aggregates across B reps 
## Fitting algorithm implemented using methodology is described in 
## Genz (1992, 1993) (for algorithm GenzBretz). 
## The default algorithm GenzBretz is randomized
BA_prob_paradox = function(year_range_index, 
                           year_target_index, 
                           setup, 
                           ncores, 
                           B = 100){
  
  do.call(rbind, mclapply(year_range_index, mc.cores = ncores, 
                          FUN = function(xx){ 
                            
                            num_test = replicate(B, prob_paradox(
                              n1 = setup$n[year_target_index],
                              n2 = setup$n[xx],
                              N1 = setup$pops[year_target_index],
                              N2 = setup$pops[xx],
                              mu1 = setup$mean_BA[year_target_index],
                              mu2 = setup$mean_BA[xx],
                              sigma1 = setup$sd_BA[year_target_index],
                              sigma2 = setup$mean_BA[xx]))
                            num_test_mean = mean(num_test)
                            num_test_se = sd(num_test) / 10

                            data.frame(year = setup$year[xx], 
                                       norm_prob = ifelse(num_test_mean < 0, 0, 
                                                          round(num_test_mean, 4)),
                                       num_test_se = num_test_se)
                          }))
}

# build setup data to be used in above function
## includes: means and standard deviations of batting averages 
## for full-time players, and n and N for all years
foo = sample_sizes %>% 
  rename(year = yearID) %>% 
  select(-thres)
para_setup = mean_sd %>% left_join(foo)

## number of years in data set
n_years = nrow(para_setup)

## set cores
ncores = detectCores() - 2
ncores

# compute paradox probabilities with respect to 2023 season
system.time({
  paradox_prob_results_2023 = BA_prob_paradox(
    year_range_index = 1:(n_years-1),
    year_target_index = n_years,
    setup = para_setup,
    ncores = ncores
  )
})
write_csv(paradox_prob_results_2023, file = "paradox_prob_results_2023.csv")

# compute paradox probabilities with respect to 2022 season
system.time({
  paradox_prob_results_2022 = BA_prob_paradox(
    year_range_index = 1:(n_years-2),
    year_target_index = n_years - 1,
    setup = para_setup,
    ncores = ncores
  )
})
write_csv(paradox_prob_results_2022, file = "paradox_prob_results_2022.csv")