
library(tidyverse)
library(orderstats)
library(Pareto)
library(parallel)
library(doParallel)
library(xtable)
library(retrosheet)
library(mvtnorm)
source('source.R')
ncores = detectCores() - 1


Npop = 1e5
Nindex = 2e2
nsize = Npop/Nindex

# means = 60:80
mus = seq(from = 60, to = 80, length = Nindex)
499/20

60 - 1/24.95

60 - 1/24.95 + 251/24.95


# sd = 1,2,5
sds = c(1, 2, 5)

X = rPareto(n = Npop, t = 1, alpha = 1.16)

round_gen = function(X, Npop, Nindex, nsize, sd){
  dat = as.data.frame(cbind(X, ind = rep(Nindex:1, each = nsize) )) %>% 
    arrange(X) 
  
  dat$u = unlist(lapply(1:Npop, function(j){
    pbeta(pPareto(dat$X[j], t = 1, alpha = 1.16), j, Npop-j+1)
  }))
  
  foo = split(dat, f = as.factor(dat$ind))
  
  bar = lapply(foo, function(xx){
    cbind(X = xx$X, u = xx$u, 
          Y = qnorm(qbeta(xx$u, shape1 = 1:nsize, shape2 = nsize:1), 
                    mean = mus[xx$ind], sd = sd))
  })
  
  baz = as.data.frame(do.call(rbind, bar)) %>% arrange(desc(Y))
  baz$ind = rep(Nindex:1, each = nsize)
  baz
}

#round_gen(X = X, Npop = Npop, Nindex = Nindex, nsize = nsize, sd = 1)

round_iter = function(dat, Npop, Nindex, nsize, sd){
  foo = split(dat, f = as.factor(dat$ind))
  bar = lapply(foo, function(xx) xx %>% arrange(X))
  baz = lapply(bar, function(xx){
    cbind(X = xx$X, u = xx$u, 
          Y = qnorm(qbeta(xx$u, shape1 = 1:nsize, shape2 = nsize:1), 
                    mean = mus[xx$ind], sd = sd))
  })
  qux = as.data.frame(do.call(rbind, baz)) %>% arrange(desc(Y))
  qux$ind = rep(Nindex:1, each = nsize)
  qux
}

#dat1 = round_gen(X = X, Npop = Npop, Nindex = Nindex, nsize = nsize, sd = 1)
#dat2 = round_iter(dat = dat1, Npop = Npop, Nindex = Nindex, nsize = nsize, sd = 1)

simulator = function(Npop, Nindex, nsize, sd){
  
  X = rPareto(n = Npop, t = 1, alpha = 1.16)
  dat = round_gen(X = X, Npop = Npop, Nindex = Nindex, nsize = nsize, sd = sd)
  index = dat$ind[which.max(dat$X)]
  round = index
  #k = 0
  
  if(index < Nindex){
    while(index < Nindex){
      #k = k + 1
      dat = round_iter(dat = dat, Npop = Npop, Nindex = Nindex, nsize = nsize, sd = sd)
      index = dat$ind[which.max(dat$X)]
      round = c(round, index)
      #print(index)
      if(length(round) == 25) index = Nindex
    }
  }
  
  round  
}

set.seed(13, "L'Ecuyer-CMRG")

B = 500
nCores = detectCores() - 2

## sd = 1
rounds = function(){
  simulator(Npop = Npop, Nindex = Nindex, nsize = nsize, sd = 1)  
}
out_sd1 = mclapply(1:B, FUN = function(j) rounds(), mc.cores = nCores)
out_sd1

nRounds_sd1 = unlist(lapply(out_sd1, length))
mean(nRounds_sd1)
sd(nRounds_sd1) / sqrt(B)
hist(nRounds_sd1)


## sd = 2
rounds = function(){
  simulator(Npop = Npop, Nindex = Nindex, nsize = nsize, sd = 2)  
}
out_sd2 = mclapply(1:B, FUN = function(j) rounds(), mc.cores = nCores)
out_sd2

nRounds_sd2 = unlist(lapply(out_sd2, length))
mean(nRounds_sd2)
sd(nRounds_sd2) / sqrt(B)
hist(nRounds_sd2)


## sd = 0.5
rounds = function(){
  simulator(Npop = Npop, Nindex = Nindex, nsize = nsize, sd = 0.5)  
}
out_sd0.5 = mclapply(1:B, FUN = function(j) rounds(), mc.cores = nCores)
out_sd0.5

nRounds_sd0.5 = unlist(lapply(out_sd0.5, length))
mean(nRounds_sd0.5)
sd(nRounds_sd0.5) / sqrt(B)
hist(nRounds_sd0.5)


tab = tibble(rounds = c(nRounds_sd2, nRounds_sd1, nRounds_sd0.5), 
             sd = rep(c(2,1,0.5), each = B))
tab

ggplot(tab) + 
  aes(x = rounds) +
  geom_histogram(binwidth = 1) + 
  facet_wrap(~sd, labeller = labeller(sd = function(labels) {
    paste0("\u03C3 = ", labels)
  })) +
  theme_minimal() + 
  labs(title = expression("Number of rounds until " ~ X[(N)] ~ " is assigned label " ~ m))

#save(out_sd0.5, out_sd1, out_sd2, tab, file = "rounds.RData")

#load(file = "rounds.RData")
