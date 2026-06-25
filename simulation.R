library(tidyverse)
condition = function(n1, n2, N1, N2, k, prec = 1e4){
  sapply(seq(from = 0.00001, to = 0.99999, length = prec), 
         FUN = function(x) {pbeta(qbeta(x, n2-k, k+1), n1-k, k+1) - 
             pbeta(qbeta(x, N2-k, k+1), N1-k, k+1)
         })
}

foo = function(k) {
  all(condition(n1 = n1, n2 = n2, k = k, N1 = N1, N2 = N2, prec = 1e3) + 1e-11 >= 0)
}


res <- c()
for (n1 in seq(10, 100, length.out = 10)) {
  for (n2 in seq(10, 100, length.out = 10)) {
    for (N1 in seq(100, 1000, length.out = 10)) {
      for (N2 in seq(100, 1000, length.out = 10)) {
          if (n1/n2 < N1/N2 & N2-N1 > n2-n1) {
            res <- rbind(res, data.frame(n1 = n1, 
                                         n2 = n2, 
                                         N1 = N1, 
                                         N2 = N2, 
                                         j= c(1:(min(n1,n2))), 
                                         ans = sapply(1:(min(n1,n2)), foo)))
        }
      }
    }
  }
}
unique(res$ans)


res2 <- c()
for (n1 in seq(10, 100, length.out = 10)) {
  for (n2 in seq(10, 100, length.out = 10)) {
    for (N1 in seq(100, 1000, length.out = 10)) {
      for (N2 in seq(100, 1000, length.out = 10)) {
        if (n1/n2 < N1/N2) {
          res2 <- rbind(res2, data.frame(n1 = n1, 
                                       n2 = n2, 
                                       N1 = N1, 
                                       N2 = N2, 
                                       j= c(1:(min(n1,n2))), 
                                       ans = sapply(1:(min(n1,n2)), foo)))
        }
      }
    }
  }
}
unique(res2$ans)
res2[(res2$N1 - res2$N2) > (res2$n1 - res2$n2), ]
tail(res2[(res2$N1 - res2$N2) > (res2$n1 - res2$n2), ], 100)


res3 <- c()
for (n1 in seq(10, 100, length.out = 10)) {
  for (n2 in seq(10, 100, length.out = 10)) {
    for (N1 in seq(100, 1000, length.out = 10)) {
      for (N2 in seq(100, 1000, length.out = 10)) {
        if (n1/n2 < N1/N2) {
          res3 <- rbind(res3, data.frame(n1 = n1, 
                                         n2 = n2, 
                                         N1 = N1, 
                                         N2 = N2, 
                                         j= c(0:(min(n1,n2)-1)), 
                                         ans = sapply(0:(min(n1,n2)-1), foo)))
        }
      }
    }
  }
}
unique(res3$ans)
res3[res3$ans == FALSE, ]

res3[(res3$N1 - res3$N2) > (res3$n1 - res3$n2), ]
tail(res3[(res3$N1 - res3$N2) > (res3$n1 - res3$n2), ], 100)

res3[((res3$N1 - res3$N2) > (res3$n1 - res3$n2)) & !res3$ans, ]

library(mvtnorm)
library(tidyverse)
library(doParallel)
library(Pareto)

set.seed(13)
n1 = 100
n2 = 200
N1 = 2000
N2 = 3000
sigma1 = 1
sigma2 = 2
mu1 = -1
mu2 = 1
simu <- function(n1, n2, N1, N2, n1j, n2k, size){
  simu_paradox <- function(n1, n2, N1, N2, n1j, n2k){
    
    Y1 <- sort(rnorm(n1, mean = mu1, sd = sigma1))
    Y2 <- sort(rnorm(n2, mean = mu2, sd = sigma2))
    
    order_pnorm <- function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4){
      p <- pnorm(q, mean = mean, sd = sd)
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    order_pnorm_vec <- function(q, mean = 0, sd = 1){
      q <- sort(q) 
      n <- length(q)
      unlist(lapply(1:n, function(j){
        order_pnorm(q[j], k = j, n = n, mean = mean, sd = sd)
      }))
    }
    
    u1 <- order_pnorm_vec(Y1, mean = mu1, sd = sigma1)
    u2 <- order_pnorm_vec(Y2, mean = mu2, sd = sigma2)
    
    order_Pareto_vec <- function(u, t = 1, alpha = 1.16, npop = 1e4){
      n <- length(u)
      if(length(npop) == 1) npop <- rep(npop, n)
      unlist(lapply(1:n, function(j){
        qPareto(qbeta(u[j], j + npop[j]-n, n + 1 - j), t = t, alpha = alpha)
      }))
    }
    
    # try it out on the percentiles obtained from order_norm_vec
    X1 <- order_Pareto_vec(u1, npop = N1)
    X2 <- order_Pareto_vec(u2, npop = N2)
    return(c(Y1[n1-n1j], Y2[n2-n2k], X1[n1-n1j], X2[n2-n2k]))
  }
  res <- do.call(rbind, mclapply(1:size, function(xx){
    simu_paradox(n1 = n1, n2 = n2, N1 = N1, N2 = N2,n1j = n1j, n2k = n2k)
  }, mc.cores = 7)) 
  res <- as.data.frame(res)
  colnames(res) <- c('Y1', 'Y2', 'X1', 'X2')
  return(nrow(res[res$Y1 < res$Y2 & res$X1 > res$X2,]))
}
result <- simu(n1 = n1, n2 = n1, N1 = N1, N2 = N2, n1j = 0, n2k = 0, size = 100000)
mean(result)/100000

# --- Sigma_11 ---
Sigma11 <- matrix(sigma2^2 / sigma1^2, n1, n1)
diag(Sigma11) <- (sigma2^2 / sigma1^2) + 1

# --- Sigma_22 ---
Sigma22 <- matrix(1, n2 - 1, n2 - 1)
diag(Sigma22) <- 2

# --- Sigma_21 and Sigma_12 ---
Sigma21 <- matrix(sigma2 / sigma1, n2 - 1, n1)
Sigma12 <- t(Sigma21)

# --- Combine blocks into Sigma ---
Sigma <- rbind(
  cbind(Sigma11, Sigma12),
  cbind(Sigma21, Sigma22)
)

upper <- c(rep((mu2-mu1)/sigma1, n1), rep(0, n2-1))
true_r <- replicate(100000, 
                    n2*pmvnorm(lower = -Inf, upper = upper, mean = rep(0, n1+n2-1), sigma = Sigma) - N2/(N1+N2))
mean(true_r)

####### simulation ########

simu_paradox <- function(n1, n2, N1, N2, 
                         mu1, mu2, sigma1, sigma2, 
                         n1j, n2k){
  
  Y1 <- sort(rnorm(n1, mean = mu1, sd = sigma1))
  Y2 <- sort(rnorm(n2, mean = mu2, sd = sigma2))
  
  order_pnorm <- function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4){
    p <- pnorm(q, mean = mean, sd = sd)
    pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
  }
  
  order_pnorm_vec <- function(q, mean = 0, sd = 1){
    q <- sort(q) 
    n <- length(q)
    unlist(lapply(1:n, function(j){
      order_pnorm(q[j], k = j, n = n, mean = mean, sd = sd)
    }))
  }
  
  u1 <- order_pnorm_vec(Y1, mean = mu1, sd = sigma1)
  u2 <- order_pnorm_vec(Y2, mean = mu2, sd = sigma2)
  
  order_Pareto_vec <- function(u, t = 1, alpha = 1.16, npop = 1e4){
    n <- length(u)
    if(length(npop) == 1) npop <- rep(npop, n)
    unlist(lapply(1:n, function(j){
      qPareto(qbeta(u[j], j + npop[j]-n, n + 1 - j), t = t, alpha = alpha)
    }))
  }
  
  # try it out on the percentiles obtained from order_norm_vec
  X1 <- order_Pareto_vec(u1, npop = N1)
  X2 <- order_Pareto_vec(u2, npop = N2)
  return(c(Y1[n1-n1j], Y2[n2-n2k], X1[n1-n1j], X2[n2-n2k]))
}

simu_config <- function(n1, n2, N1, N2,
                        mu1, mu2, sigma1, sigma2,
                        n1j = 0, n2k = 0,
                        size = 10000, cores = 4) {
  
  res <- do.call(rbind, mclapply(
    1:size,
    function(xx) {
      simu_paradox(n1 = n1, n2 = n2, N1 = N1, N2 = N2,
                        mu1 = mu1, mu2 = mu2,
                        sigma1 = sigma1, sigma2 = sigma2,
                        n1j = n1j, n2k = n2k)
    },
    mc.cores = cores
  ))
  
  res <- as.data.frame(res)
  colnames(res) <- c('Y1', 'Y2', 'X1', 'X2')
  return(nrow(res[res$Y1 < res$Y2 & res$X1 > res$X2,]))
}

true_r_config <- function(n1, n2, N1, N2,
                          mu1, mu2, sigma1, sigma2) {
  
  # --- Sigma_11 ---
  Sigma11 <- matrix(sigma2^2 / sigma1^2, n1, n1)
  diag(Sigma11) <- (sigma2^2 / sigma1^2) + 1
  
  # --- Sigma_22 ---
  Sigma22 <- matrix(1, n2 - 1, n2 - 1)
  diag(Sigma22) <- 2
  
  # --- Sigma_21 and Sigma_12 ---
  Sigma21 <- matrix(sigma2 / sigma1, n2 - 1, n1)
  Sigma12 <- t(Sigma21)
  
  # --- Combine blocks ---
  Sigma <- rbind(
    cbind(Sigma11, Sigma12),
    cbind(Sigma21, Sigma22)
  )
  
  upper <- c(rep((mu2 - mu1) / sigma1, n1),
             rep(0, n2 - 1))
  
  val <- n2 * pmvnorm(
    lower = -Inf,
    upper = upper,
    mean  = rep(0, n1 + n2 - 1),
    sigma = Sigma
  ) - N2 / (N1 + N2)
  
  as.numeric(val)
}

# candidate values
n1_vec     <- c(10, 20)
n2_vec     <- c(20, 40)
N1_vec     <- c(200, 400)
N2_vec     <- c(300, 600)
mu1_vec    <- c(-1, 0)
mu2_vec    <- c(0, 1)
sigma1_vec <- c(1, 1.5)
sigma2_vec <- c(2, 3)

param_grid <- expand.grid(
  n1 = n1_vec,
  n2 = n2_vec,
  N1 = N1_vec,
  N2 = N2_vec,
  mu1 = mu1_vec,
  mu2 = mu2_vec,
  sigma1 = sigma1_vec,
  sigma2 = sigma2_vec
)

# enforce constraints: n1/N1 <= n2/N2, sigma1 <= sigma2, mu1 <= mu2
valid_idx <- with(param_grid,
                  (n1 / N1 <= n2 / N2) &
                    (sigma1 <= sigma2) &
                    (mu1 <= mu2))

param_valid <- param_grid[valid_idx, ]
nrow(param_valid)

## number of Monte Carlo reps for each config
MC_size <- 20000  # you can increase later if needed

results_list <- vector("list", nrow(param_valid))

for (i in seq_len(nrow(param_valid))) {
  cat("Running config", i, "of", nrow(param_valid), "...\n")
  p <- param_valid[i, ]
  
  emp_prob <- simu_config(
    n1 = p$n1, n2 = p$n2,
    N1 = p$N1, N2 = p$N2,
    mu1 = p$mu1, mu2 = p$mu2,
    sigma1 = p$sigma1, sigma2 = p$sigma2,
    n1j = 0, n2k = 0,
    size = MC_size, cores = 4
  )
  
  true_val <- replicate(MC_size, true_r_config(
    n1 = p$n1, n2 = p$n2,
    N1 = p$N1, N2 = p$N2,
    mu1 = p$mu1, mu2 = p$mu2,
    sigma1 = p$sigma1, sigma2 = p$sigma2
  ))
  
  results_list[[i]] <- c(p, paradox_prob = emp_prob/MC_size,
                         true_r = mean(true_val))
}

sim_results <- do.call(rbind, results_list)
sim_results <- as.data.frame(sim_results)

print(sim_results)

set.seed(13)
library(tidyverse)
library(doParallel)
library(Pareto)
simu <- function(n1, n2, N1, N2, n1j, n2k, size){
  simu_paradox <- function(n1, n2, N1, N2, n1j, n2k){
    
    Y1 <- sort(rnorm(n1))
    Y2 <- sort(rnorm(n2))
    
    order_pnorm <- function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4){
      p <- pnorm(q, mean = mean, sd = sd)
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    order_pnorm_vec <- function(q, mean = 0, sd = 1){
      q <- sort(q) 
      n <- length(q)
      unlist(lapply(1:n, function(j){
        order_pnorm(q[j], k = j, n = n, mean = mean, sd = sd)
      }))
    }
    
    u1 <- order_pnorm_vec(Y1)
    u2 <- order_pnorm_vec(Y2)
    
    order_Pareto_vec <- function(u, t = 1, alpha = 1.16, npop = 1e4){
      n <- length(u)
      if(length(npop) == 1) npop <- rep(npop, n)
      unlist(lapply(1:n, function(j){
        qPareto(qbeta(u[j], j + npop[j]-n, n + 1 - j), t = t, alpha = alpha)
      }))
    }
    
    # try it out on the percentiles obtained from order_norm_vec
    X1 <- order_Pareto_vec(u1, npop = N1)
    X2 <- order_Pareto_vec(u2, npop = N2)
    return(c(Y1[n1-n1j], Y2[n2-n2k], X1[n1-n1j], X2[n2-n2k]))
  }
  res <- do.call(rbind, mclapply(1:size, function(xx){
    simu_paradox(n1 = n1, n2 = n2, N1 = N1, N2 = N2,n1j = n1j, n2k = n2k)
  }, mc.cores = 7)) 
  res <- as.data.frame(res)
  colnames(res) <- c('Y1', 'Y2', 'X1', 'X2')
  return(nrow(res[res$Y1 < res$Y2 & res$X1 > res$X2,]))
}


res <- c()
c = 0.5
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_0.5 <- res

res <- c()
c = 1
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_1 <- res

res <- c()
c = 0.2
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_0.2 <- res

res <- c()
c = 0.1
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_0.1 <- res

res <- c()
c = 0.3
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_0.3 <- res

res <- c()
c = 0.4
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_0.4 <- res

res <- c()
c = 0.6
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_0.6 <- res

res <- c()
c = 0.7
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_0.7 <- res

res <- c()
c = 0.8
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_0.8 <- res

res <- c()
c = 0.9
for (N1 in seq(100,2000,100)) {
  result <- replicate(10, simu(n1 = 60, n2 = 40, N1 = N1, N2 = N1*c, n1j = 1, n2k = 1, size = 1000))
  res <- c(res, mean(result)/1000)
}
res_0.9 <- res

library(ggplot2)
library(dplyr)

res <- rbind(res_0.1, res_0.2, res_0.3, res_0.4, res_0.5, res_0.6, res_0.7, res_0.8, res_0.9, res_1)
res <- as.data.frame(res)
res$c <- paste0("c = ", seq(0.1, 1.0, by = 0.1))

df_long <- res %>%
  pivot_longer(
    cols = -c,                   # All columns except 'c' are values
    names_to = "x",              # Create 'x' column from original column names
    values_to = "value"          # Store actual values in 'value'
  ) %>%
  mutate(x = as.numeric(gsub("V", "", x)))
df_long$x <- as.numeric(df_long$x)*100

ggplot(df_long, aes(x = x, y = value, color = c, fill = c)) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +  # smoothed line and band
  labs(title = "Smoothed Vectors with Confidence Bands",
       x = "N1", y = "Prob", color = "c", fill = "c") +
  theme_minimal()

set.seed(13)

library(tidyverse)
library(parallel)
library(Pareto)
library(ggplot2)
library(dplyr)

# ============================================================
# 1. Simulation function
# ============================================================

simu <- function(n1, n2, N1, N2, n1j, n2k, size, mc.cores = 7) {
  
  simu_paradox <- function(n1, n2, N1, N2, n1j, n2k) {
    
    Y1 <- sort(rnorm(n1))
    Y2 <- sort(rnorm(n2))
    
    order_pnorm <- function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4) {
      p <- pnorm(q, mean = mean, sd = sd)
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    order_pnorm_vec <- function(q, mean = 0, sd = 1) {
      q <- sort(q)
      n <- length(q)
      
      unlist(lapply(1:n, function(j) {
        order_pnorm(
          q[j],
          k = j,
          n = n,
          mean = mean,
          sd = sd
        )
      }))
    }
    
    u1 <- order_pnorm_vec(Y1)
    u2 <- order_pnorm_vec(Y2)
    
    order_Pareto_vec <- function(u, t = 1, alpha = 1.16, npop = 1e4) {
      n <- length(u)
      if (length(npop) == 1) npop <- rep(npop, n)
      
      unlist(lapply(1:n, function(j) {
        qPareto(
          qbeta(u[j], j + npop[j] - n, n + 1 - j),
          t = t,
          alpha = alpha
        )
      }))
    }
    
    X1 <- order_Pareto_vec(u1, npop = N1)
    X2 <- order_Pareto_vec(u2, npop = N2)
    
    return(c(
      Y1[n1 - n1j],
      Y2[n2 - n2k],
      X1[n1 - n1j],
      X2[n2 - n2k]
    ))
  }
  
  res <- do.call(
    rbind,
    mclapply(
      1:size,
      function(xx) {
        simu_paradox(
          n1 = n1,
          n2 = n2,
          N1 = N1,
          N2 = N2,
          n1j = n1j,
          n2k = n2k
        )
      },
      mc.cores = mc.cores
    )
  )
  
  res <- as.data.frame(res)
  colnames(res) <- c("Y1", "Y2", "X1", "X2")
  
  paradox_count <- nrow(res[res$Y1 < res$Y2 & res$X1 > res$X2, ])
  
  return(paradox_count)
}


# ============================================================
# 2. Wrapper for one achiever rank
# ============================================================
# rank = 1 means largest achiever
# rank = 2 means 2nd largest achiever
# rank = 6 means 6th largest achiever

run_one_rank <- function(rank,
                         n1 = 60,
                         n2 = 40,
                         N1_grid = seq(100, 2000, 100),
                         c_grid = seq(0.1, 1.0, by = 0.1),
                         size = 1000,
                         nrep = 10,
                         mc.cores = 7) {
  
  # Convert achiever rank to your original indexing
  # largest: n1j = 0
  # 2nd largest: n1j = 1
  # 3rd largest: n1j = 2
  offset <- rank - 1
  
  output <- list()
  
  counter <- 1
  
  for (c_val in c_grid) {
    
    for (N1 in N1_grid) {
      
      N2 <- N1 * c_val
      
      result <- replicate(
        nrep,
        simu(
          n1 = n1,
          n2 = n2,
          N1 = N1,
          N2 = N2,
          n1j = offset,
          n2k = offset,
          size = size,
          mc.cores = mc.cores
        )
      )
      
      prob <- mean(result) / size
      
      output[[counter]] <- data.frame(
        rank = rank,
        achiever = paste0(rank, ordinal_suffix(rank), " largest achiever"),
        c = c_val,
        c_label = paste0("c = ", c_val),
        N1 = N1,
        N2 = N2,
        prob = prob
      )
      
      counter <- counter + 1
    }
  }
  
  bind_rows(output)
}


# Helper function for nicer labels
ordinal_suffix <- function(x) {
  if (x == 1) return("st")
  if (x == 2) return("nd")
  if (x == 3) return("rd")
  return("th")
}


# ============================================================
# 3. Run simulations for 1st to 6th largest achievers
# ============================================================

N1_grid <- seq(100, 2000, 100)
c_grid <- seq(0.1, 1.0, by = 0.1)

all_results <- bind_rows(
  lapply(1:6, function(r) {
    cat("Running simulation for rank =", r, "\n")
    
    run_one_rank(
      rank = r,
      n1 = 60,
      n2 = 40,
      N1_grid = N1_grid,
      c_grid = c_grid,
      size = 1000,
      nrep = 10,
      mc.cores = 7
    )
  })
)

# ============================================================
# 4. Make 2 x 2 plot for 1st to 4th largest achievers only
# ============================================================

plot_results <- all_results %>%
  filter(rank %in% 1:4)

p <- ggplot(plot_results, aes(x = N1, y = prob, color = c_label, fill = c_label)) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    alpha = 0.2,
    linewidth = 0.9
  ) +
  facet_wrap(~ achiever, nrow = 2, ncol = 2) +
  labs(
    title = "Paradox Probability for 1st to 4th Largest Achievers",
    x = expression(N[1]),
    y = "Paradox Probability",
    color = expression(c == N[2] / N[1]),
    fill = expression(c == N[2] / N[1])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
  )

print(p)

ggsave(
  filename = "paradox_probability_1st_to_4th_achievers.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)



set.seed(13)
library(tidyverse)
library(doParallel)
library(Pareto)
simu <- function(n1, n2, N1, N2, n1j, n2k, size){
  simu_paradox <- function(n1, n2, N1, N2, n1j, n2k){
    
    Y1 <- sort(rnorm(n1))
    Y2 <- sort(rnorm(n2))
    
    order_pnorm <- function(q = 0, mean = 0, sd = 1, k = 1, n = 1e4){
      p <- pnorm(q, mean = mean, sd = sd)
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    order_pnorm_vec <- function(q, mean = 0, sd = 1){
      q <- sort(q) 
      n <- length(q)
      unlist(lapply(1:n, function(j){
        order_pnorm(q[j], k = j, n = n, mean = mean, sd = sd)
      }))
    }
    
    u1 <- order_pnorm_vec(Y1)
    u2 <- order_pnorm_vec(Y2)
    
    order_Pareto_vec <- function(u, t = 1, alpha = 1.16, npop = 1e4){
      n <- length(u)
      if(length(npop) == 1) npop <- rep(npop, n)
      unlist(lapply(1:n, function(j){
        qPareto(qbeta(u[j], j + npop[j]-n, n + 1 - j), t = t, alpha = alpha)
      }))
    }
    
    # try it out on the percentiles obtained from order_norm_vec
    X1 <- order_Pareto_vec(u1, npop = N1)
    X2 <- order_Pareto_vec(u2, npop = N2)
    return(c(Y1[n1-n1j], Y2[n2-n2k], X1[n1-n1j], X2[n2-n2k]))
  }
  res <- do.call(rbind, mclapply(1:size, function(xx){
    simu_paradox(n1 = n1, n2 = n2, N1 = N1, N2 = N2,n1j = n1j, n2k = n2k)
  }, mc.cores = 7)) 
  res <- as.data.frame(res)
  colnames(res) <- c('Y1', 'Y2', 'X1', 'X2')
  return(nrow(res[res$Y1 < res$Y2 & res$X1 > res$X2,]))
}

## n1 and n2 fixed
res <- matrix(0, nrow = 10, ncol = 10)
for (N1 in seq(200,2000,200)) {
  for (N2 in seq(200,2000,200)) {
    result <- replicate(10, simu(n1 = 40, n2 = 60, N1 = N1, N2 = N2, n1j = 1, n2k = 1, size = 1000))
    res[N1/200, N2/200] <- mean(result)/1000
  }
}

library(ggplot2)
library(reshape2)  # for melting the matrix
library(viridis)
library(akima)       
# Convert matrix to long format
df <- melt(res)
colnames(df) <- c('N1', 'N2', 'prob')
df$N1 <- df$N1*200
df$N2 <- df$N2*200

interp_out <- with(df, akima::interp(
  x = N1, y = N2, z = prob,
  nx = 300, ny = 300    
))

smooth_df_1 <- data.frame(
  N1   = rep(interp_out$x, times = length(interp_out$y)),
  N2   = rep(interp_out$y, each  = length(interp_out$x)),
  prob = as.vector(interp_out$z)
)

## N1 and N2 fixed
res <- matrix(0, nrow = 10, ncol = 10)
for (n1 in seq(20,200,20)) {
  for (n2 in seq(20,200,20)) {
    result <- replicate(10, simu(n1 = n1, n2 = n2, N1 = 1e3, N2 = 1e3, n1j = 1, n2k = 1, size = 1000))
    res[n1/20, n2/20] <- mean(result)/1000
  }
}

# Convert matrix to long format
df <- melt(res)
colnames(df) <- c('n1', 'n2', 'prob')
df$n1 <- df$n1*20
df$n2 <- df$n2*20
interp_out <- with(df, akima::interp(
  x = n1, y = n2, z = prob,
  nx = 300, ny = 300    
))

smooth_df_2 <- data.frame(
  n1   = rep(interp_out$x, times = length(interp_out$y)),
  n2   = rep(interp_out$y, each  = length(interp_out$x)),
  prob = as.vector(interp_out$z)
)

## n1 and N2 fixed
res <- matrix(0, nrow = 10, ncol = 10)
for (n2 in seq(20,200,20)) {
  for (N1 in seq(200,2000,200)) {
    result <- replicate(10, simu(n1 = 40, n2 = n2, N1 = N1, N2 = 1e3, n1j = 1, n2k = 1, size = 1000))
    res[n2/20, N1/200] <- mean(result)/1000
  }
}

# Convert matrix to long format
df <- melt(res)
colnames(df) <- c('n2', 'N1', 'prob')
df$n2 <- df$n2*20
df$N1 <- df$N1*200
interp_out <- with(df, akima::interp(
  x = n2, y = N1, z = prob,
  nx = 300, ny = 300    
))

smooth_df_3 <- data.frame(
  n2   = rep(interp_out$x, times = length(interp_out$y)),
  N1   = rep(interp_out$y, each  = length(interp_out$x)),
  prob = as.vector(interp_out$z)
)

## n2 and N1 fixed
res <- matrix(0, nrow = 10, ncol = 10)
for (n1 in seq(20,200,20)) {
  for (N2 in seq(200,2000,200)) {
    result <- replicate(10, simu(n1 = n1, n2 = 60, N1 = 1e3, N2 = N2, n1j = 1, n2k = 1, size = 1000))
    res[n1/20, N2/200] <- mean(result)/1000
  }
}

# Convert matrix to long format
df <- melt(res)
colnames(df) <- c('n1', 'N2', 'prob')
df$n1 <- df$n1*20
df$N2 <- df$N2*200
interp_out <- with(df, akima::interp(
  x = n1, y = N2, z = prob,
  nx = 300, ny = 300    
))

smooth_df_4 <- data.frame(
  n1   = rep(interp_out$x, times = length(interp_out$y)),
  N2   = rep(interp_out$y, each  = length(interp_out$x)),
  prob = as.vector(interp_out$z)
)
library(patchwork)

all_probs <- c(smooth_df_1$prob, smooth_df_2$prob, smooth_df_3$prob, smooth_df_4$prob)
mid_val <- median(all_probs, na.rm = TRUE)
lim_range <- range(all_probs, na.rm = TRUE)

p1 <- ggplot(smooth_df_1, aes(N1, N2, fill = prob)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(
    low = "white", mid = "lightblue", high = "blue",
    midpoint = mid_val,
    limits = lim_range,
    name = "Prob"
  ) +
  theme_minimal() +
  labs(x = "N1", y = "N2", title = '(A)')

p2 <- ggplot(smooth_df_2, aes(n1, n2, fill = prob)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(
    low = "white", mid = "lightblue", high = "blue",
    midpoint = mid_val,
    limits = lim_range,
    name = "Prob"
  ) +
  theme_minimal() +
  labs(x = "n1", y = "n2", title = '(B)')

p3 <- ggplot(smooth_df_3, aes(n2, N1, fill = prob)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(
    low = "white", mid = "lightblue", high = "blue",
    midpoint = mid_val,
    limits = lim_range,
    name = "Prob"
  )  +
  theme_minimal() +
  labs(x = "n2", y = "N1", title = '(C)')

p4 <- ggplot(smooth_df_4, aes(n1, N2, fill = prob)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_gradient2(
    low = "white", mid = "lightblue", high = "blue",
    midpoint = mid_val,
    limits = lim_range,
    name = "Prob"
  ) +
  theme_minimal() +
  labs(x = "n1", y = "N2", title = '(D)')

(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

set.seed(13)

library(tidyverse)
library(parallel)
library(Pareto)
library(ggplot2)
library(dplyr)

# ============================================================
# 1. Simulation function with Gamma achievement distribution
# ============================================================

simu <- function(n1, n2, N1, N2, n1j, n2k, size,
                 gamma_shape = 2,
                 gamma_rate = 1,
                 mc.cores = 7) {
  
  simu_paradox <- function(n1, n2, N1, N2, n1j, n2k) {
    
    # --------------------------------------------------------
    # Achievement distribution:
    # Change from Normal to Gamma
    # Gamma(shape = 2, rate = 1) is noticeably right-skewed.
    # For stronger skewness, try shape = 1 or shape = 0.75.
    # --------------------------------------------------------
    Y1 <- sort(rgamma(n1, shape = gamma_shape, rate = gamma_rate))
    Y2 <- sort(rgamma(n2, shape = gamma_shape, rate = gamma_rate))
    
    order_pgamma <- function(q = 0,
                             shape = gamma_shape,
                             rate = gamma_rate,
                             k = 1,
                             n = 1e4) {
      p <- pgamma(q, shape = shape, rate = rate)
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    order_pgamma_vec <- function(q,
                                 shape = gamma_shape,
                                 rate = gamma_rate) {
      q <- sort(q)
      n <- length(q)
      
      unlist(lapply(1:n, function(j) {
        order_pgamma(
          q[j],
          k = j,
          n = n,
          shape = shape,
          rate = rate
        )
      }))
    }
    
    u1 <- order_pgamma_vec(Y1)
    u2 <- order_pgamma_vec(Y2)
    
    order_Pareto_vec <- function(u, t = 1, alpha = 1.16, npop = 1e4) {
      n <- length(u)
      if (length(npop) == 1) npop <- rep(npop, n)
      
      unlist(lapply(1:n, function(j) {
        qPareto(
          qbeta(u[j], j + npop[j] - n, n + 1 - j),
          t = t,
          alpha = alpha
        )
      }))
    }
    
    X1 <- order_Pareto_vec(u1, npop = N1)
    X2 <- order_Pareto_vec(u2, npop = N2)
    
    return(c(
      Y1[n1 - n1j],
      Y2[n2 - n2k],
      X1[n1 - n1j],
      X2[n2 - n2k]
    ))
  }
  
  res <- do.call(
    rbind,
    mclapply(
      1:size,
      function(xx) {
        simu_paradox(
          n1 = n1,
          n2 = n2,
          N1 = N1,
          N2 = N2,
          n1j = n1j,
          n2k = n2k
        )
      },
      mc.cores = mc.cores
    )
  )
  
  res <- as.data.frame(res)
  colnames(res) <- c("Y1", "Y2", "X1", "X2")
  
  paradox_count <- nrow(res[res$Y1 < res$Y2 & res$X1 > res$X2, ])
  
  return(paradox_count)
}


# ============================================================
# 2. Wrapper for one achiever rank
# ============================================================
# rank = 1 means largest achiever
# rank = 2 means 2nd largest achiever
# rank = 3 means 3rd largest achiever
# rank = 4 means 4th largest achiever

run_one_rank <- function(rank,
                         n1 = 60,
                         n2 = 40,
                         N1_grid = seq(100, 2000, 100),
                         c_grid = seq(0.1, 1.0, by = 0.1),
                         size = 1000,
                         nrep = 10,
                         gamma_shape = 2,
                         gamma_rate = 1,
                         mc.cores = 7) {
  
  # Convert achiever rank to your original indexing
  # largest: n1j = 0
  # 2nd largest: n1j = 1
  # 3rd largest: n1j = 2
  # 4th largest: n1j = 3
  offset <- rank - 1
  
  output <- list()
  counter <- 1
  
  for (c_val in c_grid) {
    
    for (N1 in N1_grid) {
      
      N2 <- N1 * c_val
      
      result <- replicate(
        nrep,
        simu(
          n1 = n1,
          n2 = n2,
          N1 = N1,
          N2 = N2,
          n1j = offset,
          n2k = offset,
          size = size,
          gamma_shape = gamma_shape,
          gamma_rate = gamma_rate,
          mc.cores = mc.cores
        )
      )
      
      prob <- mean(result) / size
      
      output[[counter]] <- data.frame(
        rank = rank,
        achiever = paste0(rank, ordinal_suffix(rank), " largest achiever"),
        c = c_val,
        c_label = paste0("c = ", c_val),
        N1 = N1,
        N2 = N2,
        prob = prob
      )
      
      counter <- counter + 1
    }
  }
  
  bind_rows(output)
}


# ============================================================
# 3. Helper function for nicer labels
# ============================================================

ordinal_suffix <- function(x) {
  if (x == 1) return("st")
  if (x == 2) return("nd")
  if (x == 3) return("rd")
  return("th")
}


# ============================================================
# 4. Run simulations for 1st to 6th largest achievers
# ============================================================

N1_grid <- seq(100, 2000, 100)
c_grid <- seq(0.1, 1.0, by = 0.1)

# Gamma distribution parameters
# Gamma(shape = 2, rate = 1) is right-skewed.
# For stronger skewness, use gamma_shape = 1 or 0.75.
gamma_shape <- 2
gamma_rate <- 1

all_results_gamma <- bind_rows(
  lapply(1:6, function(r) {
    cat("Running Gamma simulation for rank =", r, "\n")
    
    run_one_rank(
      rank = r,
      n1 = 60,
      n2 = 40,
      N1_grid = N1_grid,
      c_grid = c_grid,
      size = 1000,
      nrep = 10,
      gamma_shape = gamma_shape,
      gamma_rate = gamma_rate,
      mc.cores = 7
    )
  })
)


# ============================================================
# 5. Make 2 x 2 plot for 1st to 4th largest achievers only
# ============================================================

plot_results_gamma <- all_results_gamma %>%
  filter(rank %in% 1:4)

p_gamma <- ggplot(
  plot_results_gamma,
  aes(x = N1, y = prob, color = c_label, fill = c_label)
) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    alpha = 0.2,
    linewidth = 0.9
  ) +
  facet_wrap(~ achiever, nrow = 2, ncol = 2) +
  labs(
    title = "Paradox Probability for 1st to 4th Largest Achievers: Gamma Achievement Distribution",
    subtitle = paste0(
      "Y ~ Gamma(shape = ",
      gamma_shape,
      ", rate = ",
      gamma_rate,
      ")"
    ),
    x = expression(N[1]),
    y = "Paradox Probability",
    color = expression(c == N[2] / N[1]),
    fill = expression(c == N[2] / N[1])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

print(p_gamma)

ggsave(
  filename = "paradox_probability_1st_to_4th_achievers_gamma.png",
  plot = p_gamma,
  width = 10,
  height = 8,
  dpi = 300
)

# ============================================================
# Sensitivity analysis for 1st--4th largest achievers
# Normal and Gamma achievement distributions
# ============================================================

set.seed(13)

library(tidyverse)
library(parallel)
library(Pareto)
library(ggplot2)
library(reshape2)
library(akima)
library(patchwork)

# ============================================================
# 0. Helper: ordinal labels
# ============================================================

ordinal_suffix <- function(x) {
  if (x == 1) return("st")
  if (x == 2) return("nd")
  if (x == 3) return("rd")
  return("th")
}


# ============================================================
# 1. Main simulation function
# ============================================================
# rank = 1 means largest achiever
# rank = 2 means second-largest achiever
# rank = 3 means third-largest achiever
# rank = 4 means fourth-largest achiever
#
# Your original code used:
#   Y1[n1 - n1j]
#   Y2[n2 - n2k]
#
# Therefore:
#   rank = 1  -> n1j = n2k = 0
#   rank = 2  -> n1j = n2k = 1
#   rank = 3  -> n1j = n2k = 2
#   rank = 4  -> n1j = n2k = 3
# ============================================================

simu <- function(n1, n2, N1, N2,
                 rank = 2,
                 size = 1000,
                 y_dist = c("normal", "gamma"),
                 gamma_shape = 2,
                 gamma_rate = 1,
                 pareto_t = 1,
                 pareto_alpha = 1.16,
                 mc.cores = 7) {
  
  y_dist <- match.arg(y_dist)
  
  n1j <- rank - 1
  n2k <- rank - 1
  
  simu_paradox <- function(n1, n2, N1, N2, n1j, n2k) {
    
    # --------------------------------------------------------
    # Generate achievement samples Y1 and Y2
    # --------------------------------------------------------
    if (y_dist == "normal") {
      Y1 <- sort(rnorm(n1))
      Y2 <- sort(rnorm(n2))
    } else if (y_dist == "gamma") {
      Y1 <- sort(rgamma(n1, shape = gamma_shape, rate = gamma_rate))
      Y2 <- sort(rgamma(n2, shape = gamma_shape, rate = gamma_rate))
    }
    
    # --------------------------------------------------------
    # CDF of achievement distribution
    # --------------------------------------------------------
    FY <- function(q) {
      if (y_dist == "normal") {
        return(pnorm(q, mean = 0, sd = 1))
      } else if (y_dist == "gamma") {
        return(pgamma(q, shape = gamma_shape, rate = gamma_rate))
      }
    }
    
    # --------------------------------------------------------
    # Order-statistic probability transform
    # Same role as your order_pnorm/order_pnorm_vec,
    # but generalized to Normal or Gamma Y.
    # --------------------------------------------------------
    order_pY <- function(q, k, n) {
      p <- FY(q)
      pbinom(k - 1, prob = p, size = n, lower.tail = FALSE)
    }
    
    order_pY_vec <- function(q) {
      q <- sort(q)
      n <- length(q)
      
      unlist(lapply(1:n, function(j) {
        order_pY(q[j], k = j, n = n)
      }))
    }
    
    u1 <- order_pY_vec(Y1)
    u2 <- order_pY_vec(Y2)
    
    # --------------------------------------------------------
    # Map order-statistic percentiles to Pareto latent aptitude X
    # --------------------------------------------------------
    order_Pareto_vec <- function(u, t = 1, alpha = 1.16, npop = 1e4) {
      n <- length(u)
      
      if (length(npop) == 1) {
        npop <- rep(npop, n)
      }
      
      unlist(lapply(1:n, function(j) {
        qPareto(
          qbeta(u[j], j + npop[j] - n, n + 1 - j),
          t = t,
          alpha = alpha
        )
      }))
    }
    
    X1 <- order_Pareto_vec(
      u1,
      t = pareto_t,
      alpha = pareto_alpha,
      npop = N1
    )
    
    X2 <- order_Pareto_vec(
      u2,
      t = pareto_t,
      alpha = pareto_alpha,
      npop = N2
    )
    
    return(c(
      Y1[n1 - n1j],
      Y2[n2 - n2k],
      X1[n1 - n1j],
      X2[n2 - n2k]
    ))
  }
  
  res <- do.call(
    rbind,
    mclapply(
      1:size,
      function(xx) {
        simu_paradox(
          n1 = n1,
          n2 = n2,
          N1 = N1,
          N2 = N2,
          n1j = n1j,
          n2k = n2k
        )
      },
      mc.cores = mc.cores
    )
  )
  
  res <- as.data.frame(res)
  colnames(res) <- c("Y1", "Y2", "X1", "X2")
  
  paradox_count <- nrow(
    res[res$Y1 < res$Y2 & res$X1 > res$X2, ]
  )
  
  return(paradox_count)
}


# ============================================================
# 2. Interpolation helper
# ============================================================

make_smooth_df <- function(df, xvar, yvar, zvar = "prob",
                           nx = 300, ny = 300) {
  
  interp_out <- akima::interp(
    x = df[[xvar]],
    y = df[[yvar]],
    z = df[[zvar]],
    nx = nx,
    ny = ny
  )
  
  smooth_df <- data.frame(
    x = rep(interp_out$x, times = length(interp_out$y)),
    y = rep(interp_out$y, each = length(interp_out$x)),
    prob = as.vector(interp_out$z)
  )
  
  colnames(smooth_df)[1:2] <- c(xvar, yvar)
  
  return(smooth_df)
}


# ============================================================
# 3. Run four sensitivity surfaces for one achiever rank
# ============================================================

run_one_rank <- function(rank = 2,
                         y_dist = c("normal", "gamma"),
                         gamma_shape = 2,
                         gamma_rate = 1,
                         size = 1000,
                         nrep = 10,
                         mc.cores = 7) {
  
  y_dist <- match.arg(y_dist)
  
  cat("====================================================\n")
  cat("Running rank =", rank, "\n")
  cat("Y distribution =", y_dist, "\n")
  if (y_dist == "gamma") {
    cat("Gamma shape =", gamma_shape, ", rate =", gamma_rate, "\n")
  }
  cat("====================================================\n")
  
  # ----------------------------------------------------------
  # (A) n1 and n2 fixed; vary N1 and N2
  # ----------------------------------------------------------
  
  cat("Running panel (A): fixed n1, n2; varying N1, N2\n")
  
  res <- matrix(0, nrow = 10, ncol = 10)
  
  for (N1 in seq(200, 2000, 200)) {
    for (N2 in seq(200, 2000, 200)) {
      
      result <- replicate(
        nrep,
        simu(
          n1 = 40,
          n2 = 60,
          N1 = N1,
          N2 = N2,
          rank = rank,
          size = size,
          y_dist = y_dist,
          gamma_shape = gamma_shape,
          gamma_rate = gamma_rate,
          mc.cores = mc.cores
        )
      )
      
      res[N1 / 200, N2 / 200] <- mean(result) / size
    }
  }
  
  df <- reshape2::melt(res)
  colnames(df) <- c("N1", "N2", "prob")
  df$N1 <- df$N1 * 200
  df$N2 <- df$N2 * 200
  
  smooth_df_1 <- make_smooth_df(df, xvar = "N1", yvar = "N2")
  
  
  # ----------------------------------------------------------
  # (B) N1 and N2 fixed; vary n1 and n2
  # ----------------------------------------------------------
  
  cat("Running panel (B): fixed N1, N2; varying n1, n2\n")
  
  res <- matrix(0, nrow = 10, ncol = 10)
  
  for (n1 in seq(20, 200, 20)) {
    for (n2 in seq(20, 200, 20)) {
      
      result <- replicate(
        nrep,
        simu(
          n1 = n1,
          n2 = n2,
          N1 = 1e3,
          N2 = 1e3,
          rank = rank,
          size = size,
          y_dist = y_dist,
          gamma_shape = gamma_shape,
          gamma_rate = gamma_rate,
          mc.cores = mc.cores
        )
      )
      
      res[n1 / 20, n2 / 20] <- mean(result) / size
    }
  }
  
  df <- reshape2::melt(res)
  colnames(df) <- c("n1", "n2", "prob")
  df$n1 <- df$n1 * 20
  df$n2 <- df$n2 * 20
  
  smooth_df_2 <- make_smooth_df(df, xvar = "n1", yvar = "n2")
  
  
  # ----------------------------------------------------------
  # (C) n1 and N2 fixed; vary n2 and N1
  # ----------------------------------------------------------
  
  cat("Running panel (C): fixed n1, N2; varying n2, N1\n")
  
  res <- matrix(0, nrow = 10, ncol = 10)
  
  for (n2 in seq(20, 200, 20)) {
    for (N1 in seq(200, 2000, 200)) {
      
      result <- replicate(
        nrep,
        simu(
          n1 = 40,
          n2 = n2,
          N1 = N1,
          N2 = 1e3,
          rank = rank,
          size = size,
          y_dist = y_dist,
          gamma_shape = gamma_shape,
          gamma_rate = gamma_rate,
          mc.cores = mc.cores
        )
      )
      
      res[n2 / 20, N1 / 200] <- mean(result) / size
    }
  }
  
  df <- reshape2::melt(res)
  colnames(df) <- c("n2", "N1", "prob")
  df$n2 <- df$n2 * 20
  df$N1 <- df$N1 * 200
  
  smooth_df_3 <- make_smooth_df(df, xvar = "n2", yvar = "N1")
  
  
  # ----------------------------------------------------------
  # (D) n2 and N1 fixed; vary n1 and N2
  # ----------------------------------------------------------
  
  cat("Running panel (D): fixed n2, N1; varying n1, N2\n")
  
  res <- matrix(0, nrow = 10, ncol = 10)
  
  for (n1 in seq(20, 200, 20)) {
    for (N2 in seq(200, 2000, 200)) {
      
      result <- replicate(
        nrep,
        simu(
          n1 = n1,
          n2 = 60,
          N1 = 1e3,
          N2 = N2,
          rank = rank,
          size = size,
          y_dist = y_dist,
          gamma_shape = gamma_shape,
          gamma_rate = gamma_rate,
          mc.cores = mc.cores
        )
      )
      
      res[n1 / 20, N2 / 200] <- mean(result) / size
    }
  }
  
  df <- reshape2::melt(res)
  colnames(df) <- c("n1", "N2", "prob")
  df$n1 <- df$n1 * 20
  df$N2 <- df$N2 * 200
  
  smooth_df_4 <- make_smooth_df(df, xvar = "n1", yvar = "N2")
  
  return(list(
    smooth_df_1 = smooth_df_1,
    smooth_df_2 = smooth_df_2,
    smooth_df_3 = smooth_df_3,
    smooth_df_4 = smooth_df_4
  ))
}


# ============================================================
# 4. Plot function for one rank
# ============================================================

plot_one_rank <- function(out,
                          rank = 2,
                          y_label = "Normal") {
  
  smooth_df_1 <- out$smooth_df_1
  smooth_df_2 <- out$smooth_df_2
  smooth_df_3 <- out$smooth_df_3
  smooth_df_4 <- out$smooth_df_4
  
  all_probs <- c(
    smooth_df_1$prob,
    smooth_df_2$prob,
    smooth_df_3$prob,
    smooth_df_4$prob
  )
  
  mid_val <- median(all_probs, na.rm = TRUE)
  lim_range <- range(all_probs, na.rm = TRUE)
  
  rank_label <- paste0(rank, ordinal_suffix(rank), " largest achiever")
  
  common_theme <- theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      )
    )
  
  p1 <- ggplot(smooth_df_1, aes(N1, N2, fill = prob)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "white",
      mid = "lightblue",
      high = "blue",
      midpoint = mid_val,
      limits = lim_range,
      name = "Prob"
    ) +
    common_theme +
    labs(
      x = expression(N[1]),
      y = expression(N[2]),
      title = "(A)"
    )
  
  p2 <- ggplot(smooth_df_2, aes(n1, n2, fill = prob)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "white",
      mid = "lightblue",
      high = "blue",
      midpoint = mid_val,
      limits = lim_range,
      name = "Prob"
    ) +
    common_theme +
    labs(
      x = expression(n[1]),
      y = expression(n[2]),
      title = "(B)"
    )
  
  p3 <- ggplot(smooth_df_3, aes(n2, N1, fill = prob)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "white",
      mid = "lightblue",
      high = "blue",
      midpoint = mid_val,
      limits = lim_range,
      name = "Prob"
    ) +
    common_theme +
    labs(
      x = expression(n[2]),
      y = expression(N[1]),
      title = "(C)"
    )
  
  p4 <- ggplot(smooth_df_4, aes(n1, N2, fill = prob)) +
    geom_raster(interpolate = TRUE) +
    scale_fill_gradient2(
      low = "white",
      mid = "lightblue",
      high = "blue",
      midpoint = mid_val,
      limits = lim_range,
      name = "Prob"
    ) +
    common_theme +
    labs(
      x = expression(n[1]),
      y = expression(N[2]),
      title = "(D)"
    )
  
  p <- ((p1 + p2) / (p3 + p4)) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(
        "Sensitivity Analysis: ",
        rank_label,
        ", ",
        y_label,
        " Achievement Distribution"
      ),
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 14
        )
      )
    )
  
  return(p)
}

# ============================================================
# 5. Run Normal case: Y ~ N(0,1)
# ============================================================

normal_results <- list()
normal_plots <- list()

for (r in 1:4) {
  
  normal_results[[r]] <- run_one_rank(
    rank = r,
    y_dist = "normal",
    size = 1000,
    nrep = 10,
    mc.cores = 7
  )
  
  normal_plots[[r]] <- plot_one_rank(
    out = normal_results[[r]],
    rank = r,
    y_label = "Normal"
  )
  
  print(normal_plots[[r]])
  
  ggsave(
    filename = paste0("sensitivity_normal_rank_", r, ".png"),
    plot = normal_plots[[r]],
    width = 10,
    height = 8,
    dpi = 300
  )
}

normal_plots <- list()

plot_one_rank <- function(out,
                          rank = 2,
                          y_label = "Normal") {
  
  smooth_df_1 <- out$smooth_df_1
  smooth_df_2 <- out$smooth_df_2
  smooth_df_3 <- out$smooth_df_3
  smooth_df_4 <- out$smooth_df_4
  
  all_probs <- c(
    smooth_df_1$prob,
    smooth_df_2$prob,
    smooth_df_3$prob,
    smooth_df_4$prob
  )
  
  mid_val <- median(all_probs, na.rm = TRUE)
  lim_range <- range(all_probs, na.rm = TRUE)
  
  library(grid)
  
  common_fill_scale <- scale_fill_gradient2(
    low = "white",
    mid = "lightblue",
    high = "blue",
    midpoint = mid_val,
    limits = lim_range,
    name = "Prob",
    guide = guide_colorbar(
      direction = "horizontal",
      barwidth = unit(6, "cm"),
      barheight = unit(0.35, "cm"),
      title.position = "top"
    )
  )
  
  rank_label <- paste0(rank, ordinal_suffix(rank), " largest achiever")
  
  common_theme <- theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      plot.tag.position = c(0.18, 1.01),
      plot.tag = element_text(face = "bold", size = 14),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
  
  p1 <- ggplot(smooth_df_1, aes(N1, N2, fill = prob)) +
    geom_raster(interpolate = TRUE) +
    common_fill_scale +
    labs(
      x = expression(N[1]),
      y = expression(N[2]),
      tag = "   (A)"
    ) +
    common_theme
  
  p2 <- ggplot(smooth_df_2, aes(n1, n2, fill = prob)) +
    geom_raster(interpolate = TRUE) +
    common_fill_scale +
    labs(
      x = expression(n[1]),
      y = expression(n[2]),
      tag = "(B)"
    ) +
    common_theme
  
  p3 <- ggplot(smooth_df_3, aes(n2, N1, fill = prob)) +
    geom_raster(interpolate = TRUE) +
    common_fill_scale +
    labs(
      x = expression(n[2]),
      y = expression(N[1]),
      tag = "   (C)"
    ) +
    common_theme
  
  p4 <- ggplot(smooth_df_4, aes(n1, N2, fill = prob)) +
    geom_raster(interpolate = TRUE) +
    common_fill_scale +
    labs(
      x = expression(n[1]),
      y = expression(N[2]),
      tag = "   (D)"
    ) +
    common_theme
  
  p <- ((p1 + p2) / (p3 + p4) / patchwork::guide_area()) +
    patchwork::plot_layout(
      guides = "collect",
      heights = c(1, 1, 0.12)
    ) +
    patchwork::plot_annotation(
      title = paste0(
        "Sensitivity Analysis: ",
        rank_label,
        ", ",
        y_label,
        " Achievement Distribution"
      ),
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 14
        )
      )
    )
  return(p)
}

for (r in 1:4) {
  
  normal_plots[[r]] <- plot_one_rank(
    out = normal_results[[r]],
    rank = r,
    y_label = "Normal"
  )
  
  print(normal_plots[[r]])
  
  ggsave(
    filename = paste0("sensitivity_normal_rank_", r, ".png"),
    plot = normal_plots[[r]],
    width = 10,
    height = 8,
    dpi = 300
  )
}

gamma_plots <- list()

for (r in 1:4) {
  
  gamma_plots[[r]] <- plot_one_rank(
    out = gamma_results[[r]],
    rank = r,
    y_label = "Gamma(shape = 2, rate = 1)"
  )
  
  print(gamma_plots[[r]])
  
  ggsave(
    filename = paste0("sensitivity_gamma_rank_", r, ".png"),
    plot = gamma_plots[[r]],
    width = 10,
    height = 8,
    dpi = 300
  )
}
