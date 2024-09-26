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

