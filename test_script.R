

test <- function(x){
  pnorm(x, lower = TRUE) - pnorm(x, mean = 1, sd = 1, lower = TRUE)
}
y <- seq(from = -4, to = 6, length = 20000)
z <-  sapply(y, test)
plot(y, z)
length(which(z < 0))
y[which.max(z)]
#plot(density(replicate(100000, rnorm(1, mean = 0.5, sd = 1/sqrt(2)) - rnorm(1, sd = 1/sqrt(2)))), 
#     ylim = c(0, 0.8))
plot(density(replicate(200000, rnorm(1, mean = 0.5, sd = 1))), 
     ylim = c(0, 0.8))
lines(y, z)


test <- function(x){
  (pnorm(x, sd = 2, lower = TRUE) - pnorm(x, mean = 4, sd = 2, lower = TRUE)) / 4
}
y <- seq(from = -8, to = 10, length = 20000)
z <-  sapply(y, test)
length(which(z < 0))
y[which.max(z)]

plot(y, sapply(y, function(x){
  dnorm(x, mean = 2, 
        sd = 1/(sqrt(2*pi) * test(2)) )
}),
ylim = c(0, 0.22)
)                   
lines(y, z, lty = 2, col = "red", lwd = 2)


test <- function(x){
  (pnorm(x, sd = 3, lower = TRUE) - pnorm(x, mean = 6, sd = 3, lower = TRUE)) / 6
}
y <- seq(from = -9, to = 15, length = 20000)
z <-sapply(y, test)
plot(y, z)
length(which(z < 0))
y[which.max(z)]

plot(y, sapply(y, function(x){
  dnorm(x, mean = 3, 
        sd = 1/(sqrt(2*pi) * test(3)) )
}),
ylim = c(0, 0.15)
)                   
lines(y, z, lty = 2, col = "red", lwd = 2)


integrate(test, lower = -20, upper = 20)
optim(0, fn = function(x) - test(x), method = "Brent", lower = -20, upper = 20)
