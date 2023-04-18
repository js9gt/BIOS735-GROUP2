library(devtools)
document("ZINB")
load_all("ZINB")

build("ZINB")

# test the package
n <- 100
X <- matrix(c(rep(1, 50), rbinom(50, 1, 0.5)), 50, 2)
y <- sample(c(rnbinom(20, 10, 0.5), rep(0, 30)))
fit.mod <- fit.zinb(X, y, add_intercept = FALSE)
fit.mod

result <- LRT1D(X, y, add_intercept=FALSE)
result

Y <- matrix(c(y, sample(y)), 50, 2)
result <- LRTnD(X, Y, add_intercept=FALSE)
result