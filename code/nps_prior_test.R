require(mvtnorm)
source(file.path("code", "nps_prior.R"))

# illustration of nps prior method
# testing on small fake data w/ noninformative prior
set.seed(0)
n <- 60

test_x <- rnorm(n)
test_X <- cbind(1, test_x)
test_beta <- c(5, 7)
test_y <- rnorm(n, test_X %*% test_beta , 3)

mu_bet <- rep(0, 2)
a <- 0
b <- 0
k0 <- n
V <- diag(1, 2)

# testing MC
nps_list <- nps_mc(4000, test_y, test_X, mu_bet, k0, V, a, b)
str(nps_list, max.level = 1)

plot(test_x, test_y, ylim = c(-15, 30))

abline(nps_list$beta_mean[1], nps_list$beta_mean[2])
points(test_x, nps_list$fitted_y_means, col = "red")
points(test_x, nps_list$fitted_y_upper, col = "blue")
points(test_x, nps_list$fitted_y_lower, col = "blue")

points(test_x, nps_list$post_pred_upper, col = "purple")
points(test_x, nps_list$post_pred_lower, col = "purple")

newx <- c(-2, -1.9, -1.8, 1.9, 2, 2.1)
newX <- cbind(1, newx)
pred <- predict(nps_list, newdata = newX)

points(newx, pred$fitted_y_means, pch = 17)
points(newx, pred$fitted_y_upper, col = "blue")
points(newx, pred$fitted_y_lower, col = "blue")

points(newx, pred$post_pred_upper, col = "purple")
points(newx, pred$post_pred_lower, col = "purple")
