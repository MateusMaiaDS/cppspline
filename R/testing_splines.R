# Generating some random data
rm(list=ls())
library(purrr)
library(tidyverse)
Rcpp::sourceCpp("src/code.cpp")
n <- 100
x <- (runif(n = n,min = -5,max = 5)) %>% sort
x_new <- seq(-4,4,length.out = 500)
y <- sin(x) + rnorm(n = n, sd = 0.1)


n_splines_ <- 10
n_post_ <- 1500
tau_b_ <- n_splines_*0.1

sampler_obj <- sum_spline_mcmc_sampler(x = x,x_new = x,n_splines = n_splines_,
                                       y = y,n_post = n_post_,n_burn = 500,tau = 1,
                                       tau_b = tau_b_)

# plot(x,y,xlim = c(-6,6),  ylim = c(-3,3))
# points(x,sampler_obj[[1]] %>% colMeans(), col = "red")
#
# for(i in 1:n_splines_){
#      points(x,sampler_obj[[5]][,i], col =(i+1),pch=20)
# }


# Each tree prediction
tree_predictions <- matrix(0, nrow = length(x), ncol = n_splines_)

for(i in 1:n_post_){
     # tree_predictions[,i] <- sampler_obj[[6]][,i]
     tree_predictions <- tree_predictions +  sampler_obj[[6]][,,i]

}
tree_predictions <- tree_predictions/n_post_
colnames(tree_predictions) <- paste0("spline.",1:n_splines_)
tree_predictions <- cbind(x,tree_predictions)
# Plotting a ggplot
tree_predictions_tidy <- tree_predictions %>% as.data.frame() %>% pivot_longer(starts_with("spline"))

ggplot()+
     geom_point(data = data.frame(x = x, y = y), mapping = aes(x = x, y = y))+
     geom_line(data = data.frame(x = x, y = sampler_obj[[1]] %>% colMeans()), mapping = aes(x = x, y = y), col = "blue")+
     geom_line(data = tree_predictions_tidy, mapping = aes(x = x, y = value, col = name), alpha = 0.5) +
     # geom_line(data = data.frame(x = x, y = sin(x)), mapping = aes(x = x, y = y ), col = "orange" )+
     ylim(c(-2,2))+
     ggtitle(paste0("Tau_b is: ", tau_b_))+
     theme_bw()

# plot(x,y)
# points(x,sampler_obj[[2]] %>% colMeans(), col = "blue")
# low_ci <- sampler_obj[[1]] %>% apply(2,function(x){quantile(x,probs = c(0.025))})
# up_ci <- sampler_obj[[1]] %>% apply(2,function(x){quantile(x,probs = c(0.975))})
# points(x,low_ci, col = "red", pch = 20)
# points(x,up_ci, col = "red", pch = 20)
# points(x_new,sin(x_new), col = "green")
# plot(sampler_obj[[4]], type = "l")

# ============
# Analsing the variance
# =============
# B <- bspline(x = x,x_obs = x)
# btb <- tcrossprod(B)
# plot(diag(btb), type = "l")
#
# plot(diag(crossprod(B,(btb + (tau_b_/100)*diag(nrow = nrow(btb))))%*%B))
#
# # Creatinga X covariate
# x_test <- runif(n = 100) %>% matrix(ncol = 2)
# crossprod(x_test) %>% diag
