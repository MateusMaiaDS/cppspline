# Generating some random data
library(purrr)
Rcpp::sourceCpp("src/code.cpp")
n <- 100
x <- runif(n = n,min = -5,max = 5) %>% sort
x_new <- seq(-4,4,length.out = 500)
y <- sin(x) + rnorm(n = n, sd = 0.1)

sampler_obj <- mcmc_sampler(x = x,x_new = x,y = y,n_post = 1500,n_burn = 500,tau = 1)

plot(x,y,xlim = c(-6,6),  ylim = c(-2,2))
points(x,sampler_obj[[1]] %>% colMeans(), col = "red")
# points(x_new,sampler_obj[[2]] %>% colMeans(), col = "blue")
# low_ci <- sampler_obj[[1]] %>% apply(2,function(x){quantile(x,probs = c(0.025))})
# up_ci <- sampler_obj[[1]] %>% apply(2,function(x){quantile(x,probs = c(0.975))})
# points(x,low_ci, col = "red", pch = 20)
# points(x,up_ci, col = "red", pch = 20)
points(x_new,sin(x_new), col = "green")
plot(sampler_obj[[4]], type = "l")

