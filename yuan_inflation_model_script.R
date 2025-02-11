library(haven)
library(ggplot2)
library(patchwork)

dat <- read_dta("replicationdata.dta")

# following GPW, assume money stock decays by 10% per year
M <- dat$nominal

M[1] <- dat$nominal[1]
for (i in 2:96){
  M[i] <- 0.9*M[i-1] + dat$nominal[i]
}

m <- log(M)
p <- log(dat$cpi[1])
y <- log(dat$pop[1])
v <- mean(log(dat$cpi) + log(dat$pop) - log(M))

simulate <- function(beta, lambda, m, p0, y0, v){
  # simulate the model with
  # m : log money stock
  # p0 : starting value of log CPI
  # y0 : starting value of log population
  # v : a constant
  # sigma : zero
  N <- length(m)
  mt <- m
  
  pt <- yt <- rep(0, N)
  pt[1] <- p0
  yt[1] <- y0
  
  for (i in 1:(N-1)){
    pt[i+1] <- (1-lambda)*pt[i] + lambda*(v + mt[i] - yt[i])
    yt[i+1] <- (1-beta)*yt[i] + beta*(v + mt[i] - pt[i])
  }
  list(p=pt, y=yt)
}

LL <- function(beta, lambda, sigma, m, pt, yt, v){
  # log-likelihood for the model
  # beta, lambda, sigma : parameters
  # m : data log money stock series
  # pt : data log CPI series
  # yt : data log population series
  # v : a constant
  # likelihood = \prod (1/sqrt(2pi))*(1/sigma) exp(-(obs-data)^2/2sigma^2)
  # where the product is over 2n terms with n = length of data series
  sim <- simulate(beta, lambda, m, pt[1], yt[1], v)
  L <- -(pt - sim$p)^2 / (2*sigma^2) - 2*log(sigma) - (yt - sim$y)^2 / (2*sigma^2)
  sum(L)
}

# find a starting point for the Gibbs sampler
pt <- log(dat$cpi)
yt <- log(dat$pop)
s <- sd(c(diff(pt), diff(yt)), na.rm=T)

b <- l <- seq(0, 1, len=100)
LLm <- matrix(0, nr=100, nc=100)
for (i in 1:100){
  for (j in 1:100){
    LLm[i,j] <- LL(b[i], l[j], s, m, pt, yt, v)
  }
}
q <- which(LLm == max(LLm), arr.ind=T)

beta_start <- b[q[1]]
lambda_start <- l[q[2]]
sigma_start <- s

start_pars <- list(beta_start=beta_start, lambda_start=lambda_start, 
              sigma_start=sigma_start)
obs_data <- list(m=m, pt=pt, yt=yt)
control <- list(proposal_distribution_length=100, # number of likelihood evaluations per parameter
                n_its=1000) # number of sampling iterations
# to get a chain which mixes properly, it would be better to take a bigger
# value of n_its and thin the output

gibbs <- function(control, obs_data, start_pars){
  
  n_its <- control$n_its
  proposal_distribution_length <- control$proposal_distribution_length
  
  beta <- lambda <- sigma <- LLout <- rep(0, n_its)
  sim_py <- matrix(0, ncol=2*2, nrow=n_its)
  colnames(sim_py) <- paste0(rep(c("p", "y"), each=2), rep(as.character(1:2), 2))
  
  m <- obs_data$m
  pt <- obs_data$pt
  yt <- obs_data$yt
  N <- length(m)
  
  beta[1] <- start_pars$beta_start
  lambda[1] <- start_pars$lambda_start
  sigma[1] <- start_pars$sigma_start
  LLout[1] <- LL(beta[1], lambda[1], sigma[1], m, pt, yt, v)
  
  for (i in 2:n_its){
    # sample beta
    beta_prop <- seq(0 , 1 , len=proposal_distribution_length)
    pr <- rep(0, proposal_distribution_length)
    for (j in 1:proposal_distribution_length){
      pr[j] <- LL(beta_prop[j], lambda[i-1], sigma[i-1], m, pt, yt, v)
    }
    pr <- exp(pr - max(pr))
    beta_new <- sample(beta_prop, 1, prob= pr/sum(pr))
    
    # sample lambda, given beta
    lambda_prop <- seq(0 , 1 , len=proposal_distribution_length)
    pr <- rep(0, proposal_distribution_length)
    for (j in 1:proposal_distribution_length){
      pr[j] <- LL(beta_new, lambda_prop[j], sigma[i-1], m, pt, yt, v)
    }
    pr <- exp(pr - max(pr))
    lambda_new <- sample(lambda_prop, 1, prob= pr/sum(pr))
    
    #sample sigma, given beta & lambda
    sim_new <- simulate(beta_new, lambda_new, m, pt[1], yt[1], v)
    A <- sum( (pt - sim_new$p)^2 / 2 + (yt - sim_new$y)^2 / 2 )
    
    # sample from inv-gamma, assuming improper prior 1/sigma
    sigma2_new <- 1/rgamma(1, N - 0.5, rate=A)
    sigma_new <- sqrt(sigma2_new)
    
    lambda[i] <- lambda_new
    beta[i] <- beta_new
    sigma[i] <- sigma_new
    
    LLout[i] <- LL(beta[i], lambda[i], sigma[i], m, pt, yt, v)
    
    # simulate inflation in next period
    sim <- simulate(beta[i], lambda[i], m[N] + rep(max(diff(m)), 2), 
                    pt[N], yt[N], v)
    sim_py[i, 1:2] <- sim$p
    sim_py[i, 3:4] <- sim$y
  }
  list(lambda=lambda, beta=beta, sigma=sigma, LL=LLout, sim_py=sim_py)
}

out <- gibbs(control, obs_data, start_pars)

Msim <- matrix(0, 1000, ncol=96)
for (i in 1:1000){
  Msim[i, ] <- simulate(out$beta[i], out$lambda[i], log(M), log(dat$cpi[1]), log(dat$pop[1]),
                        mean(log(dat$cpi) + log(dat$pop) - log(M)))$p + rnorm(96, 0, out$sigma[i])
}

# 90% credible intervals
upper <- exp(apply(Msim, 2, function(x) quantile(x, 0.95)))
lower <- exp(apply(Msim, 2, function(x) quantile(x, 0.05)))
middle <- exp(apply(Msim, 2, function(x) quantile(x, 0.5)))

pd <- data.frame(year=(1260):(1260+95), lower=lower, middle=middle,
                 upper=upper, cpi=dat$cpi)

plot1 <- ggplot(data=pd, aes(year=year, estimate=middle)) + 
  geom_line(aes(x=year, y=middle, col="middle"), linewidth=1.1, linetype="dashed") + 
  geom_ribbon(aes(x=year, ymin=lower, ymax=upper), fill="blue", alpha=0.1) + 
  geom_line(aes(x=year, y=cpi, col="cpi"), linewidth=1.1) + 
  ylab("CPI") +
  scale_color_manual(name='',
                     values = c("cpi"="red", "middle"="blue"),
                     labels = c("CPI", "fitted CPI"),
                     guide = 'legend') +
  theme(legend.position = c(0.1, 0.9)) + 
  scale_x_continuous(minor_breaks = seq(1260, 1356, 10)) +
  scale_y_continuous(minor_breaks = seq(0, 150, 10))

bl <- data.frame(x = c(out$lambda, out$beta), y=factor(rep(c("lambda","beta"), each=1000), levels=c("lambda", "beta")))
plot2 <- ggplot(bl, aes(x, fill=y, col=y)) + geom_density(alpha=0.4, linetype=0) +
  scale_fill_manual( values = c("blue3","brown2")) + 
  xlab("beta and lambda") + theme(legend.title = element_blank()) + 
  scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(minor_breaks = seq(0, 2.5, 0.5)) +
  theme(legend.position="bottom")

inflation <- (exp(out$sim_py[,2]) - exp(out$sim_py[,1]))/exp(out$sim_py[,1]) *100
plot3 <- ggplot(data.frame(x=inflation), aes(x)) + 
  geom_histogram(binwidth=10) +
  scale_x_continuous(minor_breaks = seq(0, 120, 10)) +
  xlab("simulated out-of-sample inflation (%)") +
  geom_vline(xintercept = mean(inflation), col="red")

#ggsave(plot=plot1, "yuan_gibbs_results.png", units="in", height=5, width=5)
#ggsave(plot=plot2, "yuan_beta_lambda2.png", units="in", height=4, width=5)
#ggsave(plot=plot3, "yuan_inflation.png", units="in", height=5, width=5)

## inference test
test_inference <- function(beta, lambda, sigma, pt, yt, m){
  N <- length(m)
  sim <- simulate(beta, lambda, m, pt[1], yt[1], v)
  pt <- sim$p + rnorm(N, 0, sigma)
  yt <- sim$y + rnorm(N, 0, sigma)
  obs_data <- list(m=m, pt=pt, yt=yt)
  
  s <- sd(c(diff(pt), diff(yt)), na.rm=T)
  
  b <- l <- seq(0, 1, len=100)
  LLm <- matrix(0, nr=100, nc=100)
  for (i in 1:100){
    for (j in 1:100){
      LLm[i,j] <- LL(b[i], l[j], s, m, pt, yt, v)
    }
  }
  q <- which(LLm == max(LLm), arr.ind=T)
  
  beta_start <- b[q[1]]
  lambda_start <- l[q[2]]
  sigma_start <- s
  start_pars <- list(beta_start=beta_start, lambda_start=lambda_start, 
                     sigma_start=sigma_start)
  
  g <- gibbs(control, obs_data, start_pars)
  df <- data.frame(beta_est = g$beta, lambda_est=g$lambda, sigma_est=g$sigma)[100:1000, ]
  p1 <- ggplot(df, aes(beta_est)) + geom_histogram() + 
    geom_vline(xintercept = mean(df$beta_est), col="red") +
    geom_vline(xintercept = beta, col="green")
  p2 <- ggplot(df, aes(lambda_est)) + geom_histogram() + 
    geom_vline(xintercept = mean(df$lambda_est), col="red") +
    geom_vline(xintercept = lambda, col="green")
  p3 <- ggplot(df, aes(sigma_est)) + geom_histogram() + 
    geom_vline(xintercept = mean(df$sigma_est), col="red") +
    geom_vline(xintercept = sigma, col="green")
  
  p1 + p2 / p3
}
set.seed(101)
test_inference(0.3, 0.4, 0.5, pt, yt, m)