##### Q1 #####

library(geoR)
library(rjags)
library(ggplot2)

### Preparing the data
data("gambia")
Y <- gambia$pos

X <- as.matrix(subset(x = gambia, select = -c(pos)))
X <- scale(X)

##### Q 4-8 a) #####

# I have saved the samples and the model as an
# .Rdata file, as takes long to run the MCMC 
# chain again

load('./out-4A.Rdata')
model <- outA$model
samples <- outA$samples

### JAGS model
modelString <- textConnection("model{
                              
  # Likelihood
  for(i in 1:n)
  {
    Y[i] ~ dbern(q[i])
    logit(q[i]) = inprod(X[i, ], beta[])
  }
  
  # Priors
  for(j in 1:p)
  {
    beta[j] ~ dnorm(0, 100^(-2))
  }

}")

data.list <- list(Y = Y, X = X,
                  n = length(Y),
                  p = length(X[1, ]))
model <- jags.model(file = modelString,
                    data = data.list,
                    n.chains = 2)

# Burn-in
update(model, 1e3)

# Samples

n.samples <- 1e4
params <- c("beta")

samples <- coda.samples(model = model,
                        variable.names = params,
                        n.iter = n.samples)
plot(samples)
summary(samples)

# outA <- list(model = model,
#              samples = samples)
# save(outA, file = c('out-4A.Rdata'))

##### Q 4-8 b) #####

# I have saved the samples and the model as an
# .Rdata file, as takes long to run the MCMC 
# chain again

load('./out-4B.Rdata')
model <- outB$model
samples <- outB$samples

### Grouping the data

unique.coords <- unique(x = cbind(gambia$x, gambia$y))
L <- length(unique.coords[,1])
s <- numeric(length = length(Y))

for (i in 1:L) 
{
  coords <- unique.coords[i, ]
  index <- (gambia$x == coords[1] & gambia$y == coords[2])
  s[index] <- i
}

### JAGS code

### JAGS model
modelString <- textConnection("model{
                              
  # Likelihood
  for(i in 1:n)
  {
    Y[i] ~ dbern(q[i])
    logit(q[i]) = inprod(X[i, ], beta[]) + alpha[s[i]]
  }
  
  # Priors
  for(j in 1:L)
  {
    alpha[j] ~ dnorm(0, alpha.precision)
  }
  alpha.precision ~ dgamma(0.01, 0.01) 
  tau.sq = 1/alpha.precision
  
  for(j in 1:p)
  {
    beta[j] ~ dnorm(0, 100^(-2))
  }
  
}")

data.list <- list(Y = Y, 
                  X = X,
                  s = s,
                  n = length(Y),
                  p = length(X[1,]),
                  L = L)
model <- jags.model(file = modelString,
                    data = data.list,
                    n.chains = 2)

# Burn-in
update(model, 1e3)

# Samples
n.samples <- 1e4
params <- c("beta", "alpha")

samples <- coda.samples(model = model,
                        variable.names = params,
                        n.iter = n.samples)
plot(samples)
summary(samples)

# outB <- list(model = model,
#              samples = samples)
# save(outB, file = c('out-4B.Rdata'))

samples.mat <- as.matrix(samples[[1]][,1:L])
samples.mean <- apply(samples.mat, 2, mean)

par(mfrow = c(1,1))

plot.data <- data.frame(X = unique.coords[ ,1], 
                        Y = unique.coords[ ,2],
                        alpha = samples.mean)
plt <- ggplot(plot.data, aes(X, Y))+
  geom_point(aes(color = alpha)) +
  scale_color_gradient(low = "green", high = "red") +
  labs(title = "Spatial plot of mean(alpha)")
plt





##### Q2 #####

library(MASS)
library(rjags)

data(galaxies)
Y <- galaxies
hist(Y,breaks=25)

# JAGS code
modelString <- textConnection("model{
                              
  # Likelihood
  for(i in 1:n)
  {
    Y[i] ~ dnorm(mu[i], tau[i])
    mu[i] = mu.cluster[z[i]]
    tau[i] = tau.cluster[z[i]]
    z[i] ~ dcat(pi[1:K])
  }
  
  # Priors
  for(j in 1:K)
  {
    mu.cluster[j] ~ dnorm(0, 10^(-10))
    tau.cluster[j] ~ dgamma(0.01, 0.01)
    sig2.cluster[j] = 1/tau.cluster[j]
  }     
  pi[1:K] ~ ddirch(c(1, 1, 1))
  
  for(g in 1:G)
  {
    y[g] = pi[1] * dnorm(y.grid[g], mu.cluster[1], tau.cluster[1]) + 
    pi[2] * dnorm(y.grid[g], mu.cluster[2], tau.cluster[2]) +
    pi[3] * dnorm(y.grid[g], mu.cluster[3], tau.cluster[3])
  }
  
}")

# The required grid
y.grid <- seq(from = 5000, to = 40000, by = 100)
G <- length(y.grid)

data.list <- list(Y = Y,
                  n = length(Y),
                  K = 3,
                  y.grid = y.grid,
                  G = G)

model <- jags.model(file = modelString, 
                    data = data.list, 
                    n.chains = 1)
update(model, n.iter = 1e4)

# generating samples
params <- c("mu.cluster", "sig2.cluster", "pi", "y")
samples <- coda.samples(model = model,
                        variable.names = params, 
                        n.iter = 1e4)
# plot(samples[[1]][, 1:9])

y.density.samples <- as.matrix(samples[[1]][ ,10:(G+10-1)])
y.density.median <- apply(X = y.density.samples, 
                          MARGIN = 2, FUN = median)
y.density.lower <- apply(X = y.density.samples, 
                         MARGIN = 2, 
                         FUN = function(x){quantile(x , p = 0.025)})
y.density.upper <- apply(X = y.density.samples, 
                         MARGIN = 2, 
                         FUN = function(x){quantile(x , p = 0.975)})

hist(Y, probability = T, breaks = 25,
     xlab = "Y", ylab = 'Density',
     main = "Gaussian Mixture Model (K = 3)")
lines(x = y.grid, y = y.density.upper,
      type = 'l', col = 'red', lwd = 1)
lines(x = y.grid, y = y.density.lower,
      col = 'blue', lwd = 1)
lines(x = y.grid, y = y.density.median,
      col = 'black', lwd = 2)

legend("topright", legend = c("Median", "2.5%", "97.5%"),
       lty = 1, col = c("black", "blue", "red"), lwd = 2)




##### Q3 #####
library(rjags)

Y1 <- 563
N1 <- 2820

Y2 <- 10
N2 <- 27

##### BF #####
c <- 1

num.gamma1 <- pgamma(q = 1, shape = Y1+1, rate = N1)
num.gamma2 <- pgamma(q = 1, shape = Y2+1, rate = N2)
num.gamma <- num.gamma1 * num.gamma2

denom.gamma <- pgamma(q = 1, shape = Y1+Y2+1, rate = N1+N2)

numer <- num.gamma * (1 + N2/N1)^(Y1+Y2+1)
denom <- denom.gamma * c * N2 * (N2/N1)^Y2 * choose(Y1+Y2, Y1)

BF1.12 <- numer / denom
print(BF1.12)

c <- 10

num.gamma1 <- pgamma(q = 1, shape = Y1+1, rate = N1)
num.gamma2 <- pgamma(q = 1, shape = Y2+1, rate = N2)
num.gamma <- num.gamma1 * num.gamma2

denom.gamma <- pgamma(q = 1, shape = Y1+Y2+1, rate = N1+N2)

numer <- num.gamma * (1 + N2/N1)^(Y1+Y2+1)
denom <- denom.gamma * c * N2 * (N2/N1)^Y2 * choose(Y1+Y2, Y1)

BF10.12 <- numer / denom
print(BF10.12)

##### JAGS models (c = 1) #####
c <- 1
data.list <- list(Y1 = Y1, N1 = N1,
                  Y2 = Y2, N2 = N2, 
                  c = c)

modelString1 <- textConnection("model{

  # Likelihood
  Y1 ~ dpois(N1 * lambda1)
  Y2 ~ dpois(N2 * lambda2)
  
  # Prior
  lambda1 ~ dunif(0, c)
  lambda2 ~ dunif(0, c)

}")

modelString2 <- textConnection("model{

  # Likelihood
  Y1 ~ dpois(N1 * lambda0)
  Y2 ~ dpois(N2 * lambda0)
  
  # Prior
  lambda0 ~ dunif(0, c)

}")

### Model 1
model1 <- jags.model(file = modelString1,
                     data = data.list,
                     n.chains = 2)
update(model1, n.iter = 1e4)
params <- c("lambda1", "lambda2")
samples1 <- coda.samples(model = model1,
                         n.iter = 1e4,
                         variable.names = params)

plot(samples1)
lambda1 <- samples1[[1]][ ,1]
lambda2 <- samples1[[1]][ ,2]

### Model 2
model2 <- jags.model(file = modelString2,
                     data = data.list,
                     n.chains = 2)
update(model2, n.iter = 1e4)
params <- c("lambda0")
samples2 <- coda.samples(model = model2,
                         n.iter = 1e4,
                         variable.names = params)

plot(samples2)
lambda0 <- samples2[[1]][ ,1]

##### DIC #####
loglike.m1.func <- function(iter)
{
  iter.lambda1 <- lambda1[iter]
  iter.lambda2 <- lambda2[iter]
  
  loglike.Y1 <- dpois(x = Y1, N1 * iter.lambda1, log = T)
  loglike.Y2 <- dpois(x = Y2, N2 * iter.lambda2, log = T)
  
  return(loglike.Y1 + loglike.Y2)
}

loglike.m2.func <- function(iter)
{
  iter.lambda <- lambda0[iter]
  
  loglike.Y1 <- dpois(x = Y1, N1 * iter.lambda, log = T)
  loglike.Y2 <- dpois(x = Y2, N2 * iter.lambda, log = T)
  
  return(loglike.Y1 + loglike.Y2)
}

loglike.m1 <- sapply(1:1e4, FUN = loglike.m1.func)
loglike.m2 <- sapply(1:1e4, FUN = loglike.m2.func)

deviance.m1 <- -2 * loglike.m1
deviance.m2 <- -2 * loglike.m2

Dbar.m1 <- mean(deviance.m1) 
Dbar.m2 <- mean(deviance.m2)

loglike.thetahat.m1 <- dpois(Y1, lambda = N1 * mean(lambda1), log = T) + 
  dpois(Y2, lambda = N2 * mean(lambda2), log = T)
D.thetahat.m1 <- -2 * loglike.thetahat.m1

loglike.thetahat.m2 <- dpois(Y1, lambda = N1 * mean(lambda0), log = T) + 
  dpois(Y2, lambda = N2 * mean(lambda0), log = T)
D.thetahat.m2 <- -2 * loglike.thetahat.m2

pD.m1 <- Dbar.m1 - D.thetahat.m1 
pD.m2 <- Dbar.m2 - D.thetahat.m2

DIC.m1 <- pD.m1 + Dbar.m1 
DIC.m2 <- pD.m2 + Dbar.m2

DIC.c1 <- c(DIC.m1, DIC.m2)

##### WAIC #####
loglike.m1.func <- function(iter)
{
  iter.lambda1 <- lambda1[iter]
  iter.lambda2 <- lambda2[iter]
  
  loglike.Y1 <- dpois(x = Y1, N1 * iter.lambda1, log = T)
  loglike.Y2 <- dpois(x = Y2, N2 * iter.lambda2, log = T)
  
  return(c(loglike.Y1, loglike.Y2))
}

loglike.m2.func <- function(iter)
{
  iter.lambda <- lambda0[iter]
  
  loglike.Y1 <- dpois(x = Y1, N1 * iter.lambda, log = T)
  loglike.Y2 <- dpois(x = Y2, N2 * iter.lambda, log = T)
  
  return(c(loglike.Y1, loglike.Y2))
}

loglike.m1 <- sapply(1:1e4, FUN = loglike.m1.func)
loglike.m2 <- sapply(1:1e4, FUN = loglike.m2.func)

posmeans.m1 <- apply(loglike.m1, 1, mean)
posmeans.m2 <- apply(loglike.m2, 1, mean)

posvars.m1 <- apply(loglike.m1, 1, var)
posvars.m2 <- apply(loglike.m2, 1, var)

pW.m1 <- sum(posvars.m1) 
pW.m2 <- sum(posvars.m2)

sum.means.m1 <- sum(posmeans.m1) 
sum.means.m2 <- sum(posmeans.m2)

WAIC.m1 <- -2 * sum.means.m1 + 2 * pW.m1
WAIC.m2 <- -2 * sum.means.m2 + 2 * pW.m2

WAIC.c1 <- c(WAIC.m1, WAIC.m2)








##### JAGS models (c = 10) #####
c <- 10
data.list <- list(Y1 = Y1, N1 = N1,
                  Y2 = Y2, N2 = N2, 
                  c = c)

modelString1 <- textConnection("model{

  # Likelihood
  Y1 ~ dpois(N1 * lambda1)
  Y2 ~ dpois(N2 * lambda2)
  
  # Prior
  lambda1 ~ dunif(0, c)
  lambda2 ~ dunif(0, c)

}")

modelString2 <- textConnection("model{

  # Likelihood
  Y1 ~ dpois(N1 * lambda0)
  Y2 ~ dpois(N2 * lambda0)
  
  # Prior
  lambda0 ~ dunif(0, c)

}")

### Model 1
model1 <- jags.model(file = modelString1,
                     data = data.list,
                     n.chains = 2)
update(model1, n.iter = 1e4)
params <- c("lambda1", "lambda2")
samples1 <- coda.samples(model = model1,
                         n.iter = 1e4,
                         variable.names = params)

plot(samples1)
lambda1 <- samples1[[1]][ ,1]
lambda2 <- samples1[[1]][ ,2]

### Model 2
model2 <- jags.model(file = modelString2,
                     data = data.list,
                     n.chains = 2)
update(model2, n.iter = 1e4)
params <- c("lambda0")
samples2 <- coda.samples(model = model2,
                         n.iter = 1e4,
                         variable.names = params)

plot(samples2)
lambda0 <- samples2[[1]][ ,1]

##### DIC #####
loglike.m1.func <- function(iter)
{
  iter.lambda1 <- lambda1[iter]
  iter.lambda2 <- lambda2[iter]
  
  loglike.Y1 <- dpois(x = Y1, N1 * iter.lambda1, log = T)
  loglike.Y2 <- dpois(x = Y2, N2 * iter.lambda2, log = T)
  
  return(loglike.Y1 + loglike.Y2)
}

loglike.m2.func <- function(iter)
{
  iter.lambda <- lambda0[iter]
  
  loglike.Y1 <- dpois(x = Y1, N1 * iter.lambda, log = T)
  loglike.Y2 <- dpois(x = Y2, N2 * iter.lambda, log = T)
  
  return(loglike.Y1 + loglike.Y2)
}

loglike.m1 <- sapply(1:1e4, FUN = loglike.m1.func)
loglike.m2 <- sapply(1:1e4, FUN = loglike.m2.func)

deviance.m1 <- -2 * loglike.m1
deviance.m2 <- -2 * loglike.m2

Dbar.m1 <- mean(deviance.m1) 
Dbar.m2 <- mean(deviance.m2)

loglike.thetahat.m1 <- dpois(Y1, lambda = N1 * mean(lambda1), log = T) + 
  dpois(Y2, lambda = N2 * mean(lambda2), log = T)
D.thetahat.m1 <- -2 * loglike.thetahat.m1

loglike.thetahat.m2 <- dpois(Y1, lambda = N1 * mean(lambda0), log = T) + 
  dpois(Y2, lambda = N2 * mean(lambda0), log = T)
D.thetahat.m2 <- -2 * loglike.thetahat.m2

pD.m1 <- Dbar.m1 - D.thetahat.m1 
pD.m2 <- Dbar.m2 - D.thetahat.m2

DIC.m1 <- pD.m1 + Dbar.m1 
DIC.m2 <- pD.m2 + Dbar.m2

DIC.c10 <- c(DIC.m1, DIC.m2)

##### WAIC #####
loglike.m1.func <- function(iter)
{
  iter.lambda1 <- lambda1[iter]
  iter.lambda2 <- lambda2[iter]
  
  loglike.Y1 <- dpois(x = Y1, N1 * iter.lambda1, log = T)
  loglike.Y2 <- dpois(x = Y2, N2 * iter.lambda2, log = T)
  
  return(c(loglike.Y1, loglike.Y2))
}

loglike.m2.func <- function(iter)
{
  iter.lambda <- lambda0[iter]
  
  loglike.Y1 <- dpois(x = Y1, N1 * iter.lambda, log = T)
  loglike.Y2 <- dpois(x = Y2, N2 * iter.lambda, log = T)
  
  return(c(loglike.Y1, loglike.Y2))
}

loglike.m1 <- sapply(1:1e4, FUN = loglike.m1.func)
loglike.m2 <- sapply(1:1e4, FUN = loglike.m2.func)

posmeans.m1 <- apply(loglike.m1, 1, mean)
posmeans.m2 <- apply(loglike.m2, 1, mean)

posvars.m1 <- apply(loglike.m1, 1, var)
posvars.m2 <- apply(loglike.m2, 1, var)

pW.m1 <- sum(posvars.m1) 
pW.m2 <- sum(posvars.m2)

sum.means.m1 <- sum(posmeans.m1) 
sum.means.m2 <- sum(posmeans.m2)

WAIC.m1 <- -2 * sum.means.m1 + 2 * pW.m1
WAIC.m2 <- -2 * sum.means.m2 + 2 * pW.m2

WAIC.c10 <- c(WAIC.m1, WAIC.m2)





##### Q4 #####

library(geoR)
library(rjags)

### Preparing the data
data("gambia")
Y <- gambia$pos

# Added column of 1's for beta_0
X <- as.matrix(subset(x = gambia, select = -c(x, y, pos)))
X <- scale(X)
X <- cbind(1, X)

### JAGS model
modelString <- textConnection("model{
                              
  # Likelihood
  for(i in 1:n)
  {
    Y[i] ~ dbern(q[i])
    logit(q[i]) = inprod(X[i, ], beta[])
  }
  
  # Priors
  for(j in 1:p)
  {
    beta[j] ~ dnorm(0, 100^(-2))
  }
  
  # PPD
  for(i in 1:n)
  {
    Yppd[i] ~ dbern(t[i])
    logit(t[i]) = inprod(X[i, ], beta[])
  }
  
  D[1] <- mean(Yppd[])
  D[2] <- sd(Yppd[])

}")

data.list <- list(Y = Y, X = X,
                  n = length(Y),
                  p = length(X[1, ]))
model <- jags.model(file = modelString,
                    data = data.list,
                    n.chains = 2)

# Burn-in
update(model, 5e3)

# Samples

n.samples <- 1e4
params <- c("D")

samples <- coda.samples(model = model,
                        variable.names = params,
                        n.iter = n.samples)

# out <- list(model = model,
#             samples = samples)
# save(out, file = 'out-6.Rdata')

load('./out-6.Rdata')
model <- out$model
samples <- out$samples
plot(samples)

ppd.mean <- samples[[1]][ ,1]
ppd.sd <- samples[[1]][ ,2]

plot(density(ppd.mean), main = "Mean Y", ylim = c(0, 34))
abline(v = mean(Y), col = 'red')
legend("topleft", legend = c("Mean (PPD samples)", "Mean (Data)"),
       lty = 1, col = c('black', 'red'))
pval1 <- mean(ppd.mean > mean(Y))

plot(density(ppd.sd), main = "SD Y")
abline(v = sd(Y), col = 'red')
legend("topleft", legend = c("SD (PPD samples)", "SD (Data)"),
       lty = 1, col = c('black', 'red'))
pval2 <- mean(ppd.sd > sd(Y))

p.table <- matrix(data = c(pval1, pval2), nrow = 1)
colnames(p.table) <- c("Mean", "SD")

print(p.table)





##### Q5 #####

library(datasets)
library(rjags)

data("WWWusage")
Y <- as.numeric(WWWusage)

plot(WWWusage, type = 'l', lwd = 2, 
     ylab = "Y", main = "Data (Time series)")

# JAGS model
modelString <- "model{

  # Likelihood
  for(i in 5:n)
  {
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta[1] + sum(Y[(i-L):(i-1)] * beta[2:(L+1)])
  }

  # Prior
  for(j in 1:(L+1))
  {
    beta[j] ~ dnorm(0, 100^(-2))
  }
  tau ~ dgamma(0.01, 0.01)
  sig2 = 1/tau
  
}"

# L = 1
data.list1 <- list(Y = Y, n = length(Y), L = 1)
model1 <- jags.model(file = textConnection(modelString),
                     data = data.list1,
                     n.chains = 2)
update(model1, n.iter = 1e3)
params <- c("beta", "sig2")
samples1 <- coda.samples(model = model1,
                         variable.names = params,
                         n.iter = 1e4,
                         thin = 2)
plot(samples1)

# L = 2
data.list2 <- list(Y = Y, n = length(Y), L = 2)
model2 <- jags.model(file = textConnection(modelString),
                     data = data.list2,
                     n.chains = 2, quiet = T)
update(model2, n.iter = 1e3)
params <- c("beta", "sig2")
samples2 <- coda.samples(model = model2,
                         variable.names = params,
                         n.iter = 1e4)
plot(samples2)
# dic.samples(model = model2, n.iter = 1e4)
dic2 <- 511.8

# L = 3
data.list3 <- list(Y = Y, n = length(Y), L = 3)
model3 <- jags.model(file = textConnection(modelString),
                     data = data.list3,
                     n.chains = 2, quiet = T)
update(model3, n.iter = 1e3)
params <- c("beta", "sig2")
samples3 <- coda.samples(model = model3,
                         variable.names = params,
                         n.iter = 1e4,
                         thin = 1)
plot(samples3)
# dic.samples(model = model3, n.iter = 1e4)
dic3 <- 506.5

# L = 4
data.list4 <- list(Y = Y, n = length(Y), L = 4)
model4 <- jags.model(file = textConnection(modelString),
                     data = data.list4,
                     n.chains = 2, quiet = T)
update(model4, n.iter = 1e3)
params <- c("beta", "sig2")
samples4 <- coda.samples(model = model4,
                         variable.names = params,
                         n.iter = 1e4,
                         thin = 1)
plot(samples4)
# dic.samples(model = model4, n.iter = 1e4)
dic4 <- 498.9
