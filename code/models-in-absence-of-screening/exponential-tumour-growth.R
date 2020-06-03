### OBS FUNCTIONS IN THE END ###

### Libraries ###
library("pracma") # mod, size functions
library("tseries") # ACF plots for MCMC diagnostics


############# DATA GENERATION (synthetic data) ####### 
set.seed(1)    
n <- 10000   # Number of observations  
tau <- exp(-0.165)           
gammavar <- exp(-9.602)     
y <- M_exp(n, tau, gammavar)
# Summary stats: counts in each category
sobs <- summary_stats(y)


################# ML-estimates ####################
set.seed(1)
start_time <- Sys.time()
n <- 10000
tau <- exp(-0.165)
gammavar <- exp(-9.602)
w <- c(tau, gammavar)
y <- M_exp(n, w[1], w[2])
opt <- optim(w, Ln_exp, gr=NULL,
             vcd=y, negative=TRUE, method="L-BFGS-B",
             lower=c(0.01,0.000001), upper=c(2,1))
ML_time <- Sys.time() - start_time
print(ML_time)
print(opt$par)


##################### ABC-MCMC  ######################
# Euclidean distance metric weighted by sd estimates  
# Tolerance tempering
# Normal transition kernel with sd = rho 

# V[s] estimates #
set.seed(1)
n_mad <- 1000 
smad <- numeric(24)
tmp <- matrix(0,n_mad,24)
for(i in 1:n_mad){
  y <- M_exp(n, tau, gammavar)
  tmp[i,] <- summary_stats(y)
}
smad <- apply(tmp, 2, sd)

# Hyperparameters (adjust these)
M <- 10000                    # Number of iterations 
trans_cv <- c(0.1, 0.1)      # transition kernel var
w <- 1.5*c(tau, gammavar)
epsilon <- 1.0                # Tolerance 
quantile <- 90                # Quantile for tempering     
final_acc_rate <- 0.01        # Final accetance rate

# Initialization (do not change)
set.seed(0)
start_time <- Sys.time()
d <- length(w)                   # Number of parameters 
n_stat <- length(sobs)           # Number of statistics 
w_post <- matrix(NaN,M,d)        # stores postrior values 
ssim <- numeric(n_stat)          
ssim_mat <- matrix(NaN,M,n_stat) # Stores simulated stats
sdist_vec <- rep(NaN, M)           
n_accABC <- 0                    
acc_vec <- numeric(M)             
acc_rate <- 1                    
rho <- trans_cv*w       # sd of transition kernel

for(iter in 1:M){
  
  if(mod(iter,1000)==0){
    print(iter)
    acc_rate <- mean(acc_vec[(iter-999):(iter-1)])
  } 
  if(mod(iter,500)==0){
    par(mfrow=c(d,1))
    for(i in 1:d){
      plot(w_post[,i],type='l')
      lines(w_post[,i])
    }
  }  
  
  # Sample from the transition kernel 
  wp <- rnorm(d, mean=w, sd=rho)
  if(any(wp<=0)){
    w_post[iter,] <- w
    acc_vec[iter] <- 0
    next     
  }  
  
  # Simulate from the MODEL using the proposed parameter values 
  y_sim <- M_exp(n, wp[1], wp[2])
  
  # Calculate sufficient statistics from the samples
  ssim <- summary_stats(y_sim)
  
  # Caclulate normalized euclidian distence between observed and simulated statistics 
  sdist <- sqrt(sum(((ssim-sobs)/smad)^2))/n_stat 
  
  # If the distance between the data statistics and the simulated statistics 
  # is bigger then epsilon, jump to next iteration
  if(sdist >= epsilon){
    w_post[iter,] <- w
    acc_vec[iter] <- 0
    next 
  }
  
  # Epsilon is set to the q:th quantile of accepted errors 
  # every 100th iteration   
  n_accABC <- n_accABC + 1
  ssim_mat[iter,] <- ssim 
  sdist_vec[n_accABC] <- sdist
  if(mod(n_accABC, 100) == 0 && acc_rate > final_acc_rate){
    sdist_sorted <- sort(sdist_vec[(n_accABC-99):n_accABC])
    epsilon <- sdist_sorted[quantile]
    cat(iter, "", epsilon, " \n") #Print the epsilon values
  }
  
  prior_ratio <- 1
  
  alpha <- min(1, prior_ratio) # MH acceptence probability 
  
  U <- runif(1) 
  if(alpha > U){
    w_post[iter,] <- wp
    w <- wp 
    acc_vec[iter] <- 1
  } else {
    w_post[iter,] <- w
    acc_vec[iter] <- 0
  }  
}
w_post_ABC <- w_post
abc_time <- Sys.time() - start_time
print("ABC-MCMC time")
print(abc_time)

# MIXING PLOTS
par(mfrow=c(d,1))
plot(w_post[,1],type='l', ylab=expression(hat(tau)))
lines(w_post[,1])
plot(w_post[,2],type='l', ylab=expression(hat(eta)))
lines(w_post[,2])

# POSTERIOR MEANS #
burn_in <- 6000              
w_mean <- colMeans(w_post[burn_in:M,])
print("ABC-MCMC acceptance rate:")
print(mean(acc_vec[burn_in:M]))
print("ABC-MCMC means:")
print(w_mean)
print("ABC-MCMC mean relative error in %")
print((w_mean-c(tau,gammavar))/c(tau,gammavar)*100)

# LINEAR REGRESSION ADJUSTMENT  #
w_post_adj <- vector(mode="list", length=d)
w_mean_adj <- numeric(d)
for(i in 1:d){
  lm_data = data.frame("w"=w_post[burn_in:M,i], ssim_mat[burn_in:M,])
  lm_tmp <- lm(w ~ ., data = lm_data, na.action = na.omit)
  w_post_adj[[i]] <- sum(lm_tmp$coefficients * c(1,sobs), na.rm=TRUE) + lm_tmp$residuals
  w_mean_adj[i] <- mean(w_post_adj[[i]])
}
print("ABC-MCMC regression adjusted means:")
print(w_mean_adj)
print("ABC-MCMC adjusted mean relative error in %")
print((w_mean_adj-c(tau,gammavar))/c(tau,gammavar)*100)

############### Metropolis Hastings MCMC ##################
# Normal transition kernel with sd = rho 
set.seed(0)
start_time <- Sys.time()

# Hyperparameters (adjust these)
M <- 10000                  # Number of iterations 
burn_in <- 5000             # Burn in
trans_cv <- c(0.03, 0.03)   # transition kernel sd
w <- c(1, 0.00005)          # initial parameters

# Initialization (do not change)
d <- length(w)            # Number of parameters 
w_post <- matrix(NaN,M,d) # stores posterior samples  
acc_vec <- numeric(M)   
rho <- trans_cv*w         # sd of transition kernel

for(iter in 1:M){
  
  if(mod(iter,1000)==0){
    print(iter)
  } 
  
  # Sample from the transition kernel 
  wp <- rnorm(d, mean=w, sd=rho)
  
  # Uninformative (improper) flat prior  
  prior_ratio <- 1   
  
  # Exp growth model Likelihood ratio
  Ln_wp <- Ln_exp(wp, y)
  Ln_w  <- Ln_exp(w, y)
  likelihood_ratio <- exp(Ln_wp - Ln_w)
  if(is.nan(likelihood_ratio)){
    message(iter, ": Likelihood ratio: NaN")
  }
  if(is.infinite(likelihood_ratio)){
    message(iter, ": Likelihood ratio: INF")
  }   
  
  alpha <- prior_ratio*likelihood_ratio 
  alpha <- min(1, alpha, na.rm=TRUE)  # MH acceptence probability 
  
  U <- runif(1) 
  if(alpha > U){
    w_post[iter,] <- wp
    w <- wp 
    acc_vec[iter] <- 1
  } else {
    w_post[iter,] <- w
    acc_vec[iter] <- 0
  }  
}
w_MCMC <- w_post
MCMC_time <- Sys.time() - start_time
print("MCMC time")
print(MCMC_time)

# MCMC mixing plots #
par(mfrow=c(4,1))
plot(w_post[,1],type='l', ylab=expression(hat(tau)))
lines(w_post[,1])
acf(w_post[burn_in:M,1], lag=15, main=expression(hat(tau)))
plot(w_post[,2],type='l', ylab=expression(hat(eta)))
lines(w_post[,2])
acf(w_post[burn_in:M,2], lag=15, main=expression(hat(eta)))

# MCMC posterior means #
w_mean <- colMeans(w_post[burn_in:M,])
print("ABC-MCMC acceptance rate:")
print(mean(acc_vec[burn_in:M]))
print("ABC-MCMC means:")
print(w_mean)


################ POSTERIOR PLOTS (ALL) ##################
par(mfrow=c(2,3))
plot(density(w_post_ABC[burn_in:M,1], kernel="epanechnikov", 
             adjust=1.4), xlab=expression(paste(hat(m))), main="ABC")
plot(density(w_post_adj[[1]], kernel="epanechnikov", 
             adjust=1.0), xlab=expression(paste(hat(m))), main="Adjusted ABC")
plot(density(w_MCMC[2000:M,1], kernel="epanechnikov", 
             adjust=1.2), xlab=expression(paste(hat(m))), main="MCMC")
plot(density(w_post_ABC[burn_in:M,2], kernel="epanechnikov", adjust=1.4), xlab=expression(paste(hat(eta))), main="ABC")
plot(density(w_post_adj[[2]], kernel="epanechnikov", 
             adjust=1.0), xlab=expression(paste(hat(eta))), main="Adjusted ABC")
plot(density(w_MCMC[2000:M,2], kernel="epanechnikov", 
             adjust=1.2), xlab=expression(paste(hat(eta))), main="MCMC")



################# FUNCTIONS #################

M_exp <- function(nobs, tau, gammavar){
  tau1 <- tau
  tau2 <- tau
  gam  <- gammavar
  
  grr<-rgamma(nobs,tau1,rate=tau2)  # Gamma inverse growth rates
  u<-runif(nobs,0,1)
  v0<-(0.5^3)*pi/6
  vcd<-v0-(log(1-u))/(gam*grr)
  
  return(vcd)
}

Ln_exp <- function(w, vcd, negative=FALSE){
  if(any(w<=0))
    return(-Inf)
  tau <- w[1]
  tau1 <- tau
  tau2 <- tau
  gammavar <- w[2]
  v_cell<-(0.5^3)*pi/6
  if(negative)
    L <- -sum(log(tau1)+tau1*log(tau2)+log(gammavar)+(-(tau1+1))*log(tau2+gammavar*(vcd-v_cell)))
  else
    L <- sum(log(tau1)+tau1*log(tau2)+log(gammavar)+(-(tau1+1))*log(tau2+gammavar*(vcd-v_cell)))
  return(L)
}

summary_stats <- function(vcd, N=24){
  
  #defining categories for the tumour volume
  gupperbounds<-c(1.5,2.5,seq(7.5,72.5,5),85,seq(95,155,10))
  glowerbounds<-c(0.5,1.5,seq(2.5,68.5,5),72.5,seq(85,145,10))
  
  gsizeintervals<-length(gupperbounds)
  gmid<-(gupperbounds+glowerbounds)/2
  gmid<-(gupperbounds*glowerbounds)^(1/2)  
  nobs <- length(vcd)
  dcd<-(6*vcd/pi)^(1/3)
  gsize<-matrix(-1 ,nobs,1)
  for(s in 1:gsizeintervals){
    gsize[dcd >= glowerbounds[s] & dcd < gupperbounds[s]] <- s
  } 
  stat <- numeric(24)
  for(i in 1:24){
    stat[i] <- sum(gsize == i)
  }
  return(stat)
} 

