### OBS ALL FUNCTIONS IN THE END ###
### Libraries ###
library("pracma")  # mod, size functions
library("tseries") # ACF plots for MCMC 


########### DATA GENERATION (synthetic data) #########
set.seed(1)
start_time <- Sys.time()
n <- 10000  # Number of individuals  
m <- 1.07
v <- exp(log(1.31)/1.07*m)
gammavar <- exp(-9.602) # eta (old name) 
y <- M_log(n, m, v, gammavar)        
print(Sys.time()-start_time)
# Summary stats: counts in each category
sobs <- summary_stats(y)
print(sobs)


############## ML estimates #################
start_time <- Sys.time()              
opt <- optim(c(0.5, -4), fn=Ln_log, x=y, opt=TRUE)
w_ML <- c(opt$par[1], 10^(opt$par[2]))
ML_time <- Sys.time() - start_time
print(ML_time)
print(w_ML)


#################### ABC-MCMC ######################
# Euclidian distance weighted by V[s]  
# Adaptive adjustment of tolerance  
# Normal transition kernel with sd = rho 

### Estimation of V[s] ###
set.seed(0)
start_time <- Sys.time()
n_mad <- 200
tmp <- matrix(0,n_mad,24)
for(i in 1:n_mad){
  y_mad <- M_log(n, m, v, gammavar)
  tmp[i,] <- summary_stats(y_mad)
}
smad <- apply(tmp, 2, sd)
print(Sys.time()-start_time)
set.seed(0)
start_time <- Sys.time()

# Hyperparameters (adjust these)
M <- 10000                     # Number of iterations 
burn_in <- 5000                # Burn in
trans_cv <- c(0.1, 0, 0.1)     # transition kernel var
w <- 1.5*c(m, v, gammavar)     # initial parameters
epsilon <- 0.63                # tolerance
quantile <- 0.9         # Quantile for tol. tempering
final_acc_rate <- 0.01  # Stop crit. for tempering
tempering_freq <- 100   # Decay tol every 100 acc. iter.   

# Initialization (do not change)
d <- length(w)                    
n_stat <- length(sobs)           
w_post <- matrix(NaN,M,d)        # stores postrior values 
ssim <- numeric(n_stat)          
ssim_mat <- matrix(NaN,M,n_stat) # Stores simul statistics
sdist_vec <- rep(NaN, M)         # Stores statistic distances   
n_accABC <- 0                     
acc_vec <- numeric(M)            
acc_rate <- 1                    
rho <- abs(trans_cv*w)           # sd of transition kernel

for(iter in 1:M){
  
  if(mod(iter,100)==0){
    print(iter) 
  }
  # Plot of w over iterations (mixing)
  if(mod(iter,5)==0){
    par(mfrow=c(d,1))
    for(i in 1:d){
      plot(w_post[,i],type='l')
      lines(w_post[,i])
    }
  }
  if(mod(iter,500)==0){
    print(iter)
    acc_rate <- mean(acc_vec[(iter-499):(iter-1)])
  } 
  # Sample from the transition kernel
  wp <- rnorm(d, mean=w, sd=rho) 
  if(any(is.nan(wp)) | wp[2]<=0 
     | wp[3]<=0 | any(is.infinite(wp))){
    message(cat("WARNING: ", wp))
    w_post[iter,] <- w
    acc_vec[iter] <- 0
    next     
  }
  # Simulate from the MODEL using the proposed params  
  y_sim <- M_log(n, wp[1], wp[2], wp[3])
  
  # Calculate sufficient statistics from the samples
  ssim <- summary_stats(y_sim)
  
  # Normalized euclidean distance between observed and 
  # simulated statistics 
  sdist <- sqrt(sum(((ssim-sobs)/smad)^2))/n_stat 
  
  # If the dist between the data stats and the simul stats 
  # is bigger then epsilon, jump to next iter
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
  if(mod(n_accABC, tempering_freq)==0 
     && acc_rate>final_acc_rate){
    epsilon <- quantile(sdist_vec[(n_accABC-(
      tempering_freq-1)):n_accABC], 
      quantile, names=FALSE)
    cat(iter, "", epsilon, " \n") #Print the epsilon values
  }
  
  # Uninformative (improper) flat prior 
  prior_ratio <- 1
  
  alpha <- min(1, prior_ratio) # MH acceptance probability 
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

# ABC MIXING PLOTS #
par(mfrow=c(2,1))
plot(w_post[0:M,1],type='l', ylab=expression(hat(m)))
lines(w_post[0:M,1], ylab=expression(hat(m)))
plot(w_post[0:M,3],type='l', ylab=expression(hat(eta)))
lines(w_post[0:M,3], ylab=expression(hat(eta)))

# ABC POSTERIOR MEANS #
burn_in <- 4000
w_mean <- colMeans(w_post[burn_in:M,])
print("ABC-MCMC acceptance rate:")
print(acc_rate)
print("ABC-MCMC means:")
print(w_mean)
print("ABC-MCMC mean relative error in %")
print((w_mean-c(m,v*1.5,gammavar))/c(m,v,gammavar)*100)

### LINEAR REGRESSION ADJUSTMENT ###
w_post_adj <- vector(mode="list", length=d)
w_mean_adj <- numeric(d)
for(i in c(1,3)){
  lm_data = data.frame("w"=w_post[burn_in:M,i], ssim_mat[burn_in:M,])
  lm_tmp <- lm(w ~ ., data = lm_data, na.action = na.omit)
  w_post_adj[[i]] <- sum(lm_tmp$coefficients * c(1,sobs), na.rm=TRUE) + lm_tmp$residuals
  w_mean_adj[i] <- mean(w_post_adj[[i]])
}
print("ABC-MCMC regression adjusted means:")
print(w_mean_adj)
print("ABC-MCMC adjusted mean relative error in %")
print((w_mean_adj-c(m,v*1.5,gammavar))/c(m,v,gammavar)*100)


############## Metropolis Hastings MCMC #################
# Normal transition kernel with sd = rho 
set.seed(0)
start_time <- Sys.time()

# Hyperparameters (adjust these)
M <- 10000                  # Number of iterations 
burn_in <- 5000             # Burn in
trans_cv <- c(0.05, 0.05)   # transition kernel var
w <- 1.5*c(m, gammavar)     # initial parameters

# Initialization (do not change)
d <- length(w)                # Number of parameters 
w_post <- matrix(NaN,M,d)     # stores postrior values 
acc_vec <- numeric(M) 
rho <- trans_cv*c(m,gammavar) # sd of transition kernel

for(iter in 1:M){
  
  if(mod(iter,10)==0){
    print(iter)
  }
  # Mixing
  if(mod(iter,10)==0){
    par(mfrow=c(d,1))
    for(i in 1:d){
      plot(w_post[,i],type='l')
      lines(w_post[,i])
    }
  }  
  # Sample from the transition kernel 
  wp <- rnorm(d, mean=w, sd=rho)
  
  prior_ratio <- 1 # Flat prior
  
  # Exp growth model Likelihood ratio
  Ln_wp <- Ln_log(wp, y)
  Ln_w  <- Ln_log(w, y)
  likelihood_ratio <- exp(Ln_w - Ln_wp)
  if(is.nan(likelihood_ratio)){
    message(iter, ": Likelihood ratio: NaN")
  }
  if(is.infinite(likelihood_ratio)){
    message(iter, ": Likelihood ratio: INF")
  }   
  
  alpha <- prior_ratio*likelihood_ratio 
  alpha <- min(1, alpha, na.rm=TRUE)  # MH acceptance prob 
  
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

# MCMC MIXING PLOTS #
par(mfrow=c(2,2))
plot(w_post[,1],type='l', ylab=expression(hat(m)))
lines(w_post[,1])
acf(w_post[burn_in:M,1], lag=50, main=expression(hat(m)))
plot(w_post[,2],type='l', ylab=expression(hat(eta)))
lines(w_post[,2])
acf(w_post[burn_in:M,2], lag=50, main=expression(hat(eta)))
w_MCMC <- w_post

# MCMC POSTERIOR MEAN #
w_mean_MCMC <- colMeans(w_post[burn_in:M,])
print("ABC-MCMC acceptance rate:")
print(mean(acc_vec[burn_in:M]))
print("ABC-MCMC means:")
print(w_mean_MCMC)


########## ALL POSTERIOR PLOTS ############
par(mfrow=c(2,3))
par(cex.axis=1.3, cex.lab=1.5)
plot(density(w_post_ABC[burn_in:M,1], kernel="epanechnikov", 
             adjust=1.4), xlab=expression(paste(hat(m))), main="ABC")
plot(density(w_post_adj[[1]], kernel="epanechnikov", 
             adjust=1.0), xlab=expression(paste(hat(m))), main="Adjusted ABC")
plot(density(w_MCMC[burn_in:M,1], kernel="epanechnikov", 
             adjust=1.2), xlab=expression(paste(hat(m))), main="MCMC")
plot(density(w_post_ABC[burn_in:M,3], kernel="epanechnikov", 
             adjust=1.4), xlab=expression(paste(hat(eta))), main="ABC")
plot(density(w_post_adj[[3]], kernel="epanechnikov", 
             adjust=1.0), xlab=expression(paste(hat(eta))), main="Adjusted ABC")
plot(density(w_MCMC[burn_in:M,2], kernel="epanechnikov", 
             adjust=1.2), xlab=expression(paste(hat(eta))), main="MCMC")


###############  FUNCTIONS  #####################
M_log <- function(nobs, m, v, gammavar){
  v <- exp(log(1.31)/1.07*m)
  location <- log(m^2 / sqrt(v + m^2))
  shape <- sqrt(log(1 + (v / m^2)))
  gr <- rlnorm(n=nobs, location, shape)
  # Numerical inverse sampling 
  volcd <- apply(as.array(gr),
                 MARGIN=1, FUN=inv_sampl, 
                 gammavar=gammavar) 
  return(volcd)   # Volumes at symptomatic detection 
}

inv_sampl <- function(r, gammavar){
  u <- runif(1)
  vmax<-(128^3)*pi/6
  root <- uniroot(f=vsd_cdf, interval=c(0, vmax), r, gammavar, u)$root
  return(root)
}

summary_stats <- function(vcd, N=24){
  gupperbounds<-c(1.5,linspace(2.5,46,17), 48, 50, 55, 70, 90, 128)
  glowerbounds<-c(0.5,linspace(1.5,45,17), 46, 48, 50, 55, 70, 90)
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

# pdf of V_sd|R=r
vsd_pdf <- function(r, x, gammavar){
  vcell<-(0.5^3)*pi/6
  vmax<-(128^3)*pi/6
  a <- (vmax/vcell)^0.25 - 1
  b <- 0.25*r
  c <- vmax 
  M<-numeric(length(x))
  dM_dx<-numeric(length(x))
  if(all(x>0) && all(x<vmax)){
    M <- (log(a/(1-(x/c)^(0.25))) - log(a+1)  
          - 1/3*(x/c)^(0.75) - 1/2*(x/c)^(0.5) - (x/c)^(0.25) 
          - (11*a^3+27*a^2+18*a)/(6*(a+1)^3) + 11/6) 
    
    dM_dx <- 0.25/(c*(1-(x/c)^(0.25)))
  }
  return(gammavar*c/b*exp(-gammavar*c/b*M)*dM_dx*as.numeric(r>0)*as.numeric(gammavar>0))
}

# cdf of V_sd|R=r
vsd_cdf <- function(x, r, gammavar, u=0){
  n <- length(r)
  vcell<-(0.5^3)*pi/6
  vmax<-(128^3)*pi/6
  a <- (vmax/vcell)^0.25 - 1
  b <- 0.25*r
  c <- vmax 
  cdf<-numeric(n)-u
  if(x>0 && x<vmax){
    M <- (log(a/(1-(x/c)^(0.25))) - log(a+1)  
          - 1/3*(x/c)^(0.75) - 1/2*(x/c)^(0.5) - (x/c)^(0.25) 
          - (11*a^3+27*a^2+18*a)/(6*(a+1)^3) + 11/6) 
    cdf <- (1-exp(-gammavar*c/b*M))*as.numeric(r>0)*as.numeric(gammavar>0)-u 
  } else if(x>=vmax) {
    cdf <- 1-u
  } 
  return(cdf)
}

L_integr <- function(r, x, mv, gammavar){
  v <- exp(log(1.31)/1.07*mv)
  location <- log(mv^2 / sqrt(v + mv^2))
  shape <- sqrt(log(1 + (v / mv^2))) 
  return(vsd_pdf(r,x,gammavar) * dlnorm(r,location, shape))
}

L_loghelp <- function(x, m, gammavar){
  integrate(L_integr, lower=0, upper=Inf, x=x, mv=m, gammavar=gammavar, subdivisions=5000, rel.tol=1e-8, stop.on.error=FALSE)$value
}

L_log <- function(x, m, gammavar){
  apply(as.array(x), MARGIN=1, FUN=L_loghelp, m=m, gammavar=gammavar)
} 

# negative log likelihood  
Ln_log <- function(w, x, optim=FALSE){
  m <- w[1]
  if(optim){
    gammavar <- 10^(w[2])
  } else{
    gammavar <- w[2]
  }
  Ln_vec <- apply(as.array(x), MARGIN=1, FUN=L_log, m=m, gammavar=gammavar)
  return(-sum(log(Ln_vec)))  
}
