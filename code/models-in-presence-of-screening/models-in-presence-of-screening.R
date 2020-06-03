### OBS! All functions are included in the END! ###
### Libraries 
library("pracma") # mod, size functions

### DATA GENERATION (synthetic data) 
set.seed(1)    
n <- 1400000 # Number of individuals  
tau1 = 2.36 
tau2 = 4.16
eta=exp(-8.63)
beta1=-4.75
beta2=0.56
w_obs <- c(tau1, tau2, eta, beta1, beta2)
d <- length(w_obs)
sobs <- as.vector(simulate_model(n, tau1=tau1, tau2=tau2, 
                                 eta=eta, beta1=beta1, beta2=beta2))
n_stat <- length(sobs)

### Histograms of V|(A_d, M_d)
par(mfrow=c(2,2))
set.seed(0)
a <- simulate_model(n, tau1=tau1, tau2=tau2, eta=eta,
                    beta1=beta1, beta2=beta2)
rownames(a) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")
barplot((t(a)/colSums(a))[1,], ylim=c(0,0.20), col="blue", xlab="Tumour size bins", main="Sympt. det. age 40-42")
barplot((t(a)/colSums(a))[4,], ylim=c(0,0.20), col="blue", xlab="Tumour size bins", main="Screen det. age 42" )
barplot((t(a)/colSums(a))[2,], ylim=c(0,0.20), col="blue", xlab="Tumour size bins", main="Sympt. det. age 42-44" )
barplot((t(a)/colSums(a))[3,], ylim=c(0,0.20), col="blue", xlab="Tumour size bins", main="Sympt. det. age 44-48" )

### Estimation of standard deviation of summary statistics
set.seed(1)
start_time <- Sys.time()
n_mad <- 100
tmp <- matrix(1, n_mad, n_stat)
for(i in 1:n_mad){
  tmp[i,] <- as.vector(simulate_model(n,tau1=tau1,tau2=tau2, 
                                      eta=eta, beta1=beta1, beta2=beta2))
}
ssd <- apply(tmp, 2, sd)
print(Sys.time()-start_time)

### Finding an initial epsilon 
set.seed(1)
n_pilot <- 100
tmp <- numeric(n_pilot) 
for(i in 1:n_pilot){
  ssim <- as.vector(simulate_model(n, tau1=tau1*1.3, 
                                   tau2=tau2*1.3, eta=eta*1.3, 
                                   beta1=beta1*1.3, beta2=beta2*1.3))
  tmp[i] <- sqrt(sum(((ssim-sobs)/ssd)^2))/n_stat
}
quantile(tmp, 0.01)

### ABC-MCMC
# Euclidean distance metric weighted by sd estimates  
# Adaptive adjustment of tolerance  
# Flat improper prior
# Normal (symmetric) transition kernel
set.seed(0)
start_time <- Sys.time()

# Hyperparameters (adjust these)
M <- 20000                          # Number of iterations 
burn_in <- 8000                     # Burn in
trans_cv <- c(0.1,0.1,0.1,0.1,0.1)  # trans kernel variance
w <- 1.3*w_obs                      # initial parameters  
epsilon <- 0.40                     # Tolerance
quantile <- 80                      # For tolerance decay  
final_acc_rate <- 0.01              # Final tolerance 

# Initialization (do not change) 
w_post <- matrix(NaN,M,d)        # stores postrior values 
ssim <- numeric(n_stat)          # Stores temporary simulstat
ssim_mat <- matrix(NaN,M,n_stat) # Simul stats for all iter.
sdist_vec <- rep(NaN, M)         # Stores statistic distances   
n_accABC <- 0                    # Counter for ABC accptance 
acc_vec <- numeric(M)            # 1 if iter was acc 0 else  
acc_rate <- 1                       
rho <- abs(trans_cv*w_obs)       # sd of transition kernel

for(iter in 1:M){
  
  # print acceptance rates every 1000 iter
  if(mod(iter,1000)==0){
    print(iter)
    acc_rate <- mean(acc_vec[(iter-999):(iter-1)])
  } 
  # Plots of w over iterations (mixing)
  if(mod(iter,10)==0){
    par(mfrow=c(3,2))
    for(i in 1:d){
      plot(w_post[,i],type='l')
      lines(w_post[,i])
    }
  }
  # Sample from the transition kernel 
  wp <- rnorm(d, mean=w, sd=rho)
  
  # If tau2<0 or eta<0 reject
  if(wp[2]<=0 || wp[3]<=0){
    w_post[iter,] <- w
    acc_vec[iter] <- 0
    next    
  }
  # Simulate from the MODEL using the proposed parameter values 
  ssim <- as.vector(simulate_model(n, tau1=wp[1], tau2=wp[2], eta=wp[3], 
                                   beta1=wp[4], beta2=wp[5]))
  
  # Caclulate normalized euclidian distence between observed and simulated statistics 
  sdist <- sqrt(sum(((ssim-sobs)/ssd)^2))/n_stat
  
  # If the distance between the data statistics and the simulated statistics 
  # is bigger then epsilon, jump to next iteration
  if(sdist >= epsilon){
    w_post[iter,] <- w
    acc_vec[iter] <- 0
    next 
  }
  
  # Epsilon is set to the q:th quantile of accepted 
  # errors every 100th iteration   
  n_accABC <- n_accABC + 1
  ssim_mat[iter,] <- ssim 
  sdist_vec[n_accABC] <- sdist
  if(mod(n_accABC, 100) == 0 && acc_rate > final_acc_rate){
    sdist_sorted <- sort(sdist_vec[(n_accABC-99):n_accABC])
    epsilon <- sdist_sorted[quantile]
    cat(iter, "", epsilon, " \n") #Print the epsilon values
  }
  
  # # Uninformative (improper) flat prior 
  prior_ratio <- 1 
  
  if(is.nan(prior_ratio)){
    warning("Prior ratio: NAN")
  }
  if(is.infinite(prior_ratio)){
    warning("Prior ratio: INF")
  }
  
  alpha <- min(1, prior_ratio) # MH acceptence prob 
  
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
abc_time <- Sys.time() - start_time
print("ABC-MCMC time")
print(abc_time)

### MIXING PLOTS
# tau1, tau2, eta
par(mfrow=c(3,1))
plot(w_post[,1],type='l', ylab=expression(hat(tau)[1]))
lines(w_post[,1])
plot(w_post[,2],type='l', ylab=expression(hat(tau)[2]))
lines(w_post[,2])
plot(w_post[,3],type='l', ylab=expression(hat(eta)))
lines(w_post[,3])
# beta1, beta2
par(mfrow=c(2,1))
plot(w_post[,4],type='l', ylab=expression(hat(beta)[1]))
lines(w_post[,4])
plot(w_post[,5],type='l', ylab=expression(hat(beta)[2]))
lines(w_post[,5])

### ABC POSTERIOR MEANS
w_mean <- colMeans(w_post[burn_in:M,])
print("ABC-MCMC acceptance rate:")
print(mean(acc_vec[burn_in:M]))
print("ABC-MCMC means:")
print(w_mean)
print("ABC-MCMC mean relative error in %")
print((w_mean-w_obs)/w_obs*100)

### LINEAR REGRESSION ADJUSTMENT  
w_post_adj <- vector(mode="list", length=d)
w_mean_adj <- numeric(d)
for(i in 1:d){
  lm_data = data.frame("w"=w_post[burn_in:M,i], 
                       ssim_mat[burn_in:M,])
  lm_tmp <- lm(w ~ ., data = lm_data, na.action = na.omit)
  w_post_adj[[i]] <- sum(lm_tmp$coefficients * c(1,sobs), 
                         na.rm=TRUE) + lm_tmp$residuals
  w_mean_adj[i] <- mean(w_post_adj[[i]])
}
print("ABC-MCMC regression adjusted means:")
print(w_mean_adj)
print("ABC-MCMC regression adjusted means rel error in %")
print((w_mean_adj-w_obs)/w_obs*100)
w_ABC_adj <- w_post_adj
w_ABC <- w_post[burn_in:M,]

### POSTERIOR DESNITY PLOTS 
# tau1, tau2, eta
par(mfrow=c(3,2))
plot(density(w_ABC[,1], kernel="epanechnikov", adjust=1.5), 
     xlab=expression(paste(hat(tau)[1])), main="ABC")
plot(density(w_ABC_adj[[1]],kernel="epanechnikov",adjust=1), 
     xlab=expression(paste(hat(tau)[1])), main="Adjusted ABC")
plot(density(w_ABC[,2], kernel="epanechnikov", adjust=1.5), 
     xlab=expression(paste(hat(tau)[2])), main="ABC")
plot(density(w_ABC_adj[[2]],kernel="epanechnikov",adjust=1), 
     xlab=expression(paste(hat(tau)[2])), main="Adjusted ABC")
plot(density(w_ABC[,3], kernel="epanechnikov", adjust=1.5), 
     xlab=expression(paste(hat(eta))), main="ABC")
plot(density(w_ABC_adj[[3]],kernel="epanechnikov",adjust=1), 
     xlab=expression(paste(hat(eta))), main="Adjusted ABC")
# beta1, beta2
par(mfrow=c(2,2))
plot(density(w_ABC[,4],kernel="epanechnikov",adjust=1.5), 
     xlab=expression(paste(hat(beta)[1])), main="ABC")
plot(density(w_ABC_adj[[4]],kernel="epanechnikov",adjust=1), 
     xlab=expression(paste(hat(beta)[1])), main="Adjusted ABC")
plot(density(w_ABC[,5], kernel="epanechnikov", adjust=1.5), 
     xlab=expression(paste(hat(beta)[2])),main="ABC")
plot(density(w_ABC_adj[[5]],kernel="epanechnikov",adjust=1), 
     xlab=expression(paste(hat(beta)[2])),main="Adjusted ABC")


############## FUNCTIONS ###############
age_onset_cdf <-  function(t, u=0, A=-0.075, B=1.1*10^(-4), 
                           delta=0.5){
  return(1 - ((B-A)*exp(B*t)/(B*exp((B-A)*t)-A))^delta - u) 
}

inv_sampl <- function(u, A, B, delta){
  root <- uniroot(f=age_onset_cdf, interval=c(1e-5, 10000), u, A, B, delta)$root
  return(root)
}

ages_screen_model <- function(n){
  return(42)
}

age_onset_model <- function(n, A=-0.075, B=1.1*10^(-4), 
                            delta=0.5){
  dt <- 0.005
  t <- seq(0, 250, dt)
  dens <- A*B*delta*(B-A)^delta*exp(delta*B*t)*(1-exp((B-A)*t))/(
    B*exp((B-A)*t)-A)^(delta+1)
  onset_ages <- sample(t, n, replace=TRUE, p=dens)
  return(onset_ages)
}

growth_rate_model <- function(n, tau1=2.36, tau2=4.16){
  rates <- rgamma(n, tau1, rate=tau2)
  return(rates)
}

age_sd_model <- function(growth_rates, ages_onset, 
                         eta=exp(-8.63)){
  n <- length(growth_rates)
  u <- runif(n,0,1)
  v0 <- (0.5^3)*pi/6
  vdet <- v0 - (log(1-u))/(eta*growth_rates)
  tdet <- growth_rates*log(vdet/v0)
  age_det <- ages_onset + tdet 
  return(age_det)
}

screen_out_model <- function(growth_rates, ages_onset, 
                             ages_screen, beta1=-4.75, beta2=0.56){
  n <- length(growth_rates)
  S <- numeric(n)
  v0 <- (0.5^3)*pi/6
  a0 <- ages_onset
  vol <- v0*exp((ages_screen-a0)/growth_rates)
  diam <- (6*vol/pi)^(1/3)
  screen_sens <- 1/(1+exp(-(beta1+beta2*diam)))*(
    (sign(ages_screen-a0)+1)/2)
  u <- runif(n)
  S[a0<ages_screen & screen_sens>u] <- 1 
  return(S)
}

age_det_model <- function(screen_outs, ages_scr, ages_sd){
  age_scrdet <- ifelse(screen_outs, 42, Inf)
  ages_det <- pmin(age_scrdet, ages_sd)
  return(ages_det)
}

det_mode_model <- function(ages_det, ages_sd){
  det <- as.numeric(ages_sd == ages_det)
  return(det)
}

vol_det_model <- function(ages_onset,ages_det,growth_rates){
  v0 <- (0.5^3)*pi/6
  vol <- v0*exp((ages_det-ages_onset)/growth_rates)
  return(vol)
} 

cat_table <- function(vols_det, ages_det, det_modes, 
                      ages_screen){
  ind_coh <- ages_det >= 40 & ages_det <= 48
  vols <- vols_det[ind_coh]
  ages <- ages_det[ind_coh]
  detmodes <- det_modes[ind_coh]
  agesscr  <- 42
  vol_breaks <- c(0, 5, 80, 300, 700, 1500, 3000, 5000, 8000, 
                  12500,20000, 33000, 55000, 120000, Inf)
  n_volbreaks <- length(vol_breaks)  
  age_breaks <- c(40, 42, 44, 48)
  n_agebreaks <- length(age_breaks)
  n_screens <- length(agesscr)  
  n_agecat <- n_agebreaks-1+n_screens  
  n_volcat <- n_volbreaks-1
  table <- matrix(0, n_volcat, n_agecat)  
  for(j in 1:n_agecat){
    for(i in 1:n_volcat){
      if(j < n_agebreaks){
        tmp <- sum(detmodes == 1 &
                     vols >= vol_breaks[i] & 
                     vols <  vol_breaks[i+1] & 
                     ages >= age_breaks[j] & 
                     ages <  age_breaks[j+1]) 
      } else{
        tmp <- sum(vols >= vol_breaks[i] & 
                     vols <  vol_breaks[i+1] & 
                     ages == agesscr[j-n_agebreaks+1])
      }
      table[i,j] <- tmp
    }
  }
  return(table)
}

simulate_model <- function(nobs, tau1=2.36, tau2=4.16, 
                           eta=exp(-8.63), beta1=-4.75, beta2=0.56, 
                           A=-0.075, B=(1.1*10^(-4)), delta = 0.5){
  ages_screen <- ages_screen_model(nobs)
  growth_rates <- growth_rate_model(nobs, tau1, tau2) 
  ages_onset <- age_onset_model(nobs, A, B, delta)
  ages_sd <- age_sd_model(growth_rates, ages_onset, eta)
  screen_outs <- screen_out_model(growth_rates, ages_onset, 
                                  ages_screen, beta1, beta2)
  ages_det <- age_det_model(screen_outs,ages_screen,ages_sd)
  det_modes <- det_mode_model(ages_det, ages_sd) 
  vols_det <- vol_det_model(ages_onset, ages_det, 
                            growth_rates) 
  tab <- cat_table(vols_det, ages_det, det_modes, 
                   ages_screen)  
  return(tab)
}

