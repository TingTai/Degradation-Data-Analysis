library(nlme)
library(ggplot2)
library(dplyr)
library(invgamma)
library(stats4)
library(MASS)
library(reshape2)
data("Fatigue")

Fatigue$relLength <- Fatigue$relLength*0.9

cycles2 <- Fatigue$cycles^2

Fatigue <- cbind(Fatigue,cycles2) %>% as.data.frame()

cycle_incre = Fatigue$cycles + 1e-10

Fatigue <- Fatigue %>% mutate(cycle_incre)

fit_s <- nlme::lmList(relLength ~ cycles + I(cycles^2) | Path, data = Fatigue)
fit_as <- summary(fit_s)
fit_as$coefficients
theta_q <- fit_as$coefficients %>% as.data.frame() %>% 
  subset(select=c('Estimate.(Intercept)', 'Estimate.cycles','Estimate.I(cycles^2)'))
colnames(theta_q) = c('Intercept', 'cycles', 'cycles2')


res <- residuals(fit_s)

sd(res) #0.01336842

#method1 dataset

f = function(t, beta2, beta1, beta0){
  return(beta2*t^2+beta1*t+beta0-1.6)
}

f2 = function(dat){
  root = uniroot(f, c(0,5), beta2=dat[3], beta1=dat[2], 
                 beta0=dat[1] )
  return(root$root)
}

lifetime_est = apply(theta_q, 1, f2)

lifetime_1 = matrix(lifetime_est[1:12], ncol = 1, nrow = 12)
lifetime_2 = subset(Fatigue, Path %in% c(13:21) & cycles %in% c(0.12))

T = matrix(c(lifetime_1, lifetime_2$cycles), ncol=1, nrow=21)
Path = seq(1,21,1)
relLength_1 = rep(1.6,12)
relLength_2 = c(1.52, 1.45, 1.49, 1.40, 1.38, 1.35, 1.31, 1.29, 1.27)
matrix_1 = matrix(relLength_1, nrow = 12, ncol=1)
matrix_2 = matrix(relLength_2, nrow =9, ncol = 1)
relLength = rbind(matrix_1, matrix_2)
T = cbind(T, Path) %>% as.data.frame()

obj_ = function(x){
  if(x==0.12) {return(0)}
  else {return(1)}
}

status = sapply(T, obj_)

Fati = cbind(T, status) %>% as.data.frame()

colnames(Fati) <- c('T', 'status')



#method 1
#MLE
log_all <- function(dat, pars){
  
  t = dat$T
  D = 1.6
  beta0 = exp(pars[1])
  beta1 = exp(pars[2])
  lamb = 38.880418
  k = 2.481892
  
  right = ((-beta1 * t - beta0 + D) / (lamb * t^2))^k
  log_right = log(1-exp(-right))*(1-dat$status)
  left_1 = -k* (((-beta1 * t - beta0 +D)/lamb * t^2)^(k-1))
  left_2 = (beta1*t+2*beta0-2*D)/(lamb*t^3)
  left_3 = exp(-(((-beta1 * t - beta0 + D) / (lamb * t^2))^k))
  log_left = (log(left_1 * left_2 * left_3))*(dat$status)
  log_wei = sum(log_right+log_left)
  
  if(is.na(log_wei)==TRUE) {return(1e+10)}
  else  {return(log_wei)}
}

est_mle_1=optim(c(0.1, 0.1), fn=log_all, dat = Fati, hessian =T, method = 'BFGS',control=list(fnscale=-1))
#a = est_mle$hessian

a_1 = est_mle_1$par %>% exp()

c_1 = c(38.880418, 2.481892)

-----------------------------------------------------------------------------------------------------------------
#method 2
#dataset

lifetime = matrix(lifetime_est, ncol=1, nrow = 21) %>% as.data.frame()
colnames(lifetime)<- c('cycles')
relLength = matrix(1.6, ncol=1, nrow = 21)
lifetime = cbind(Path, lifetime, relLength)
Fatigue_out = rbind(Fatigue, lifetime[c(13:21),])

#mle

log_out <- function(dat, pars){
  
  t = dat$cycles
  D = 1.6
  beta0 = exp(pars[1])
  beta1 = exp(pars[2])
  lamb = 38.880418#exp(pars[3])
  k = 2.481892 #exp(pars[4]) 
  
  left_1 = -k* (((-beta1 * t - beta0 +D)/lamb * t^2)^(k-1))
  left_2 = (beta1*t+2*beta0-2*D)/(lamb*t^3)
  left_3 = exp(-(((-beta1 * t - beta0 + D) / (lamb * t^2))^k))
  log_left = log(left_1 * left_2 * left_3)
  log_wei = sum(log_left)
  
  return(log_wei)
}

est_mle_2=optim(c(0.1, 0.1), fn=log_out, 
              dat = lifetime, hessian =F, method='BFGS',
              control=list(fnscale=-1))

b_1 = est_mle_2$par %>% exp() #1.457673 1.182454

log_out(lifetime, b_1)

c(a_1, b_1) %>% round(4)

-----------------------------------------------------------------------------------------------------------------
#method 3

like <- function(dat, pars){
  
  beta0 = exp(pars[1])
  beta1 = exp(pars[2])
  lamb = exp(pars[3])
  k = exp(pars[4])
  sigma = exp(pars[5])
  
  f = function(beta2){
    
    norm_f = dnorm(y, mean = beta0 + beta1 * t + beta2 * t^2, sd = sigma)#sd 0.01336842
    weibull = dweibull(beta2, shape = k, scale = lamb)
    f_all = prod(norm_f)*weibull
    return(f_all)
  }
  
  Li = array(1:21)
  
  for(i in 1:21){
    dt = subset(Fatigue, Path == i)
    
    y = dt$relLength
    
    t = dt$cycles
    
    L = integrate(Vectorize(f),0,60)$value
    
    L = -log(L)
    
    Li[i] = L
    
  }
  Li = Li %>% unlist
  # print(Li)
  L_all = sum(Li)
  return(L_all)
}

thetas = function(dat){
  d = c(0.1, 0.1, 0.1,0.1,0.1)
  for(i in 1:10){
    est_mle=optim(d, fn=like, 
                  dat = dat, hessian = T )
    d = est_mle$par
  }
  return(d %>% exp)
}

#-logLikelihood
d = thetas(Fatigue) #0.9060, 1.9497, 2.2208, 35.0164, 0.01756
d = log(d)



#
logL=like(Fatigue, d)
aic_3 = 2*5+2*logL
logL
aic_3

#cdf plot/ mean curve

len <- matrix(rep(NA,21*13), nrow=21 )

for(i in 1:21){
  index <- which(Fatigue$Path==i)
  len[i,1:length(index)] <- Fatigue[index,3]
}

rownames(len) <- c(1:21)
colnames(len) <- seq(from=0, to=0.12, by=0.01)


#mean curve
plot(x=seq(from=0, to=0.12, by=0.01), len[1,], type="l", 
     xlim=c(0,0.12), ylim=c(0.8,1.8),
     xlab="Millions of Cycles", ylab="Crack Length (inches)", main="data(Fatigue)")

for(i in 2:21){
  lines(x=seq(from=0, to=0.12, by=0.01), y=len[i,])
}
abline(h=1.6)
abline(v=0.12)

a_1 #mle for method 1
b_1 #mle for method2
c_1 #mle for weibull lambda and k

lambda = 38.880418
k = 2.481892

wei_mean = lambda * gamma(1+1/k)

wei_mean_1 = 35.0164 * gamma(1+1/2.2208)


t = seq(0, 0.12, by=0.01)

y_mean_1 = a_1[1]+a_1[2]*t + wei_mean * t^2

y_mean_2 = b_1[1]+b_1[2]*t + wei_mean * t^2

y_mean_3 = d[1]+d[2]*t + wei_mean_1 * t^2


#mean curve
plot(x=seq(from=0, to=0.12, by=0.01), len[1,], type="l", 
     xlim=c(0,0.12), ylim=c(0.8,1.8),
     xlab="Millions of Cycles", ylab="Crack Length (inches)", main="mean curve for different methods")

for(i in 2:21){
  lines(x=seq(from=0, to=0.12, by=0.01), y=len[i,])
}
abline(h=1.6)
abline(v=0.12)

lines(x = t, y = y_mean_1, type="l", col=2, lwd=2)
lines(x = t, y = y_mean_2, type="l", col='blue', lwd=2)
lines(x = t, y = y_mean_3, type="l", col='darkgreen', lwd=2)
legend(x = 'topleft', legend = c('method 1', 'method 2', 'method 3'), col = c('red', 'blue', 'darkgreen'), lty = 1)


#cdf plot
life.true <- c(0.088, 0.1, 0.101, 0.103, 0.103, 0.106, 0.106, 0.109,
               0.113, 0.115, 0.118, 0.118, 0.129, 0.133, 0.138, 0.144, 0.146,
               0.151, 0.160, 0.167, 0.170)

F_x <- seq(from=1/21, to=1, by=1/21)  
t1 <- seq(0, 0.20, length=10000)
F_t_1 = exp(-((-a_1[2] * t1 - a_1[1] + 1.6) / (lambda * t1^2))^k)
F_t_2 = exp(-((-b_1[2] * t1 - b_1[1] + 1.6) / (lambda * t1^2))^k)
F_t_3 = exp(-((-d[2] * t1 - d[1]  + 1.6) / (d[4] * t1^2))^ d[3])


plot(x=c(0,life.true,0.18), y=c(0,F_x,1), type="s", xlim=c(0,0.18),
     main="kaplan-meier and estimated F(t)", xlab="t(millions cycles)", ylab="F(t)")
points(x=life.true, y=F_x, pch=16, cex=0.5)
lines(t1, F_t_1, type="l", col=2)
lines(t1, F_t_2,type="l", col='blue')
lines(t1, F_t_3,type="l", col='darkgreen')

legend(x = 'topleft', legend = c('method 1', 'method 2', 'method 3'), col = c('red', 'blue', 'darkgreen'), lty = 1)


#final comparison

#cdf plot
plot(x=c(0,life.true,0.18), y=c(0,F_x,1), type="s", xlim=c(0,0.18),
     main="kaplan-meier and estimated F(t)", xlab="t(millions cycles)", ylab="F(t)")
points(x=life.true, y=F_x, pch=16, cex=0.5)
lines(t1, F_t_1, type="l", col='red',lwd=1.5)
lines(t1, F_t_2,type="l", col='blue',lwd=1.5)
lines(t_seq, t_cdf_1, type="l", col='darkgreen', lwd=1.5) #Model_1
lines(t_seq, t_cdf_2, type="l", col='purple', lwd=1.5) #Model_2
lines(t_seq, t_cdf_3, type="l", col='orange', lwd=1.5) #Model_3
legend(x = 'topleft', legend = c('regression 1,D=0.16',
                                 'regression 2,D=0.10', 
                                 'stochastic 1,D=0.61', 
                                 'stochastic 2,D=0.54', 
                                 'stochastic 3,D=0.63'), 
       col = c('red', 'blue','darkgreen', 'purple', 'orange' ), lty = 1, cex = 0.8)


#mean curve
plot(x=seq(from=0, to=0.12, by=0.01), data[1,], xlim=c(0,0.12), ylim=c(0.8,1.8), 
     type="b", col="gray", cex=0.2, pch=16,    
     xlab="Millions of Cycles", ylab="Crack Length (inches)", main="Mean Curve for Different Models")
for(i in 2:21){
  lines(x=seq(from=0, to=0.12, by=0.01), y=data[i,], 
        type="b", col="gray", cex=0.2, pch=16)
}
abline(h=1.6)
abline(v=0.12)

t_seq <- seq(0, 0.12, by=0.01)
lines(x = t, y = y_mean_1, type="l", col='red', lwd=2)
lines(x = t, y = y_mean_2, type="l", col='blue', lwd=2)
lines(t_seq, 0.9+(par_opt_1[1]*t_seq+par_opt_1[2]*t_seq^2)*par_opt_1[3], 
      lwd=2, type="l", col='darkgreen') #Model_1
lines(t_seq, 0.9+(par_opt_2[1]*t_seq^par_opt_2[2])*par_opt_2[3],
      lwd=2, type="l", col='purple') #Model_2
lines(t_seq, 0.9+(-0.28)+(par_opt_3[1]*par_opt_3[2]^t_seq)*par_opt_3[3], 
      lwd=2, type="l", col='orange') #Model_3

legend(x = 'topleft', legend = c('regression 1', 'regression 2', 'stochastic 1', 'stochastic 2', 'stochastic 3'), 
       col = c('red', 'blue','darkgreen', 'purple', 'orange' ), lty = 1)


