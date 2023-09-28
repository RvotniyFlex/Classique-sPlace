library("Rlab")    
library(ggplot2)
library(reshape2)
library(pracma)
library(actuar)
library(progress)
############################################################################
# Fuctions to generate process

#Counting process
generateCountingProcess <- function(T=100, lambda=5){

  xi <- rexp(ceiling(T*lambda*2), lambda)
  ts <- cumsum(xi)
  ts <- ts[ts<=T]
  return(ts)
}


# CPP
compPoissonTrend <- function(Nt, gamma=0, T, p=0.01,
                             ac_dist=rnorm(Nt(T), 0, 1),
                             disc_dist=rztpois(Nt(T), 15)) {
  # Nt - value vector of the counting process
  # gamma - the component responsible for the trend component of the process
  # T - time grid size
  # p (0 <= p <=1)- probability of the discrete jump
  # ac_dist - distribution class from which absolutely continuous rdist(Nt(T), *params) component is generated
  # disc_dist - distribution class from which the discrete component is generated

  jump_mask <- rbern(Nt(T), prob=p)
  Y_ac <- ac_dist
  Y_d <- disc_dist
  Yi <- Y_d*(jump_mask==1) + Y_ac*(jump_mask==0)
  Xt <- numeric(length=T+1)

  Xt[1] <- 0
  for(i in 1:T){
    Xt[i+1] <- gamma*i + sum(Yi[1:Nt(i)])
  }
  
  return(Xt)
}
###########################################################################
# Modeling (Example)



set.seed(123)
T <- 1000
lambda <- 1
mu <- 1
alpha <- 0.1

ts <- generateCountingProcess(T=T, lambda=10)
df <- data.frame(Nt=1:length(ts), t=ts)
df.melted <- melt(df, id.vars = "t")
Nt <- stepfun(ts, 0:length(ts))

# Visualize CPP
ggplot(df.melted, aes(x=t, y=value, color=variable)) + 
  geom_step() + 
  xlab("t") + 
  ylab("N(t)") + 
  scale_color_manual(values=c("red", "blue", "green", "purple"))



X <- compPoissonTrend(Nt, T=T)
df <- data.frame(t=0:T, X=X)
df.melted <- melt(df, id.vars = "t")

mean(tail(X, n=100))
ggplot(df.melted, aes(x=t, y=value, color=variable)) + 
  geom_step() + 
  #ggtitle("Sostavnoy protochniy process Poissona") + 
  xlab("t") + 
  ylab("X(t)") + 
  scale_color_manual(values=c("purple"))

############################################################################## 
# Вычислительная часть для параметра Lambda

make_simulations <- function(n_simulations= c(100, 500, 1000),
                             estimation_func,
                             n=100){
  # empty matrix to store estimated parameters
  
  data <- matrix(0, ncol= length(n_simulations), nrow=n)
  
  c <- 1
  for (j in n_simulations){
    
    out <-  numeric(length=n)
    for (i in 1:n){
      T <- j
      lambda <- 3
      mu <- 1
      alpha <- 0.1
      
      ts <- generateCountingProcess(T=T, lambda=lambda)
      df <- data.frame(Nt=1:length(ts), t=ts)
      df.melted <- melt(df, id.vars = "t")
      Nt <- stepfun(ts, 0:length(ts))
      
      X <- compPoissonTrend(Nt, T=T)
      Zk <- diff(X)
      
      
      out[i] <- estimation_func(Zk)
      
    }
    data[, c] <- out
    c <- c + 1
  }
  return(data)
}


makeExtendedSimulations <- function(lambdas=c(0.5, 1, 2.5, 5, 10),
                                     n_simulations= c(100, 500, 1000),
                                     estimation_func,
                                     n=100){
  errors <- matrix(0, ncol= length(n_simulations), nrow=n*length(lambdas))
  lambda_error <- matrix(0, ncol=length(n_simulations), nrow=length(lambdas))
  
  pb <- progress_bar$new(
    format = "[:bar] :percent :eta",
    total = length(n_simulations)
  )
 
  c <- 1
  for (j in n_simulations){
    i <- 1
    r <- 1
    window <- 1
    out <-  numeric(length=n*length(lambdas))
    for (lambda in lambdas){
      for (t in 1:n){
        T <- j
        lambda <- lambda
        mu <- 1
        alpha <- 0.1
        
        ts <- generateCountingProcess(T=T, lambda=lambda)
        df <- data.frame(Nt=1:length(ts), t=ts)
        df.melted <- melt(df, id.vars = "t")
        Nt <- stepfun(ts, 0:length(ts))
        
        X <- compPoissonTrend(Nt, T=T)
        Zk <- diff(X)
        out[i] <- abs(estimation_func(Zk) - lambda)
        i <- i + 1
      }
      
      t <- r*n 
     lambda_error[r,c] <- mean(out[window:t]/lambda*100)
     window <- t + 1
     r <- r + 1
    }
    errors[, c] <- out
    c <- c + 1
    pb$tick()
  }
  return(list("errod" = errors, "lamda_error" = lambda_error))
}


Zk <- diff(X)
#delta <- 3
#Zk <- Zk[seq(1, length(Zk), by=delta)]
#Kernel-based approximation

emp_char <- function(u) log(mean(exp(1i * u * Zk)))*(1 - u**2)**3
e_comp <- function(u) exp(-1i*u)
under_int <- function(u) emp_char(u)*e_comp(u)

#new <- function(x) integrate(Re(under_int(u, x)), -1, 1)$value
#new(u=1)

#ans <- -1/(2*pi)*integrate(new, -0.0005, 0.0005)$value
hn = 0.0005
M = (0.1221*hn)**(-1)
epsn = 2*pi*hn
Hn = epsn*M
#1/(2*pi*delta)*quad2d(function(u, x) Re(under_int(u, x)), -1, 1, -epsn, -Hn) + 1/(2*pi*delta)*quad2d(function(u, x) Re(under_int(u, x)), -1, 1, epsn, Hn) 
  

ans <- integrate(function(x) { 
  sapply(x, function(x) {
    integrate(function(u) Re(under_int(u)), -1, 1)$value
  })}, -Hn, -epsn, subdivisions=500)$value +
  integrate(function(x) { 
    sapply(x, function(x) {
      integrate(function(u) Re(under_int(u)), -1, 1)$value
    })}, epsn, Hn, subdivisions=500)$value
-ans/(2*pi)



KernelBasedEst <- function(Zk){
  emp_char <- function(u) log(mean(exp(1i * u * Zk)))*(1 - u**2)**2
  e_comp <- function(u, x) exp(1i*u*x)
  under_int <- function(u, x) emp_char(u)*e_comp(u, x)
  
  #new <- function(x) integrate(Re(under_int(u, x)), -1, 1)$value
  #new(u=1)
  
  #ans <- -1/(2*pi)*integrate(new, -0.0005, 0.0005)$value
  hn = 0.0005
  M = (0.1221*hn)**(-1)
  epsn = 2*pi*hn
  Hn = epsn*M
  #1/(2*pi*delta)*quad2d(function(u, x) Re(under_int(u, x)), -1, 1, -epsn, -Hn) + 1/(2*pi*delta)*quad2d(function(u, x) Re(under_int(u, x)), -1, 1, epsn, Hn) 
  
  
  ans <- integrate(function(x) { 
    sapply(x, function(x) {
      integrate(function(u) Re(under_int(u,x)), -1, 1)$value
    })}, -Hn, -epsn, subdivisions=500)$value +
    integrate(function(x) { 
      sapply(x, function(x) {
        integrate(function(u) Re(under_int(u,x)), -1, 1)$value
      })}, epsn, Hn, subdivisions=500)$value
  return(-ans)
}

c = 0.0156379237
KernelBasedEstGamm <- function(Zk){
  hn = 0.0005
  
  M = (0.1221*hn)**(-1)
  epsn = 2*pi*hn
  Hn = epsn*M
  emp_char <- function(u) log(mean(exp(1i * u * Zk)))*(1 - u**2)**2
  e_comp <- function(u, x) exp(1i*u*x)
  under_int <- function(u, x) x*emp_char(u)*e_comp(u, x)
  ans <- -integrate(function(x) { 
    sapply(x, function(x) {
      integrate(function(u) Re(under_int(u,x)), -1, 1)$value
    })}, -hn, hn, subdivisions=500)$value 
  return(ans)
}
ans_gamma <- KernelBasedEstGamm(Zk)/(c)
ans_gamma



KernelBasedEst(Zk)




dat <- make_simulations(estimation_func=KernelBasedEst)


dat1 <- makeExtendedSimulations(estimation_func = KernelBasedEst, n=50)
dat1$errod

save(simmulations, file="data.Rda")

  dat1
simmulations <- data.frame(dat)
simmulations$names
colnames(simmulations) <- c('100', '500', '1000')
save(simmulations, file="data.Rda")

save(dat1, file="dat1.Rda")


vec <- c(rep(c(0.5, 1, 2.5, 5, 10), each = 50))
est_errors <- data.frame(dat1$errod)
colnames(est_errors) <- c('100', '500', '1000')

est_errors$lambda <- vec
boxplot(dat1$errod, main = "Boxplot of Estimated Parameters", xlab = "T", ylab = "Absolute Error")




boxplot(sim_100 ~ lambda,        # Base R boxplot with manual y-axis
        est_errors[1:250, c(1,2,3, 4)],
        ylim = c(- 0, 5))

boxplot(sim_500 ~ lambda,        # Base R boxplot with manual y-axis
        est_errors[1:250, c(1,2,3, 4)],
        ylim = c(- 0, 5))

boxplot(sim_1000 ~ lambda,        # Base R boxplot with manual y-axis
        est_errors[1:250, c(1,2,3, 4)],
        ylim = c(- 0, 5))

boxplot(est_errors[1:50, c(1,2,3)], main = "Lambda = 0.5", xlab = "T", ylab = "Absolute Error", ylim = c(- 0, 2))

boxplot(est_errors[51:100, c(1,2,3)], main = "Lambda = 1", xlab = "T", ylab = "Absolute Error", ylim = c(- 0, 2))

boxplot(est_errors[101:150, c(1,2,3)], main = "Lambda = 2.5", xlab = "T", ylab = "Absolute Error", ylim = c(- 0, 5))

boxplot(est_errors[151:200, c(1,2,3)], main = "Lambda = 5", xlab = "T", ylab = "Absolute Error", ylim = c(- 0, 5))

boxplot(est_errors[201:250, c(1,2,3)], main = "Lambda = 10", xlab = "T", ylab = "Absolute Error", ylim = c(- 0, 5))

errors <- data.frame(dat1$lamda_error)
colnames(errors) <- c('100', '500', '1000')
rownames(errors) <- c(0.5, 1, 2.5, 5, 10)
errors


# Построение boxplot
ggplot(est_errors, aes(x = lambda)) +
  geom_boxplot(aes(y = sim_100, fill = "X")) +
  geom_boxplot(aes(y = sim_500, fill = "Y")) +
  geom_boxplot(aes(y = sim_1000, fill = "Z")) +
  xlab("Category") + ylab("Values") +
  scale_fill_manual("", values = c("X" = "red", "Y" = "blue", "Z" = "green"),
                    labels = c("X", "Y", "Z")) +
  theme_classic()


ggplot(est_errors, aes(x = "", y = sim_100, fill = "X")) +
  geom_boxplot() +
  geom_boxplot(aes(x = "", y = sim_500, fill = "Y")) +
  geom_boxplot(aes(x = "", y = sim_1000, fill = "Z")) +
  facet_wrap(~ lambda, ncol = 1) +
  xlab("") + ylab("Values") +
  scale_fill_manual("", values = c("X" = "red", "Y" = "blue", "Z" = "green"),
                    labels = c("X", "Y", "Z")) +
  theme_classic()
