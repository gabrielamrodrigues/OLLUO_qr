rm(list=ls(all=TRUE))

#Models GAMLSS
source("OLLUO.q25.r")
source("OLLUO.q50.r")
source("OLLUO.q75.r")

library('gamlss')

# Function to generate regression model with 1 covariate

#tau=0.25
reg_OLLO.q25 <- function(n,beta10,beta11,beta20,beta21,beta30){
  x1 <- runif(n,0,1)
  mu <- (exp(beta10+beta11*x1))/(1+exp(beta10+beta11*x1))
  sigma <- exp(beta20+beta21*x1)
  nu <- exp(beta30)
  y <- rOLL.Omega.q25(n=n,mu=mu,sigma=sigma,nu=nu)
  return(list(y=y,x1=x1,mu=mu,sigma=sigma,nu=nu))
}

#tau=0.50
reg_OLLO.q50 <- function(n,beta10,beta11,beta20,beta21,beta30){
  x1 <- runif(n,0,1)
  mu <- (exp(beta10+beta11*x1))/(1+exp(beta10+beta11*x1))
  sigma <- exp(beta20+beta21*x1)
  nu <- exp(beta30)
  y <- rOLL.Omega.q50(n=n,mu=mu,sigma=sigma,nu=nu)
  return(list(y=y,x1=x1,mu=mu,sigma=sigma,nu=nu))
}

#tau=0.75
reg_OLLO.q75 <- function(n,beta10,beta11,beta20,beta21,beta30){
  x1 <- runif(n,0,1)
  mu <- (exp(beta10+beta11*x1))/(1+exp(beta10+beta11*x1))
  sigma <- exp(beta20+beta21*x1)
  nu <- exp(beta30)
  y <- rOLL.Omega.q75(n=n,mu=mu,sigma=sigma,nu=nu)
  return(list(y=y,x1=x1,mu=mu,sigma=sigma,nu=nu))
}


##Initial values 
beta10 <- 2.0; beta11 <- 1.5; beta20 <- 0.5; beta21 <- 0.9; beta30 <- 0.6
initial <- c(beta10,beta11,beta20,beta21,beta30)

#tau=0.25
n0 <- c(50,100,200,300,400,500,600,700,800,900,1000)
thetas.25 <- list()
rqs.25 <- list()
for(k in 1:length(n0)){
  iter <- 0
  r <- 1000
  theta <- matrix(0,r,5)
  res.rq = matrix(0,r,n0[k])
  i <- 1
  set.seed(1311)
  while(i<=r){
    dados1 <- reg_OLLO.q25(n0[k],beta10,beta11,beta20,beta21,beta30)
    x <- dados1$x
    y <- dados1$y
    mu0 <- (exp(beta10+beta11*x))/(1+exp(beta10+beta11*x))
    sigma0 <- exp(beta20+beta21*x)
    nu0 <- exp(beta30)
    ajuste1=try(gamlss(y~x,sigma.formula = ~x, family = OLL.Omega.q25,n.cyc=200,data=dados1,sigma.start = sigma0, mu.start = mu0,nu.start = nu0,trace=F))
    if((class(ajuste1)[1] != "try-error")==T){
      teste <- ajuste1$converged
      if(teste == TRUE){
        betas <- c(mu=ajuste1$mu.coefficients ,sigma=ajuste1$sigma.coefficients,nu=ajuste1$nu.coefficients)
        theta[i,] <-  betas
        res.rq[i,] = ajuste1$residuals
        i <- i+1
      }
      else{i <- i}
    }
    else{i <- i}
  }
  thetas.25[[k]] <- theta
  rqs.25[[k]] <- res.rq
}


#tau=0.50
n0 <- c(50,100,200,300,400,500,600,700,800,900,1000)
theta.50 <- list()
rqs.50 <- list()
for(k in 1:length(n0)){
  iter <- 0
  r <- 1000
  theta <- matrix(0,r,5)
  res.rq = matrix(0,r,n0[k])
  i <- 1
  set.seed(1311)
  while(i<=r){
    dados1 <- reg_OLLO.q50(n0[k],beta10,beta11,beta20,beta21,beta30)
    x <- dados1$x
    y <- dados1$y
    mu0 <- (exp(beta10+beta11*x))/(1+exp(beta10+beta11*x))
    sigma0 <- exp(beta20+beta21*x)
    nu0 <- exp(beta30)
    ajuste1=try(gamlss(y~x,sigma.formula = ~x, family = OLL.Omega.q50,n.cyc=200,data=dados1,sigma.start = sigma0, mu.start = mu0,nu.start = nu0,trace=T))
    if((class(ajuste1)[1] != "try-error")==T){
      teste <- ajuste1$converged
      if(teste == TRUE){
        betas <- c(mu=ajuste1$mu.coefficients ,sigma=ajuste1$sigma.coefficients,nu=ajuste1$nu.coefficients)
        theta[i,] <-  betas
        res.rq[i,] = ajuste1$residuals
        i <- i+1
      }
      else{i <- i}
    }
    else{i <- i}
  }
  thetas.50[[k]] <- theta
  rqs.50[[k]] <- res.rq
}

#tau=0.75
n0 <- c(50,100,200,300,400,500,600,700,800,900,1000)
thetas.75 <- list()
rqs.75 <- list()
for(k in 1:length(n0)){
  iter <- 0
  r <- 1000
  theta <- matrix(0,r,5)
  res.rq = matrix(0,r,n0[k])
  i <- 1
  set.seed(1311)
  while(i<=r){
    dados1 <- reg_OLLO.q75(n0[k],beta10,beta11,beta20,beta21,beta30)
    x <- dados1$x
    y <- dados1$y
    mu0 <- (exp(beta10+beta11*x))/(1+exp(beta10+beta11*x))
    sigma0 <- exp(beta20+beta21*x)
    nu0 <- exp(beta30)
    ajuste1=try(gamlss(y~x,sigma.formula = ~x, family = OLL.Omega.q75,n.cyc=200,data=dados1,sigma.start = sigma0, mu.start = mu0,nu.start = nu0,trace=F))
    if((class(ajuste1)[1] != "try-error")==T){
      teste <- ajuste1$converged
      if(teste == TRUE){
        betas <- c(mu=ajuste1$mu.coefficients ,sigma=ajuste1$sigma.coefficients,nu=ajuste1$nu.coefficients)
        theta[i,] <-  betas
        res.rq[i,] = ajuste1$residuals
        i <- i+1
      }
      else{i <- i}
    }
    else{i <- i}
  }
  thetas.75[[k]] <- theta
  rqs.75[[k]] <- res.rq
}


# Calculate results
results <- function(theta,n,initial){
  means <- apply(thetas.25[[1]], 2, mean)
  var <- apply(thetas.25[[1]], 2,var)
  MSEs <- c(var[1]+(means[1]-initial[1])^2,
            var[2]+(means[2]-initial[2])^2,
            var[3]+(means[3]-initial[3])^2,
            var[4]+(means[4]-initial[4])^2,
            var[5]+(means[5]-initial[5])^2)
  Biases <- c((means[1]-initial[1]),
              (means[2]-initial[2]),
              (means[3]-initial[3]),
              (means[4]-initial[4]),
              (means[5]-initial[5]))
  parameter <- rep(c('beta10', 'beta11','beta20', 'beta21','nu'))
  size <- rep(n,5)
  result <- cbind(initial,size,parameter, data.frame(means, Biases,MSEs))
  colnames(result) <- c('True value','n','Parameters', 'AEs','Biases', 'MSEs')
  return(result)
}

#For tau=0.25 and n=50
results(thetas.25[[1]],50,initial)
