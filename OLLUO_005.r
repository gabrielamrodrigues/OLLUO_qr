require(numDeriv)


OLL.Omega.q005 <- function (mu.link = "logit", sigma.link="log", nu.link = "log")
{
  mstats <- checklink(   "mu.link", "OLL.Omega.q005", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "OLL.Omega.q005", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink(   "nu.link", "OLL.Omega.q005", substitute(nu.link),    
                         c("inverse", "log", "identity", "own"))
  structure(
    list(family = c("OLL.Omega.q005", "OLL Omega"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
         nopar = 3, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         nu.link = as.character(substitute(nu.link)), 
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         nu.dr = vstats$mu.eta,
         
         dldm = function(y,mu,sigma,nu){ #----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiOLL.Omega.q005(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           dldm
         },
         d2ldm2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiOLL.Omega.q005(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           d2ldm2 <- -dldm * dldm
           d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
           d2ldm2 
         },     
         dldd = function(y,mu,sigma,nu){#----------------------------------------------------- ok  
           lpdf<-function(t,mu,x,nu){log(dauxiOLL.Omega.q005(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,method='simple')
           dldd
         } ,
         d2ldd2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,mu,x,nu){log(dauxiOLL.Omega.q005(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,method='simple')
           d2ldd2 <- -dldd*dldd
           d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
           d2ldd2
         },   
         dldv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,mu,sigma,x){log(dauxiOLL.Omega.q005(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,method='simple')
           dldv
         },
         d2ldv2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok 
           lpdf<-function(t,mu,sigma,x){log(dauxiOLL.Omega.q005(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,method='simple')
           d2ldv2<- -dldv * dldv
           d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                   
           d2ldv2
         },
         d2ldmdd = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiOLL.Omega.q005(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           lpdf<-function(t,mu,x,nu){log(dauxiOLL.Omega.q005(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu)
           d2ldmdd = -(dldm * dldd)
           d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
           d2ldmdd                 
         },
         d2ldmdv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiOLL.Omega.q005(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           lpdf<-function(t,mu,sigma,x){log(dauxiOLL.Omega.q005(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu)
           d2ldmdv = -(dldm * dldv)
           d2ldmdv				
         },
         d2ldddv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,mu,x,nu){log(dauxiOLL.Omega.q005(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,method='simple')
           lpdf<-function(t,mu,sigma,x){log(dauxiOLL.Omega.q005(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu)
           d2ldddv = -(dldd * dldv)
           d2ldddv	
         },
         #----------------------------------------------------- ok
         G.dev.incr  = function(y,mu,sigma,nu,...) 
         { 
           -2*dOLL.Omega.q005(y,mu,sigma,nu,log=TRUE)
         } ,                     
         rqres = expression(   
           rqres(pfun="pOLL.Omega.q005", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)) ,
         #mu.initial = expression(mu <- (y+mean(y))/2), 
         #mu.initial = expression(mu <- y+mean(y)),
         #mu.initial = expression(mu <- rep(2, length(y))),
         mu.initial = expression(mu <- rep(median(y), length(y))), 
         sigma.initial = expression(sigma <- rep(1, length(y))), #OK
         nu.initial = expression(nu <- rep(1,length(y))), #)k
         mu.valid = function(mu) all(mu > 0 & mu < 1),
         sigma.valid = function(sigma) all(sigma > 0),
         nu.valid = function(nu) all(nu > 0),
         y.valid = function(y)  all(y > 0 & y < 1) 
    ),
    class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dOLL.Omega.q005 <- function(x, mu = 0.2, sigma = 2, nu = 0.5, log = FALSE){
  if (any(mu < 0)|any(mu> 1))  stop(paste("mu must be between 0 and 1", "\n", ""))    
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  
  q1 <- 0.05
  w <- (q1^(1/nu)/((1-q1)^(1/nu)+ q1^(1/nu)))
  w1 <- ((1-w)^(-1/sigma)-1)/((1-w)^(-1/sigma)+1)
  mu1 <- log(w1)/log(mu)
  
  G <- 1-((1+x^mu1)/(1-x^mu1))^(-sigma)
  g <- ((2*sigma*mu1*x^(mu1-1))/(1-x^(2*mu1)))*(((1+x^mu1)/(1-x^mu1))^(-sigma))
  f <- (nu*g*(G*(1-G))^(nu-1))/(((G^nu)+(1-G)^nu)^2)  
  
  if(log==FALSE) fy<-f else fy<-log(f)
  fy
}    
#-----------------------------------------------------------------  
pOLL.Omega.q005 <- function(q, mu = 0.2, sigma = 2, nu = 0.5, lower.tail = TRUE, log.p = FALSE){  
  if (any(mu < 0)|any(mu> 1))  stop(paste("mu must be between 0 and 1", "\n", ""))    
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  
  q1 <- 0.05
  w <- (q1^(1/nu)/((1-q1)^(1/nu)+ q1^(1/nu)))
  w1 <- ((1-w)^(-1/sigma)-1)/((1-w)^(-1/sigma)+1)
  mu1 <- log(w1)/log(mu)
  
  G <- 1-((1+q^mu1)/(1-q^mu1))^(-sigma)
  FF1 <- (G^nu)/((G^(nu)+(1-G)^nu))
  
  if(lower.tail==TRUE) cdf<-FF1 else cdf<- 1-FF1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#-----------------------------------------------------------------  
qOLL.Omega.q005 <-  function(p, mu=0.2, sigma=1, nu=1, lower.tail = TRUE, log.p = FALSE){   
  if (any(mu < 0)|any(mu> 1))  stop(paste("mu must be between 0 and 1", "\n", ""))    
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  
  q1 <- 0.05
  w <- (q1^(1/nu)/((1-q1)^(1/nu)+ q1^(1/nu)))
  w1 <- ((1-w)^(-1/sigma)-1)/((1-w)^(-1/sigma)+1)
  mu1 <- log(w1)/log(mu)
  
  l1 <- 1/nu
  u1 <- (p^(l1))/(((1-p)^l1)+p^l1)
  q <- (((1-u1)^(-1/sigma)-1)/((1-u1)^(-1/sigma)+1))^(1/mu1)
  q
}

#-----------------------------------------------------------------  
rOLL.Omega.q005 <- function(n, mu=0.2, sigma=1, nu=1){
  if (any(mu <= 0) | any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", "")) 
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  #  if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))
  uni<- runif(n = n,0,1)
  #r <-(log(nu)-(log(log((exp(nu)*(1-uni)+uni)))))^(1/sigma)/mu
  r <- qOLL.Omega.q005(uni, mu=mu, sigma=sigma, nu=nu)  
  r
}
#----------------------------------------------------------------- 
dauxiOLL.Omega.q005 <- function(t,mu,sigma,nu){ 
  q1 <- 0.05
  
  w <- (q1^(1/nu)/((1-q1)^(1/nu)+ q1^(1/nu)))
  w1 <- ((1-w)^(-1/sigma)-1)/((1-w)^(-1/sigma)+1)
  mu1 <- log(w1)/log(mu)
  
  G <- 1-((1+t^mu1)/(1-t^mu1))^(-sigma)
  g <- ((2*sigma*mu1*t^(mu1-1))/(1-t^(2*mu1)))*(((1+t^mu1)/(1-t^mu1))^(-sigma))
  fy1 <- (nu*g*(G*(1-G))^(nu-1))/(((G^nu)+(1-G)^nu)^2)  
  
  fy1}

