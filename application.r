rm(list=ls(all=TRUE))

library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(readxl)
library(gamlss)
library(moments)
library(xtable)
library(knitr)

#Models GAMLSS

source("OLLUO.q05.r")
source("OLLUO.q10.r")
source("OLLUO.q25.r")
source("OLLUO.q50.r")
source("OLLUO.q75.r")
source("OLLUO.q90.r")
source("OLLUO.q95.r")

#Dataset

dta_gini_2 <- read_excel('data_gini.xlsx')

#Descriptive analysis

means <- mean(dta_gini_2$HDI)
desv <- sd(dta_gini_2$HDI)
min <- min(dta_gini_2$HDI)
max <- max(dta_gini_2$HDI)
skew <-skewness(dta_gini_2$HDI)
median <- median(dta_gini_2$HDI)
kurtosi <- kurtosis(dta_gini_2$HDI)
cv <- 100*desvios/medias

descriptive <- cbind(Mean=means,s.d.=desv,Median=median,Min.=min,Max=max, Skewness=skew, VC=cv,Kurtosis=kurtosi)
kable(descriptive, digits = 4)

x11(width = 5, height = 4.5)
ggplot(dta_gini_2, aes(x=HDI)) +
  geom_histogram(aes(y=..density..),color="black", position="identity", alpha=0.5)+xlab('x')+ylab('Density')+theme_bw()+theme(axis.text=element_text(size=11),strip.text.x = element_text(size = 11), axis.title=element_text(size=13))



## Model M0
fit.90_m0  <- gamlss(HDI~1, family = OLL.Omega.q90,n.cyc=800,data=dta_gini_2)
fit.75_m0  <- gamlss(HDI~1, family = OLL.Omega.q75,n.cyc=800,data=dta_gini_2)
fit.50_m0  <- gamlss(HDI~1, family = OLL.Omega.q50,n.cyc=800,data=dta_gini_2)
fit.25_m0  <- gamlss(HDI~1, family = OLL.Omega.q25,n.cyc=800,data=dta_gini_2)
fit.10_m0  <- gamlss(HDI~1, family = OLL.Omega.q10,n.cyc=800,data=dta_gini_2)
fit.005_m0  <- gamlss(HDI~1, family = OLL.Omega.q005,n.cyc=800,data=dta_gini_2)
fit.95_m0  <- gamlss(HDI~1, family = OLL.Omega.q95,n.cyc=800,data=dta_gini_2)

## Model M1
fit.90_m1  <- gamlss(HDI~GINI, family = OLL.Omega.q90,n.cyc=800,data=dta_gini_2)
fit.75_m1  <- gamlss(HDI~GINI, family = OLL.Omega.q75,n.cyc=800,data=dta_gini_2)
fit.50_m1  <- gamlss(HDI~GINI, family = OLL.Omega.q50,n.cyc=800,data=dta_gini_2)
fit.25_m1  <- gamlss(HDI~GINI, family = OLL.Omega.q25,n.cyc=800,data=dta_gini_2)
fit.10_m1  <- gamlss(HDI~GINI, family = OLL.Omega.q10,n.cyc=800,data=dta_gini_2)
fit.005_m1  <- gamlss(HDI~GINI, family = OLL.Omega.q005,n.cyc=800,data=dta_gini_2)
fit.95_m1  <- gamlss(HDI~GINI, family = OLL.Omega.q95,n.cyc=800,data=dta_gini_2)


## Model M2
fit.95_m2  <- gamlss(HDI~1,sigma.formula = ~GINI, family = OLL.Omega.q95,n.cyc=800,data=dta_gini_2)
fit.90_m2  <- gamlss(HDI~1,sigma.formula = ~GINI, family = OLL.Omega.q90,n.cyc=800,data=dta_gini_2)
fit.75_m2  <- gamlss(HDI~1,sigma.formula = ~GINI, family = OLL.Omega.q75,n.cyc=800,data=dta_gini_2)
fit.50_m2  <- gamlss(HDI~1,sigma.formula = ~GINI, family = OLL.Omega.q50,n.cyc=800,data=dta_gini_2)
fit.25_m2  <- gamlss(HDI~1,sigma.formula = ~GINI, family = OLL.Omega.q25,n.cyc=800,data=dta_gini_2)
fit.10_m2  <- gamlss(HDI~1,sigma.formula = ~GINI, family = OLL.Omega.q10,n.cyc=800,data=dta_gini_2)
fit.005_m2  <- gamlss(HDI~1,sigma.formula = ~GINI, family = OLL.Omega.q005,n.cyc=800,data=dta_gini_2)


## Model M3
fit.95_m3  <- gamlss(HDI~GINI,sigma.formula = ~GINI, family = OLL.Omega.q95,n.cyc=800,data=dta_gini_2)
fit.90_m3  <- gamlss(HDI~GINI,sigma.formula = ~GINI, family = OLL.Omega.q90,n.cyc=800,data=dta_gini_2)
fit.75_m3  <- gamlss(HDI~GINI,sigma.formula = ~GINI, family = OLL.Omega.q75,n.cyc=800,data=dta_gini_2)
fit.50_m3  <- gamlss(HDI~GINI,sigma.formula = ~GINI, family = OLL.Omega.q50,n.cyc=800,data=dta_gini_2)
fit.25_m3  <- gamlss(HDI~GINI,sigma.formula = ~GINI, family = OLL.Omega.q25,n.cyc=800,data=dta_gini_2)
fit.10_m3  <- gamlss(HDI~GINI,sigma.formula = ~GINI, family = OLL.Omega.q10,n.cyc=800,data=dta_gini_2)
fit.005_m3  <- gamlss(HDI~GINI,sigma.formula = ~GINI, family = OLL.Omega.q005,n.cyc=800,data=dta_gini_2)


#AIC
c(AIC(fit.005_m0),AIC(fit.10_m0),AIC(fit.25_m0),AIC(fit.50_m0),AIC(fit.75_m0),AIC(fit.90_m0),AIC(fit.95_m0))
c(AIC(fit.005_m1),AIC(fit.10_m1),AIC(fit.25_m1),AIC(fit.50_m1),AIC(fit.75_m1),AIC(fit.90_m1),AIC(fit.95_m1))
c(AIC(fit.005_m2),AIC(fit.10_m2),AIC(fit.25_m2),AIC(fit.50_m2),AIC(fit.75_m2),AIC(fit.90_m2),AIC(fit.95_m2))
c(AIC(fit.005_m3),AIC(fit.10_m3),AIC(fit.25_m3),AIC(fit.50_m3),AIC(fit.75_m3),AIC(fit.90_m3),AIC(fit.95_m3))


#Model M3 summary
summary(fit.005_m3, type='qr')
summary(fit.10_m3, type='qr')
summary(fit.25_m3, type='qr')
summary(fit.50_m3, type='qr')
summary(fit.75_m3, type='qr')
summary(fit.90_m3, type='qr')
summary(fit.95_m3, type='qr')


reg <- function(x, beta10,beta11){
  1/(1+exp(-(beta10+beta11*x)))
}

ggplot(data = dta_gini_2, aes(GINI, HDI)) +   
  geom_point() +xlab("GINI") +  ylab("HDI")+
  stat_function(fun = reg,args = list(0.75509483,0.04013317),aes(colour='0.05'),size = 1) +
  stat_function(fun = reg,args = list(0.7370969, 0.2408806),aes(colour='0.10'),size = 1) +
  stat_function(fun = reg,args = list(0.6906058, 0.5805297),aes(colour='0.25'),size = 1) +
  stat_function(fun = reg,args = list(0.6241878,0.9578462),aes(colour='0.50'),size = 1) +
  stat_function(fun = reg,args = list(0.5614646,1.3132873),aes(colour='0.75'),size = 1) +
  stat_function(fun = reg,args = list(0.5250049, 1.5905412),aes(colour='0.90'),size = 1) +
  stat_function(fun = reg,args = list(0.5154997 , 1.7317459),aes(colour='0.95'),size = 1)+
  labs(colour="Quantil", values=1)  +
  theme_bw()+
  theme(axis.title = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 11,),
        legend.title = element_text(size = 11),
        legend.background=element_blank(),legend.margin = margin(0,0,0,0))

#Residuals M3
qqnorm(fit.005_m3$residuals,pch=20)
abline(0,1,col='red',lwd=2)
qqnorm(fit.10_m3$residuals,pch=20)
abline(0,1,col='red',lwd=2)
qqnorm(fit.90_m3$residuals,pch=20)
abline(0,1,col='red',lwd=2)
qqnorm(fit.50_m3$residuals,pch=20)
abline(0,1,col='red',lwd=2)
qqnorm(fit.25_m3$residuals,pch=20)
abline(0,1,col='red',lwd=2)
qqnorm(fit.75_m3$residuals,pch=20)
abline(0,1,col='red',lwd=2)
qqnorm(fit.95_m3$residuals,pch=20)
abline(0,1,col='red',lwd=2)

## Index x residuals  
res005 <- fit.005_m3$residuals
res010 <- fit.10_m3$residuals
res025 <- fit.25_m3$residuals
res050 <- fit.50_m3$residuals
res075 <- fit.75_m3$residuals
res090 <- fit.90_m3$residuals
res095 <- fit.95_m3$residuals

#005
index <- seq(1:length(res005))
d1r <- data.frame(index,res005)
ggplot(d1r, aes(x=index, y=res005)) + ylim(-4.4,4.4)+
  geom_point()+ylab('Quantile residual')+
  xlab('Index')+
  geom_hline(yintercept = 3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = -3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = 0,color='darkblue',size=1)+
  theme(legend.position="none",panel.background = element_rect(fill = "white", colour = "grey50"),axis.title = element_text( size = (14)),axis.text = element_text(size=13))


#010
index <- seq(1:length(res010))
d1r <- data.frame(index,res010)
ggplot(d1r, aes(x=index, y=res010)) + ylim(-4.4,4.4)+
  geom_point()+ylab('Quantile residual')+
  xlab('Index')+
  geom_hline(yintercept = 3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = -3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = 0,color='darkblue',size=1)+
  theme(legend.position="none",panel.background = element_rect(fill = "white", colour = "grey50"),axis.title = element_text( size = (14)),axis.text = element_text(size=13))

#025
index <- seq(1:length(res025))
d1r <- data.frame(index,res025)
ggplot(d1r, aes(x=index, y=res025)) + ylim(-4.4,4.4)+
  geom_point()+ylab('Quantile residual')+
  xlab('Index')+
  geom_hline(yintercept = 3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = -3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = 0,color='darkblue',size=1)+
  theme(legend.position="none",panel.background = element_rect(fill = "white", colour = "grey50"),axis.title = element_text( size = (14)),axis.text = element_text(size=13))


#075
index <- seq(1:length(res075))
d1r <- data.frame(index,res075)
ggplot(d1r, aes(x=index, y=res075)) + ylim(-4.4,4.4)+
  geom_point()+ylab('Quantile residual')+
  xlab('Index')+
  geom_hline(yintercept = 3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = -3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = 0,color='darkblue',size=1)+
  theme(legend.position="none",panel.background = element_rect(fill = "white", colour = "grey50"),axis.title = element_text( size = (14)),axis.text = element_text(size=13))


#090
index <- seq(1:length(res090))
d1r <- data.frame(index,res090)
ggplot(d1r, aes(x=index, y=res090)) + ylim(-4.4,4.4)+
  geom_point()+ylab('Quantile residual')+
  xlab('Index')+
  geom_hline(yintercept = 3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = -3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = 0,color='darkblue',size=1)+
  theme(legend.position="none",panel.background = element_rect(fill = "white", colour = "grey50"),axis.title = element_text( size = (14)),axis.text = element_text(size=13))


#095
index <- seq(1:length(res095))
d1r <- data.frame(index,res095)
ggplot(d1r, aes(x=index, y=res095)) + ylim(-4.4,4.4)+
  geom_point()+ylab('Quantile residual')+
  xlab('Index')+
  geom_hline(yintercept = 3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = -3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = 0,color='darkblue',size=1)+
  theme(legend.position="none",panel.background = element_rect(fill = "white", colour = "grey50"),axis.title = element_text( size = (14)),axis.text = element_text(size=13))
