load("C:/Users/Pulkit Jain/Desktop/ISEN__/626 Time Series/Project/Junho/reprojectround3/dat0705.RData")

#############presentation1#############
### log price centering
for(i in 1:8){
  logdat[,i] = logdat[,i]- mean(logdat[,i])
}

tslog = ts(logdat, frequency = 12, start = c(1987,5))
tsdf = ts(difflogdat, frequency = 12, start = c(1987,6))

#### CCF for difflogricce v.s. other variables(not include grains)
par(mfrow=c(2,3))
for(i in c(2,4:8)){
ccf(difflogdat$dfrice, difflogdat[,i], main = paste("ccf for dflogrice v.s. ",colnames(difflogdat)[i]))
}
par(mfrow=c(1,1))


par(mfrow=c(2,3))
for(i in c(1,4:8)){
  plot(logdat[,i], main = paste(colnames(logdat)[i]), type= 'l')
}
par(mfrow=c(1,1))


#### significant positive lags? corn(+1), wheat(+2), soybean(+2),uread(+6), dap(+1)
##############
###########################################
#############presentation2#################
###########################################
library(forecast)
library(tseries)
library(strucchange)

#######seasonality ~ periodogram######

n = length(logdat$logrice)
I = abs(fft(logdat$logrice))^2/n # the periodogram
P = (4/n)*I[1:((n+1)/2)] # the scaled periodogram
f = 0:((n-1)/2)/n # frequencies
plot(f, P, type="l", xlab="Frequency", ylab="Scaled Periodogram", main = "Periodogram for log(rice)")
### frequency is too small -> period is too long (~15years)


#automatic arima fitting result
fit.aic <- auto.arima(logdat$logrice, ic="aic")
fit.aic
fit.bic <- auto.arima(logdat$logrice, ic="bic")
fit.bic

####both suggests ARIMA(0,1,1), which is quiet obvious result###
#Series: dat$logrice 
#ARIMA(0,1,1)                    

#Coefficients:
#  ma1
#0.4143
#s.e.  0.0459

#sigma^2 estimated as 0.003228:  log likelihood=522.05
#AIC=-1040.1   AICc=-1040.07   BIC=-1032.33
####suggest ARIMA(0,1,1) which is quiet obvious result###
#########################################################

##residual distribution and MSE

qqnorm(fit.aic$residuals, main = "ARIMA(0,1,1) fitting residual QQplot")
qqline(fit.aic$residuals)
#mean(fit.aic$residuals)
#[1] 0.001265024


##### nonlinear fitting 

##### unit root test (augmented Dickey-Fuller test)
#adf.test(dat$logrice)
#adf.test(diff(dat$logrice))


##### change in structure after regression ####

#### significant positive lags? corn(+1), wheat(+2), soybean(+2),urea(+6), dap(+1)
dat_modif = ts.intersect(dfrice = tsdf[,1],dfcorn1 = lag(tsdf[,4],k=1),
                         dfwheat2 = lag(tsdf[,5],k=2), dfsoybean2 = lag(tsdf[,6],k=2),
                         dfurea6 = lag(tsdf[,7],k=6), dfdap1 = lag(tsdf[,8],k=1))

###first fitting linear regression to see what variables are significant
lm1 = lm(dfrice~., data = dat_modif)
#summary(lm1)
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.001111   0.003160   0.352   0.7253    
#dfcorn1      0.053261   0.056275   0.946   0.3446    
#dfwheat2    -0.016660   0.057577  -0.289   0.7725    
#dfsoybean2   0.119118   0.063090   1.888   0.0598 .  
#dfurea6     -0.081506   0.036181  -2.253   0.0249 *  
#dfdap1       0.288854   0.051110   5.652 3.32e-08 ***

temp <- breakpoints(dfrice~dfurea6+dfdap1, data = dat_modif, h = 36)
## change in structure h means minimum number of sample for each cluster
##(dont want (too) short term changes)
## this regression include intercept due to the inflation

#summary(temp)
round(coef(temp, breaks = 2),3) #assume two breakpoints and print where
  #                    Intercept) dfurea6 dfdap1
#1987(6) - 2004(7)        0.001   0.023  0.159
#2004(8) - 2008(4)        0.023  -0.361  0.410
#2008(5) - 2016(10)      -0.006   0.018  0.270

plot(tslog[,1], ylim = c(-0.6, 2.2), ylab="", main = "Change in Structure")
lines(tsdf[,1]+1.5, col=2)
abline(h=c(0,1.5), lty=c(2,2),col=c(1,2))
#lines(ts2[,2]+1)
#lines(fitted(temp, breaks=2)+1.5, col = 4)
lines(confint(temp, breaks = 2))

###prediction

###forcast without regression
dfrice=dat_modif[,1]
qt = as.numeric(time(dat_modif))
nbreak = which(qt==(2004+7/12)) ## 2005(8) first breakpoint

armamodel1 <- arima(dat_modif[1:nbreak,1],order = c(0,0,1), include.mean = T)

npred = length(dfrice)-nbreak-1
fcast1 <- forecast(armamodel1,h=npred)
plot(fcast1, ylim=c(-0.3,0.4),main = "Forecasts from MA(1)")

qq = as.numeric(dat_modif[,1])
abline(h=0,lty=2)
lines(qq,col=2)
lines(qq[1:nbreak],col=1)

###forecast after regression
armamodel2 <- arima(dat_modif[1:nbreak,1], xreg= dat_modif[1:nbreak,5:6],
                    order = c(0,0,1), include.mean = T)
fcast2 <- forecast(armamodel2, xreg=dat_modif[(nbreak+1):(nbreak+1+npred),5:6],h=npred)
plot(fcast2, ylim=c(-0.3,0.4),main = "Forecasts from MA(1) after regression")

qq = as.numeric(dat_modif[,1])
abline(h=0,lty=2)
lines(qq,col=2)
lines(qq[1:nbreak],col=1)

###forecast after fitting different regression

fcast3 = fcast1
nbreak2 = which(qt==(2008+4/12)) ## 2005(8) first breakpoint
cc = coef(temp, breaks = 2) #assume two breakpoints and print where
tt = rep(NA, npred)

tt[1:(nbreak2-nbreak)] = cc[2,1] + 
  cc[2,2]*dat_modif[(nbreak+1):(nbreak2),5] +
  cc[2,3]*dat_modif[(nbreak+1):(nbreak2),6]

tt[(nbreak2-nbreak+1):(npred)] = cc[3,1] + 
  cc[3,2]*dat_modif[(1+nbreak2):(npred+nbreak),2] +
  cc[3,3]*dat_modif[(1+nbreak2):(npred+nbreak),3]


tt = ts(tt, start=nbreak+1, end = nbreak+length(tt))
tq = fcast3$mean + tt
fcast3$mean = tq
tt2 = cbind(tt,tt)
fcast3$lower = fcast3$lower + tt2
fcast3$upper = fcast3$upper + tt2


plot(fcast3, ylim=c(-0.3,0.4),main = "Forecasts from MA(1) after piecewise regression")
qq = as.numeric(dat_modif[,1])
abline(h=0,lty=2)
lines(qq,col=2)
lines(qq[1:nbreak],col=1)

#### residuals

resid1 = qq[(nbreak+1):(nbreak+npred)] - fcast1$mean
resid2 = qq[(nbreak+1):(nbreak+npred+1)] - fcast2$mean
resid3 = qq[(nbreak+1):(nbreak+npred)] - fcast3$mean


fitarma1 = arima(resid1, order = c(0,0,1), include.mean = T)
fitarma2 = arima(resid2, order = c(0,0,1), include.mean = T)
fitarma3 = arima(resid3, order = c(0,0,1), include.mean = T)

qqnorm(fitarma1$residuals, main = "MA(1) only QQplot")
qqline(fitarma1$residuals)


qqnorm(fitarma2$residuals, main = "MA(1) after regerssion QQplot")
qqline(fitarma2$residuals)

qqnorm(fitarma3$residuals, main = "MA(1) after piecewise regerssion QQplot")
qqline(fitarma3$residuals)


###SSE
sum(resid1^2)
sum(resid2^2)
sum(resid3^2)



###########################################
#############presentation3#################
###########################################

#GARCH fitting.

library(fGarch)
resid3 = as.vector(resid3)
library("astsa")

acf2(resid3^2, main = "Series : residual^2 from M3")

fitgarch1 = garchFit(~garch(1,1), resid3, include.mean = F)
fitgarch2 = garchFit(~garch(2,1), resid3, include.mean = F)

summary(fitgarch1)
summary(fitgarch2)


u = fitgarch1@sigma.t

pred = as.vector(fcast3$mean)
lb = as.vector(fcast3$lower[,2])
ub = as.vector(fcast3$upper[,2])

plot(pred,type='l', ylim = c(-0.4,0.45), ylab = "", main = "GARCH(1,1) 95% Confidence Interval")
lines(pred+1.96*u, lty=2, col = 4, lwd =2)
lines(pred-1.96*u, lty=2, col = 4, lwd = 2)
lines(qq[(nbreak+1):352],col=2, type='l', lwd = 3)
legend("topright", legend = c("real value","95% C.I."), col = c(2,4), lwd = c(3,2), lty = c(1,2))

plot(pred,type='l', ylim = c(-0.4,0.45), ylab = "", main = "MA(1) 95% Confidence Interval")
lines(lb, lty=2, col = 4, lwd =2)
lines(ub, lty=2, col = 4, lwd = 2)
lines(qq[(nbreak+1):352],col=2, type='l', lwd = 3)
legend("topright", legend = c("real value","95% C.I."), col = c(2,4), lwd = c(3,2), lty = c(1,2))

plot(fitarma3$residuals, main="Residuals after Piecewise regression fit")
abline(h=0, lty=2, col="red")
acf(fitarma3$residuals, lag.max = 20, main = "ACF of Residuals after Piecewise regression fit")
acf(fitarma3$residuals, lag.max = 20, main = "PACF of Residuals after Piecewise regression fit", type="partial")
acf(fitarma3$residuals^2, type="partial")
acf(resid3^2, type="partial")

acf(fitarma3$residuals^2)
acf(resid3^2)

