#dat = read.csv("C:/Users/junhoyang/Desktop/stat626/project/data/dat(final).csv")
load("C:/Users/Pulkit Jain/Desktop/ISEN__/626 Time Series/Project/Data/dat.RData")
rm(qq)


library("tseries")
library("astsa")
library("stats")
library(ggplot2)
library(forecast)
library(strucchange)
# install.packages("fGarch")
library(fGarch)


######presentation1#############
### centering
dat$logrice = dat$logrice - mean(dat$logrice)
dat$logoil = dat$logoil - mean(dat$logoil)

difflag = function(x,lag=1){
  n = length(x)
  temp = x[(1+lag):n]-x[1:(n-lag)]
  return(temp)
}

### difference lag=1
ricelag1 = difflag(dat$logrice,lag=1)
oillag1 = difflag(dat$logoil,lag=1)

ts1 = ts(dat[,c("rice","oil_real","logrice","logoil")], 
         frequency = 12, start = c(1987,4))
ts2 = ts(data.frame(ricelag1,oillag1), 
         frequency = 12, start = c(1987,5))

### raw data and data transformation
par(mfrow=c(3,1))
plot(ts1[,1],type='l', xlab="Year",
     ylab="Price ($ / metric ton)", main = "rice price")

abline(v=2008, col="red", lty=2)
text(2007.5, 800, "2008", srt=90, col="red") 

plot(ts1[,3],type='l', ylab="", main = "log(rice price) (centered)")
abline(h=0, lty=2)
plot(ts2[,1],type='l', ylab="", main = "Difference of log(rice price)")
abline(h=0, lty=2, col="red")
par(mfrow=c(1,1))
graphics.off()

#### moving average plot
plot(ts1[,3], type='l', ylab="", 
     main = "log(rice price) with Moving Averages")
for(w in c(3, 10, 25)) lines(filter(ts1[,3],filter=rep(1/w, w)), col=w, lwd=2)
legend("topleft",legend=c("MA3","MA10","MA25"), col=c(3,10,25), lwd=2)

#### ACF and CCF
par(mfrow=c(1,2))
acf(ricelag1, main = "ACF for diff_log(rice)")
tt = ccf(ricelag1, oillag1, ylab="CCF", main="CCF for diff_log(rice) and diff_log(oil)")
tt$lag[which(tt$acf == max(tt$acf))]  #max ccf lag = -15
tt$lag[which(tt$acf == min(tt$acf))]  #min ccf lag = -10

par(mfrow=c(1,1))

###########################################
#############presentation2#################
###########################################
library(forecast)
library(tseries)
library(strucchange)

####seasonality ~ periodogram

n = length(dat$logrice)
# fft is for free fourier transform
I = abs(fft(dat$logrice))^2/n # the periodogram
P = (4/n)*I[1:((n+1)/2)] # the scaled periodogram
f = 0:((n-1)/2)/n # frequencies
plot(f, P, type="l", xlab="Frequency", ylab="Scaled Periodogram", main = "Periodogram for log(rice)")
### frequency is too small -> period is too long (~15years)

mvspec( dat$logrice , log="no")
abline(v=2*0.0015, lty="dotted",col="red")
# this gives 188 months or 15 years
rice_fit <- lm(dat$logrice ~ time(dat$logrice))

mvspec( resid(rice_fit), log="no")
# same result

#automatic arima fitting result
fit.aic <- auto.arima(dat$logrice, ic="aic")
fit.aic
fit.bic <- auto.arima(dat$logrice, ic="bic")
fit.bic

# why not diff lag?

auto.arima(ricelag1, ic="bic")

par(mfrow=c(2,1))
acf(ricelag1, main = "ACF for diff_log(rice)")
pacf(ricelag1, main = "PACF for diff_log(rice)")
graphics.off()

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

qqnorm(fit.aic$residuals, 
       main = "ARIMA(0,1,1) fitting residual QQplot")
qqline(fit.aic$residuals)
#mean(fit.aic$residuals)
#[1] 0.001251258

##### nonlinear fitting 

##### unit root test (augmented Dickey-Fuller test)
#adf.test(dat$logrice)
#adf.test(diff(dat$logrice))

##### change in structure after regression ####
dat_modif = ts.intersect(dfrice = ts2[,1],
                         dfoil10 = lag(ts2[,2], k=10),
                         dfoil15 = lag(ts2[,2], k=15))
# ccf plot indicate diff(logrice) are correlated with diff(logoil) 
# lag 10 and 15

#lm1 = lm(dfrice~1, data= dat_modif)
#lm1
### dfrice is NOT set to zero maybe due to the inflation?

temp <- breakpoints(dfrice ~ dfoil10 + dfoil15, 
                    data = dat_modif, h = 36)
## change in structure h means minimum number of sample
# for each cluster
##(dont want (too) short term changes)
## this regression include intercept due to the inflation

#summary(temp)
coef(temp, breaks = 2) #assume two breakpoints and print where
#(Intercept)    dfoil10    dfoil15
#1987(5) - 2001(4) -0.001480141 -0.0546086 0.27978759
#2001(5) - 2008(5)  0.021944386 -0.1142983 0.05949307
#2008(6) - 2016(1) -0.010940847 -0.1166024 0.10751625

plot(ts1[,3], ylim = c(-0.6, 2.2), ylab="", 
     main = "Change in Structure")
lines(ts2[,1]+1.5, col=2)
abline(h=c(0,1.5), lty=2)
#lines(ts2[,2]+1)
#lines(fitted(temp, breaks=2)+1.5, col = 4)
lines(confint(temp, breaks = 1))

###prediction

###forcast without regression
# probably trying to predict values after breakpoint 
# using values before breakpoint
dfrice=dat_modif[,1]
qt = as.numeric(time(dat_modif))
nbreak = which(qt==(2001+4/12)) ## 2001(4) first breakpoint

armamodel1 <- arima(dat_modif[1:nbreak,1], order = c(0,0,1), 
                    include.mean = T)

forecast(arima(dat_modif[1:nbreak-1+i,1], order = c(0,0,1), 
               include.mean = T), h=1)

npred = length(dfrice) -nbreak -1
fcast1 <- forecast(armamodel1, h = npred)
plot(fcast1, ylim = c(-0.3,0.4), main = "Forecasts from MA(1)")

qq = as.numeric(dat_modif[,1])
abline(h=0, lty=2)
lines(qq, col=2)
lines(qq[1:nbreak], col=1)

###forecast after regression
armamodel2 <- arima(dat_modif[1:nbreak,1], 
                    xreg= dat_modif[1:nbreak,2:3],
                    order = c(0,0,1), include.mean = T)

fcast2 <- forecast(armamodel2, 
                   xreg = dat_modif[(nbreak+1):(nbreak+1+npred),2:3],
                   h = npred)
plot(fcast2, ylim=c(-0.3,0.4),
     main = "Forecasts from MA(1) after regression")

qq = as.numeric(dat_modif[,1])
abline(h=0, lty=2)
lines(qq, col=2)
lines(qq[1:nbreak], col=1)

###forecast after fitting different regression

fcast3 = fcast1
nbreak2 = which(qt==(2008+5/12)) ## 2008(5) second breakpoint
cc = coef(temp, breaks = 2) #assume two breakpoints and print where
tt = rep(NA, npred)

tt[1:(nbreak2-nbreak)] = cc[2,1] + 
  cc[2,2]*dat_modif[(nbreak+1):(nbreak+1+nbreak2-nbreak),2] +
  cc[2,3]*dat_modif[(nbreak+1):(nbreak+1+nbreak2-nbreak),3]

tt[(nbreak2-nbreak+1):(npred)] = cc[3,1] + 
  cc[3,2]*dat_modif[(nbreak+1+nbreak2-nbreak):(npred),2] +
  cc[3,3]*dat_modif[(nbreak+1+nbreak2-nbreak):(npred),3]

tt = ts(tt, start=170, end = 170+175-1)
tq = fcast3$mean + tt
fcast3$mean = tq
tq = fcast3$lower
tt2 = cbind(tt,tt)
fcast3$lower = fcast3$lower + tt2
fcast3$upper = fcast3$upper + tt2

plot(fcast3, ylim=c(-0.3,0.4),
     main = "Forecasts from MA(1) after piecewise regression")
qq = as.numeric(dat_modif[,1])
abline(h=0, lty=2)
lines(qq, col=2)
lines(qq[1:nbreak], col=1)

#### residuals

resid1 = qq[(nbreak+1):(nbreak+npred)] - fcast1$mean
resid2 = qq[(nbreak+1):(nbreak+npred+1)] - fcast2$mean
resid3 = qq[(nbreak+1):(nbreak+npred)] - fcast3$mean



fitarma1 = arima(resid1, order = c(0,0,1), include.mean = T)
sarima(resid1, 0, 0, 1)

fitarma2 = arima(resid2, order = c(0,0,1), include.mean = T)
sarima(resid2, 0, 0, 1)

fitarma3 = arima(resid3, order = c(0,0,1), include.mean = T)

qqnorm(fitarma1$residuals, main = "MA(1) only QQplot")
qqline(fitarma1$residuals)

qqnorm(fitarma2$residuals, main = "MA(1) after regerssion QQplot")
qqline(fitarma2$residuals)

qqnorm(fitarma3$residuals, 
       main = "MA(1) after piecewise regerssion QQplot")
qqline(fitarma3$residuals)

mean(resid1^2)
mean(resid2^2)
mean(resid3^2)

# My code

rice_prod <- read.csv("./Data/rice_prod.csv")
rice_prod <- data.frame(rice_prod)
colnames(rice_prod) <- c("Year", "Start_Stock", "End_Stock")
rice_prod <- rice_prod[1:56,]

rice_prod[,2] <- as.character(rice_prod[,2])
rice_prod$Start_Stock = gsub(",","",rice_prod$Start_Stock)
rice_prod[,2] <- as.integer(rice_prod[,2])

rice_prod[,3] <- as.character(rice_prod[,3])
rice_prod$End_Stock = gsub(",","",rice_prod$End_Stock)
rice_prod[,3] <- as.integer(rice_prod[,3])
rice_prod <- cbind(rice_prod, 1961:2016)

plot(rice_prod[,4], rice_prod$Start_Stock/1000, type="l", ylim = c(0,550),
     ylab = "Rice quantity (Billion metric ton)",
     xlab = "Year",
     main = "Starting Stock vs Ending Stock")
lines(rice_prod[,4],
      rice_prod[,3]/1000, col="red", type="l")
legend("topleft",legend=c("Starting Stock","Ending Stock"),
       col=c("Black", "Red"), lwd=2)

# use ts2[,1]

auto.arima(dat$logrice)
sarima(dfrice, 0, 0, 1)

# fit_ar2 = arima(ts2[,1], order = c(2,0,0), include.mean = T)
fit_ar2 <- sarima(dfrice, 2, 0, 0)
acf2(fit_ar2$fit$residuals^2)

# fit_ma1 = arima(dat$logrice, order = c(0,0,1), include.mean = T)
fit_ma1 <- sarima(dat$logrice, 0, 1, 1)
acf2(fit_ma1$fit$residuals)
acf2(fit_ma1$fit$residuals^2)

arch_ar2 <- garchFit(~arma(2,0)+garch(1,0), dfrice)
summary(arch_ar2)


garch_ma1 <- garchFit(~arma(0,1)+garch(1,1), dfrice)
summary(garch_ma1)
u2 <- garch_ma1@sigma.t
plot(window(dfrice, start=(2001+4/12), end=2016), 
     main="Plot of diff(log(rice))", 
     ylab="diff(log(rice))",
     xlab="Year", ylim=c(-0.5,1))
lines(window(dfrice-2*u2), lty=2, col=4)
lines(window(dfrice+2*u2), lty=2, col=4)

pred_2 <- matrix(0,1,npred)
for(i in 1:npred){
pred_2[i] <- forecast(arima(dat_modif[1:nbreak-1+i,1],
                            order = c(0,0,1), 
                            include.mean = T), 
                      h=1)$mean[1]
}

pred_2 <- ts(pred_2)



plot(ts(dfrice[nbreak:length(dfrice)]),
     ylim=c(-0.5, 0.8),
     main="Comparison with real thing and MA(1) Prediction")
lines(window(pred_2), col="red")
legend("topright",legend=c("Real Value","MA(1) Prediction"), 
       col=c("black", "red"), lwd=2)

plot(window(pred_2), 
     main="Combined prediction of MA(1) and GARCH(1,1))", 
     ylab="diff(log(rice))",
     xlab="Year", ylim=c(-0.5,0.8),
     col="red")
lines(window(pred_2[[1]] - 2*u2[(nbreak+2):length(u2)]),
      lty=2, col=4)
lines(window(pred_2[[1]] + 2*u2[(nbreak+2):length(u2)]),
      lty=2, col=4)
legend("topright",legend=c("MA(1) Prediction","GARCH(1,1) Prediction"), 
       col=c("red", "blue"), lwd=2)

plot(ts(dfrice[nbreak:length(dfrice)]),
     ylim=c(-0.5, 0.8),
     main="Real thing and GARCH(1,1)")
lines(window(pred_2[[1]] - 2*u2[(nbreak+2):length(u2)]),
      lty=2, col=4)
lines(window(pred_2[[1]] + 2*u2[(nbreak+2):length(u2)]),
      lty=2, col=4)
# lines(window(pred_2), col="red")
legend("topright",legend=c("Real Value","GARCH(1,1) Prediction"), 
       col=c("black", "blue"), lwd=2)



lines(window(dfrice+2*u2), lty=2, col=4)


# garch 2,2
# 
# 

garch2_ma1 <- garchFit(~arma(0,1)+garch(2,2), dfrice)
summary(garch2_ma1)
u2_2 <- garch2_ma1@sigma.t
plot(window(dfrice, start=(2001+4/12), end=2016), 
     main="Plot of diff(log(rice))", 
     ylab="diff(log(rice))",
     xlab="Year", ylim=c(-0.5,1))
lines(window(dfrice-2*u2_2), lty=2, col=4)
lines(window(dfrice+2*u2_2), lty=2, col=4)

pred2_2 <- matrix(0,1,npred)
for(i in 1:npred){
  pred2_2[i] <- forecast(arima(dat_modif[1:nbreak-1+i,1],
                              order = c(0,0,1), 
                              include.mean = T), 
                        h=1)$mean[1]
}

pred2_2 <- ts(pred2_2)



plot(window(pred_2), 
     main="Combined prediction of MA(1) and GARCH(2,2))", 
     ylab="diff(log(rice))",
     xlab="Year", ylim=c(-0.5,0.8),
     col="red")
lines(window(pred_2[[1]] - 2*u2_2[(nbreak+2):length(u2_2)]),
      lty=2, col=4)
lines(window(pred_2[[1]] + 2*u2_2[(nbreak+2):length(u2_2)]),
      lty=2, col=4)
legend("topright",legend=c("MA(1) Prediction","GARCH(2,2) Prediction"), 
       col=c("red", "blue"), lwd=2)

plot(ts(dfrice[nbreak:length(dfrice)]),
     ylim=c(-0.5, 0.8),
     main="Real thing and GARCH(2,2)")
lines(window(pred_2[[1]] - 2*u2_2[(nbreak+2):length(u2_2)]),
      lty=2, col=4)
lines(window(pred_2[[1]] + 2*u2_2[(nbreak+2):length(u2_2)]),
      lty=2, col=4)
# lines(window(pred_2), col="red")
legend("topright",legend=c("Real Value","GARCH(1,1) Prediction"), 
       col=c("black", "blue"), lwd=2)



lines(window(dfrice+2*u2_2), lty=2, col=4)






# fitarma_ar2ma1 = arima(ts2[,1], order = c(2,0,1), include.mean = T)
fitarma_ar2ma1 <- sarima(dfrice, 2, 0, 1)
acf2(fitarma_ar2ma1$fit$residuals)
