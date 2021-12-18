## Copyright (c) 2020 Ian Hall
## See LICENCE for licensing information

## Growth rate estimates for confirmed cases in Europe and for different metrics in Italy using GAM
## Figure 1 (main text) and figures S1 and S2 (electronic supplementary material) of:
## 
## Pellis L, Scarabel F, Stage HB, Overton CE, Chappell LHK, Fearon E, Bennett E, 
## University of Manchester COVID-19 Modelling Group, Lythgoe KA, House TA and Hall I, 
## "Challenges in control of COVID-19: short doubling time and long delay to effect of interventions", 
## Philosophical Transactions of the Royal Society B (2021)


DoubleTime <- function(dat, timev, npts=200, meth="GCV.Cp", FE='None', plt=FALSE, thck=1.5, figtitle="", origin="2020-01-09"){
  
  res <- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts),doub=rep(0,npts),doubup=rep(0,npts),doublow=rep(0,npts))
  Tv <- timev
  
  if(FE=='None'){
    MGAM <- gam(dat~s(Tv), family=quasipoisson, method=meth)
  }else{
    DW <- weekdays(as.Date(Tv, origin = origin))
    if(FE=='WE'){
      DW <- ifelse(DW=='Sunday','WE',ifelse(DW=='Saturday','WE','WD'))
    }
    MGAM <- gam(dat~s(Tv)+DW, family=quasipoisson, method=meth)
  }
  
  xv<-seq(min(Tv),max(Tv), length=npts)
  if(FE=='None'){
    newd <- data.frame(Tv=xv)
  }else{
    dow <- weekdays(as.Date(xv, origin = origin))
    if(FE=='WE'){
      dow <- ifelse(dow=='Sunday','WE',ifelse(dow=='Saturday','WE','WD'))
    }
    newd <- data.frame(Tv=xv, DW=dow)
  }
  X0 <- predict(MGAM, newd,type="lpmatrix")
  
  eps <- 1e-7 ## finite difference interval
  xv<-seq(min(Tv),max(Tv), length=npts)+eps
  if(FE=='None'){
    newd <- data.frame(Tv=xv)
  }else{
    dow <- weekdays(as.Date(xv, origin = origin))
    if(FE=='WE'){
      dow <- ifelse(dow=='Sunday','WE',ifelse(dow=='Saturday','WE','WD'))
    }
    newd <- data.frame(Tv=xv, DW=dow)
  }
  X1 <- predict(MGAM, newd,type="lpmatrix")
  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
  
  off <- ifelse(FE=='None',1,ifelse(FE=='WE',2,7))  
  Xi <- Xp*0 
  Xi[,1:9+off] <- Xp[,1:9+off] ## weekend Xi%*%coef(MGAM) = smooth deriv i
  df <- Xi%*%coef(MGAM)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  ## derivative calculation, pers comm S. N. Wood, found in mgcv:  Mixed  GAM  Computation  Vehicle  with  automatic  smoothness  estimation.  R  packageversion 1.8-31 (2019) https://CRAN.R-project.org/package=mgcv.
  
  res$sdt <- df
  res$sdtup <- df+2*df.sd
  res$sdtlow <- df-2*df.sd
  res$doub <- ifelse(res$sdt < 0, 100, log(2)/res$sdt)
  res$doubup <- ifelse(res$sdtup < 0, 100, log(2)/res$sdtup)
  res$doublow <- ifelse(res$sdtlow < 0, 100, log(2)/res$sdtlow)
  
  MGLM <- glm(dat~(Tv), family=quasipoisson)
  Doubling<-c(MGLM$coefficients[2],confint(MGLM)[2,1],confint(MGLM)[2,2])
  
  if(plt==TRUE){
    par(mar = c(5,4,4,4) + 0.1)
    par(mfrow=c(1,2))
    plot(as.Date(xv, origin = origin),df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)), ylab='Instantaneous growth rate', xlab='Time', main=figtitle, lwd=2*thck)
    lines(as.Date(xv, origin = origin),df+2*df.sd,lty=2, lwd=thck);
    lines(as.Date(xv, origin = origin),df-2*df.sd,lty=2, lwd=thck)
    abline(h=0, col=4)
    text(as.Date(43930, origin = origin),max(df),ifelse(df[length(df)]-2*df.sd[length(df)]>0,'Growth',ifelse(df[length(df)]+2*df.sd[length(df)]<0,'Decay','Plateau')))
    axis(4,at=c(-log(2)/7, 0,log(2)/7,log(2)/4,log(2)/3,log(2)/2), labels=c(-7,'Infinite',7,4,3,2))
    mtext("Doubling time", side = 4, line = 3)
    
    plot(as.Date(Tv, origin = origin), dat, main='Fit compared with model', ylab='Number', xlab='Time', pch=16, col=4)
    lines(as.Date(Tv, origin = origin), fitted(MGAM), lwd=2*thck)
    
    p <- predict(MGAM, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGAM$family$linkinv(upr)
    lwr <- MGAM$family$linkinv(lwr)
    lines(as.Date(Tv, origin = origin), upr, col=1, lty=2, lwd=thck)
    lines(as.Date(Tv, origin = origin), lwr, col=1, lty=2, lwd=thck)
    
  }
  
  res
}

library(mgcv)
