library(matrixcalc)
library(copula)
library(TMB)
library("VGAM")
require(corrplot)
library(glmmTMB)
library(optimr)
#import data
Data = read.csv("anonymous.csv",sep = ";", header = TRUE)
# for a new dataset, where we count the number of bankpuptsies
YearFiled <- rep(1980:2017,each=12)
MonthFiled <- rep(1:12,38)
blue <- c("WA","OR","CA","NV","NM","MN","IA","MO","IL","WI","MI","OH","PA","MD","DE","NJ","NY","CT","RI","MA","VT","NH","ME","FL","HI")
red <- c("ID","UT","AZ","MT","MY","CO","TX","ND","SD","NE","KS","OK","AR","LA","MS","AL","IN","KY","TN","GA","SC","NC","VA","MV","AK")

#creating new variable blue, if state blue 1 if red 0
Data$HeadStAtFiling <- as.character(Data$HeadStAtFiling)
funkblue <- function(x){
  return(ifelse(x%in%blue,1,0))
}
newvar <- sapply(Data$HeadStAtFiling,funkblue)
Data$blue <- newvar
data <- Data[order(Data$YearFiled),]
# new counting of defaults, now based on blue or red
blueStates <- numeric(length(YearFiled))
redStates <- numeric(length(YearFiled))
for(i in 1:nrow(Data)){
  year <- YearFiled[i]
  month <- MonthFiled[i]
  for(j in 1:nrow(Data)){
    year2 <- Data$YearFiled[j]
    State <- Data$blue[j]
    month2 <- Data$MonthFiled[j]
    if(year==year2 && month==month2 && State==1){
      blueStates[i] <- blueStates[i]+1
    }
    else if(year==year2 && month==month2 && State==0){
      redStates[i] <- redStates[i]+1
    }
  }
}
redStates <- redStates[1:(length(redStates)-3)]
blueStates <- blueStates[1:(length(blueStates)-3)]
#The desired bivariate count time series is the following
def_st<- as.matrix(cbind(redStates,blueStates))
#---------Now that we all the data loaded in, we proceed fiting various models--------------
#----------more effectiv to use tmb objects as follows----------

#compile the tmb functions, skip if these have been compiled before

compile("linear.cpp") 
compile("log_linear.cpp")  
compile("log_lin_cov.cpp")
compile("bivpoisson.cpp")
compile("bivpoissonlog.cpp")
compile("defaultTmb.cpp")
compile("covariateModel.cpp")
compile("lin_cov.cpp")

# load the a bove dynamically 
dyn.load(dynlib("linear"))
dyn.load(dynlib("log_linear"))
dyn.load(dynlib("log_lin_cov"))
dyn.load(dynlib("lin_cov"))
dyn.load(dynlib("bivpoisson"))
dyn.load(dynlib("bivpoissonlog"))
dyn.load(dynlib("defaultTmb"))
#dyn.load(dynlib("covariateModel"))

#function to fit univariate model with covariates
# given data and covarite data,
calculate.parlog_linear.cov <- function(data,covdata) ##BS: changed from lambda to ny in version_2
{
  size <- length(data)
  m1 = glmmTMB(data[2:size] ~ log(data[1:(size - 1)] + 1)+covdata[1:(size-1)], family=poisson,zi=~covdata[1:(size-1)]+data[1:(size - 1)])
  start <- coef(summary(m1))[[1]][1:3]
  #p = c(start[1],0, start[2],start[3])
  #start = glm(data[2:size] ~ log(data[1:(size - 1)] + 1)+covdata[1:(size-1)], family = poisson)$coef
  dat <- list(y=data,x=covdata)
  parm <- list(d=start[1],a=0,b=start[2],c=start[3])
  umodcov <- MakeADFun(dat,parm,DLL="log_lin_cov",silent = TRUE) 
  #results2 = optim(umodcov$par,umodcov$fn, umodcov$gn)$par
  result = optim(umodcov$par,umodcov$fn, umodcov$gn)
  results2 =result$par
  
  ####Calculation of \hat{\lambda}
  ny=rep(NA,size)
  ny[1]=log(mean(data))
  
  #lambda[1] = log(lambda[1]) ##Try this! Currently, that means lambda really is ny
  #ny=rep(NA, length(data.test))
  #ny[1]=mean(data.test)
  
  for (t in 2:size)
  {
    ny[t]=results2[1] + results2[2]*ny[t-1] + results2[3]*log(data[t-1]+1) + results2[4]*covdata[t]
  }
  #lambda <- exp(lambda)
  
  return(list(results2, ny,result))
}

#Function to deliver starting values for bivariate- cov-model,
#Given data and covariate data
StartingParam_nondiag.cov <- function(defdata,covdata) ##here is the change from Kostas.r file!!! 
{
  data <- defdata
  step1=calculate.parlog_linear.cov(data[,1],covdata) #denne funksjonen er gitt rett nedenfor
  step2=calculate.parlog_linear.cov(data[,2],covdata)
  pr.par1=step1[[1]]
  pr.par2=step2[[1]]
  ny1=step1[[2]]
  ny2=step2[[2]]
  ny=cbind(ny1,ny2)
  size <- nrow(defdata)
  outlm1=vglm(data[2:size,] ~ny[1:(size-1),]+log(data[1:(size-1),]+1)+covdata[2:size], poissonff)
  hatd=coef(outlm1)[1:2] #This is even better!!!
  hatA=coef(outlm1)[3:6]
  hatB=coef(outlm1)[7:10]
  hatC=coef(outlm1)[11:12]
  p0=c(hatd, hatA, hatB,hatC) #startingvalues
  # TMB needs to recognize the parameters
  p0=c(p0[[1]],p0[[2]],p0[[3]],p0[[4]],p0[[5]],p0[[6]],p0[[7]],p0[[8]],p0[[9]],p0[[10]],p0[[11]],p0[[12]])
  par0 <- list(d1=p0[1],d2=p0[2],a11=p0[3],a12=p0[4],a21=p0[5],a22=p0[6],b11=p0[7],b12=p0[8],b21=p0[9],b22=p0[10],c1=p0[11],c2=p0[12])
  return(par0)
}
#--------------Now the linear model---------------------
#Function to fit univariate linear model
# uses TMB
calculate.parlinear <- function(data)
{ 
  ###Estimation
  data.test=data
  r1=arima(data.test, order=c(1,0,1), method="CSS", optim.control=list(maxit=100))  
  theta=max(min(r1$coef[2],-0.001),-0.99)
  phi=max(r1$coef[1],0.001-theta)
  mu=sigma2=max(r1$coef[3],0.01)
  start=c(mu*(1-phi), sigma2*(1-phi), -theta, (phi+theta))
  dat <- list(y=data.test)
  parm <- list(d=start[2],a=start[3],b=start[4])
  #parm <- c(start[2],a=start[3],start[4])
  umodlin <- MakeADFun(dat,parm,DLL="linear",silent = TRUE) 
  result = optim(umodlin$par,umodlin$fn, umodlin$gn)
  rprt <- sdreport(umodlin,par.fixed = result$par)
  SE <- summary(rprt,"report")
  #result <- optim(parm,liklinear.poisson, scorelinear.poisson,data=data.test)
  results2 <- result$par
  #results2=constrOptim(theta=start[2:4], f=liklinear.poisson , grad=scorelinear.poisson,  data=data.test, uinew, cinew, outer.iterations = 100, outer.eps = 1e-05, method="BFGS")$par
  ####Calculation of \hat{\lambda}
  lambda=rep(NA, length(data.test))
  lambda[1]=mean(data.test)
  for (t in 2:length(data.test))
  {
    lambda[t]=results2[1]+results2[2]*lambda[t-1]+results2[3]*data.test[t-1]  ##BS: loglinear
  }
  return(list(results2,lambda,result,SE))
}

#Function to obtain starting values for the bivariate linear model
# Given the data(2-columns)
StartingParam_nondiag <- function(data) 
{
  size <- nrow(data)
  step1=calculate.parlinear(data[,1])
  step2=calculate.parlinear(data[,2])
  pr.par1=step1[[1]]
  pr.par2=step2[[1]]
  lambda1=step1[[2]]
  lambda2=step2[[2]]
  lambda=cbind(lambda1,lambda2)
  outlm1=lm(data[2:size,] ~lambda[1:(size-1),]+data[1:(size-1),]) ##Starting value from here!
  hatd=c(pr.par1[1], pr.par2[1])
  hatA=coef(outlm1)[2:3,]
  hatB=coef(outlm1)[4:5,]
  p0=c(hatd, vec(hatA), vec(hatB))
  p0 <- list(d1=p0[1],d2=p0[2],a11=p0[3],a12=p0[4],a21=p0[5],a22=p0[6],b11=p0[7],b12=p0[8],b21=p0[9],b22=p0[10])
  return(list(p0,lambda))
}
# Function to fit a univariate log-linear model
calculate.parloglinear <- function(data) ##BS: changed from lambda to ny in version_2
{
  ###Estimation
  data.test=data
  ##BS: new startvalues!
  #size=length(data.test)
  #m1 = glmmTMB(data[2:size] ~ log(data[1:(size - 1)] + 1), family = poisson(link="log"),ziformula = ~.)
  #start <- coef(summary(m1))[[1]][1:3]
  start = glm(data.test[2:length(data.test)] ~ log(data.test[1:(length(data.test) - 1)] + 1), family = poisson)$coef
  y <- list(y=data.test)
  parm <- list(d=start[[1]],a=0,b=start[[2]])
  #parm <- c(start[[1]],0,start[[2]])
  unimod <- MakeADFun(y,parm,DLL="log_linear",silent = TRUE) 
  result <- optim(unimod$par,unimod$fn,unimod$gr,method = "BFGS")
  #result <- optim(parm,lik_loglinear.poisson,score_loglinear.poisson,data=data.test,method = "BFGS")
  results2<- result$par
  rprt2 <- sdreport(unimod,par.fixed = results2)
  SE <- summary(rprt2,"report")
  
  #results2 = optim(lik_loglinear.poisson, p = c(start[1],0, start[2]), score_loglinear.poisson, data = data.test,method = "BFGS")$par
  
  ####Calculation of \hat{\lambda}
  ny=rep(NA, length(data.test))
  ny[1]=log(mean(data.test))
  
  #lambda[1] = log(lambda[1]) ##Try this! Currently, that means lambda really is ny
  #ny=rep(NA, length(data.test))
  #ny[1]=mean(data.test)
  
  for (t in 2:length(data.test))
  {
    ny[t]=results2[3]+results2[1]*ny[t-1]+results2[2]*log(data.test[t-1]+1)
  }
  #lambda <- exp(lambda)
  
  return(list(results2, ny,result,SE))
}

#Function to obtain starting values for the bivariate log-linear model
#Given data (2-columns)
StartingParam_nondiag_log <- function(data) ##here is the change from Kostas.r file!!! 
{
  size = nrow(data)
  step1=calculate.parloglinear(data[,1]) #denne funksjonen er gitt rett nedenfor
  step2=calculate.parloglinear(data[,2])
  pr.par1=step1[[1]]
  pr.par2=step2[[1]]
  ny1=step1[[2]]
  ny2=step2[[2]]
  ny=cbind(ny1,ny2)
  ny= replace(data, is.na(ny), 1) # in case NAN produced
  outlm1=vglm(data[2:size,] ~ny[1:(size-1),]+log(data[1:(size-1),]+1), poissonff)
  hatd=coef(outlm1)[1:2] #This is even better!!!
  hatA=coef(outlm1)[3:6]
  hatB=coef(outlm1)[7:10]
  p0_log=c(hatd, hatA, hatB) #startingvalues
  param_log <- list(d1=p0_log[[1]],d2=p0_log[[2]],a11=p0_log[[3]],a12=p0_log[[4]],
                    a21=p0_log[[5]],a22=p0_log[[6]],b11=p0_log[[7]],b12=p0_log[[8]],
                    b21=p0_log[[9]],b22=p0_log[[10]])
  return(list(param_log,data))
}
#_________________________________________________________________________________________________
#function to calculate in-sample predictoin and RMSE, provided estimated theta 
#corresponding to each univariate fit
predrmse <- function(par,mod,covar=0){
  pred <- numeric(length(redStates))
  pred[1] <- 0
  if(length(par)>3){
    for(t in 2:length(redStates)){
      pred[t] <- par[1] +par[2]*pred[t-1]+par[3]*log(counts$freq[t-1]+1) + par[4]*covar[t]
    }
    RMSE <- sqrt(sum((counts$freq-exp(pred))^2))
    return(list(pred=pred,rmse=RMSE))
  }else if(mod=="log"){
    for(t in 2:length(redStates)){
      pred[t] <- par[1] +par[2]*pred[t-1]+par[3]*log(counts$freq[t-1]+1)
    }
    RMSE <- sqrt(sum((counts$freq-exp(pred))^2))
    return(list(pred=pred,rmse=RMSE))
  }else if(mod=="lin"){
    for(t in 2:length(redStates)){
      pred[t] <- par[1] +par[2]*pred[t-1]+par[3]*counts$freq[t-1]
    }
    RMSE <- sqrt(sum((counts$freq-pred)^2))
    return(list(pred=pred,rmse=RMSE))
  }else{
    print("please check function aruments")
  }
}
#______________________________Multivariate model selection____________________

#Function to fit a given bivarite model
bivmodel <- function(data=def_st,covdata,model){
  if(model=="cov"){
    p <- StartingParam_nondiag.cov(data,covdata)
    dat <- list(Y=data,X=covdata)
    mod <- MakeADFun(dat,p,DLL="defaultTmb",silent = TRUE)
    fit <- optim(mod$par,mod$fn,mod$gr,method = "BFGS")
    #adrep.cov <- sdreport(mod,par.fixed = fit)
    #summary(adrep.cov,"report")
    return(fit$par) # NB! change for out-of sample pred
  }else if(model=="log"){
    par0 <- StartingParam_nondiag_log(data)[[1]]
    Y <- list(Y=data)
    model_log <- MakeADFun(Y,par0,DLL="bivpoissonlog",silent = TRUE)
    fit_log <- optim(model_log$par,model_log$fn,model_log$gr,method = "BFGS")
    #adrep2 <- sdreport(model_log, par.fixed = fit_log$par) 
    #summary(adrep2, "report")
    return(fit_log$par)
  }else if(model=="lin"){
    param <- StartingParam_nondiag(data)[[1]]
    Y <- list(Y=data)
    model_lin <- MakeADFun(Y,param,DLL="bivpoisson",silent = TRUE)
    fitlin <- optim(param, model_lin$fn, model_lin$gr,method = "BFGS") #
    #adrep <- sdreport(model_lin, par.fixed = fitlin$par) 
    #summary(adrep, "report") # estimates and std.errors
    return(fitlin$par)
  }
}

#modres(bivmodel(model = "lin"))
#modres(bivmodel(model = "log"))
#bivpred(bivmodel(model = "lin"),mod = "lin")

# function to calculate bivariate in-sample prediction and the corresponding RMSE
bivpred <- function(parvec,cov,mod){
  pred2 <- matrix(0, ncol = length(blueStates), nrow = 2)
  d=matrix(parvec[1:2], byrow=T, nrow=2)
  A=matrix(parvec[3:6], byrow=T, nrow=2)
  B=matrix(parvec[7:10], byrow=T, nrow=2)
  C=matrix(parvec[11:12], byrow=T, nrow=2)
  pred2[,1]<- c(1,1)
  #in-sample-predictions for the log-model including covariates
  if(mod=="cov"){
    for(t in 2:length(blueStates)){
      pred2[,t] <- d + A%*%pred2[,t-1] + B%*%(log(t(def_st)[,t-1]+1)) + C%*%cov[t]
    }
    pred2 <- t(exp(pred2))
    RMSE <- sqrt(sum((def_st[,1]-pred2[,1])^2+(def_st[,2]-pred2[,2])^2))
    return(list(pred=pred2,rmse=RMSE))
  }
  if(mod=="log"){
    for(t in 2:length(blueStates)){
      pred2[,t] <- d + A%*%pred2[,t-1] + B%*%(log(t(def_st)[,t-1]+1)) 
    }
    pred2 <- t(exp(pred2))
    RMSE <- sqrt(sum((def_st[,1]-pred2[,1])^2+(def_st[,2]-pred2[,2])^2))
    return(list(pred=pred2,rmse=RMSE))
  }
  if(mod=="lin"){
    for(t in 2:length(blueStates)){
      pred2[,t] <- d + A%*%pred2[,t-1] + B%*%(t(def_st)[,t-1]) 
    }
    pred2 <- t(pred2)
    RMSE <- sqrt(sum((def_st[,1]-pred2[,1])^2+(def_st[,2]-pred2[,2])^2))
    return(list(pred=pred2,rmse=RMSE))
  }
}

uni_selection2 <- vector("list",ncol(covar)-1)
names(covar)
for(i in 2:ncol(covar)){
  print(i)
  mod <- bivmodel(def_st,covdata=covar[,i],model = "cov")
  #variable=names(covar)[i]
  #aic=modres(mod)[[1]]
  #bic=modres(mod)[[2]]
  #rmse <- bivpred(mod$par,covar[,i])[[2]]
  #convergence=mod$convergence
  #uni_selection2[[i]] <- c(variable,aic,bic,rmse,convergence,mod$par)
  uni_selection2[[i]] <- mod
}
uni_selection2
#------------------------out-of-sample prediction------------------------------
m <- 269
len <- length(redStates)-1
#first rolling fit of the model with cov sp500ret
est.cov <- vector("list", 453)
for(k in (len+1):m){
  print(k)
  est.cov[[k]] <- bivmodel(data = def_st[1:k,],covar$SP500ret[1:k],model = "cov")
}

# Rolling fit for the linear model
#NB it can only go back as far as 271, m needs to be changed but since 
# log model with and with-out cov are better fit we proceed with those, settin m=260
est.lin <- vector("list", 453)
for(k in (len+1):m){
  print(k)
  est.lin[[k]] <- bivmodel(data = def_st[1:k,],model = "lin")
  print(est.lin[[k]])
}
#try <- StartingParam_nondiag_log(def_st[1:451,])[[1]]
est.log <- vector("list", 453)
for(k in (len+1):m){
  print(k)
  est.log[[k]] <- bivmodel(data = def_st[1:k,],model = "log")
}

# Now that we have the estimates of theta, we can calculate the one-step a head forcast
# such that lambda_{t+1} <- lambda_{t+1}(theta_{t})
forcast.cov <- matrix(0, ncol = nrow(def_st), nrow = 2)
forcast.log <- matrix(0, ncol = nrow(def_st), nrow = 2)

#forcast.lin <- matrix(0, ncol = nrow(def_st), nrow = 2)
#function to extract parameter matrices 
getpar <- function(x){
  if(length(x)>10){
    d=matrix(x[1:2], byrow=T, nrow=2)
    A=matrix(x[3:6], byrow=T, nrow=2)
    B=matrix(x[7:10], byrow=T, nrow=2)
    C=matrix(x[11:12], byrow=T, nrow=2)
    return(list(d=d,A=A,B=B,C=C))
  }
  else{
    d=matrix(x[1:2], byrow=T, nrow=2)
    A=matrix(x[3:6], byrow=T, nrow=2)
    B=matrix(x[7:10], byrow=T, nrow=2)
    return(list(d=d,A=A,B=B))
  }
}
#test getpar function
#try <- getpar(est.cov[[434]]) 
forcast.cov[,m] <- c(0,0)
forcast.log[,m] <- c(1,1)
#forcast.lin[,m] <- c(0,0)
for(t in m:len){
  print(t)
  p1 <- getpar(est.cov[[t]])
  p2 <- getpar(est.log[[t]])
  #p3 <-  getpar(est.lin[[t]])
  forcast.cov[,t+1] <- p1$d + p1$A%*%forcast.cov[,t] +p1$B%*%(log(t(def_st)[,t]+1)) +(p1$C)%*%(covar$SP500ret[t+1])
  forcast.log[,t+1] <- p2$d + p2$A%*%forcast.log[,t] +p2$B%*%(log(t(def_st)[,t]+1))
  # forcast.lin[,t+1] <- p3$d + p3$A%*%forcast.lin[,t] +p3$B%*%(t(def_st)[,t])
}
forcast.cov <- t(exp(forcast.cov))
forcast.log <- t(exp(forcast.log))
#forcast.lin <- t(forcast.lin)
par(mfrow=c(1,1))
plot(blueStates,type="l",col="grey")
lines(forcast.cov[,2],col="green")
lines(forcast.log[,2],col="red")
#lines(forcast.lin[,2],col="blue")
# function to calculate prediction errors
errfunc <- function(predicted,err){
  MSFE <- rep(0,len+1)
  FS <- rep(0,len+1)
  if(err=="MSFE"){
    for(t in (m+1):(len)){
      sum1 <- 0
      for(s in m:t){
        sum1 <- sum1+(redStates[s+1]-predicted[,1][s+1])^2 +(blueStates[s+1]-predicted[,2][s+1])^2
      }
      MSFE[t] <- 1/(t-m)*sum1
    }
    return(MSFE)
  }
  else if(err=="FS"){
    for(t in (m+1):(len)){
      sum2 <- 0
      for(s in m:t){
        sum2 <- sum2+(redStates[s+1]*log(predicted[,1][s+1])-predicted[,1][s+1])^2 +(blueStates[s+1]*log(predicted[,2][s+1])-predicted[,2][s+1])^2
      }
      FS[t] <- 1/(t-m)*sum2
    }
    return(FS)
  }
  else{
    print("error measure not well specified")
  }
}
msfe.cov <- errfunc(forcast.cov,"MSFE")[(m+1):len]
msfe.log <- errfunc(forcast.log,"MSFE")[(m+1):len]
#msfe.lin <- errfunc(forcast.lin,"MSFE")[(m+1):len]
fs.cov <- errfunc(forcast.cov,"FS")[(m+1):len]
fs.log <- errfunc(forcast.log,"FS")[(m+1):len]
fs.lin <- errfunc(forcast.lin,"FS")[(m+1):len]

pdf("oospred.pdf",width = 9,height = 6)

par(mfrow=c(2,1),mar=c(1,4,2,1),omi=c(0.4,0,0,1))
plot(msfe.cov,type="l",col="green",ylab="MSFE",xlab="Time",xaxt = "n",cex.lab=1,cex.axis=1)
axis(1, at=c(8,68,128), labels=c("2003","2008","2012"),cex.axis=1)
lines(msfe.log,col="red")
lines(msfe.lin,col="blue")


plot(fs.cov,type="l",col="green",xaxt = "n",ylab="FS",xlab="Time",cex.lab=1,cex.axis=1,ylim=c(-1,25))
axis(1, at=c(8,68,128), labels=c("2003","2008","2012"),cex.axis=1)
lines(fs.log,col="red")
dev.off()
lines(fs.lin,col="blue")
#ggplot graphics for the above
p1 <- autoplot(ts( cbind(msfe.cov,msfe.log)  , start = c(2003,7), frequency = 12 ),facets = FALSE) + labs(color="models")
p2 <- autoplot(ts( cbind(fs.cov,fs.log)  , start = c(2003,7), frequency = 12 ),facets = F) + labs(color="models")
p3 <- grid.arrange(p1,p2,ncol=1)
ggsave("oospred.pdf",plot = p3,width = 8,height = 6,units = "cm")
#------------------Cummulative periodogram for each model-------------
#function to plot cummulative periodograms
library(ggfortify)
library(ggplot2)
library(gridExtra)
cumpgram <- function(pred){
  residred <- (redStates-pred[,1])/(sqrt(pred[,1]))
  residblue <- (blueStates-pred[,2])/(sqrt(pred[,2]))
  p1 <- ggcpgram(residred)
  p2 <- ggcpgram(residblue)
  residred <- as.ts(residred)
  residblue <- as.ts(residblue)
  p3 <- autoplot(acf(residred,plot = FALSE,lag.max = 60))
  p4 <- autoplot(acf(residblue,plot = FALSE,lag.max = 60))
  grid.arrange(p1,p2,p3,p4,ncol=2)
}


