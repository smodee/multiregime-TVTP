rm(list=ls())
library(HMMGAS)
library(statmod)
B=200


par.true<-c(-1,1,0.5) #c(mu1, mu2, sigma2), in empirical, it should be estimated from data

##########GAS


######## 
yemp=rnorm(100,0,1)
GQ<-gauss.quad.prob(30,"normal", mu=median(yemp), sigma=sd(yemp))  ##yemp is empirical data series
weights<-GQ$weights
nodes<-GQ$nodes



dataGASCD= function (M,N) {
  par<-c(par.true,0.5,0.5,0,0,0.9,0.9)  ##this par is the parameter estimated from empirical data
  data=matrix(0,M,N)
  
  mu0=par[1];
  mu1=par[2];
  sigma2=par[3];
  sim.y<-function(par)
  {
    mu1<-par[1]
    sigma2<-par[2]
    mean.y<-mu1
    rnorm(1, mean.y,sqrt(sigma2))}
  
  par1=c(mu0,sigma2)
  par2=c(mu1,sigma2)
  
  
  omega1_LR=log(par[4]/(1-par[4]));
  omega2_LR=log(par[5]/(1-par[5]));
  
  A1=par[6];
  A2=par[7];
  
  B1=par[8];
  B2=par[9];
  
  omega1=omega1_LR;
  omega2=omega2_LR;
  
  for (i in 1: M)
    
  {
    
    eta =as.numeric(2*N)
    logLik = as.numeric(N)
    f1= as.numeric(N)
    f2= as.numeric(N)
    p1= as.numeric(N)
    p2= as.numeric(N)
    p11= as.numeric(N)
    p22= as.numeric(N)
    X_tlag= as.numeric(2*N)
    X_t= as.numeric(2*N)
    score_SCAL = as.numeric(2*N)
    I = as.numeric(N)
    y.sim=as.numeric(N)
    S=as.numeric(N)
  
  f1[1]=omega1_LR;
  f2[1]=omega2_LR;
  
  
  p11[1]=1/(1+exp(-f1[1]));
  p22[1]=1/(1+exp(-f2[1]));
  
  p1[1]=(1-p22[1])/(2-p11[1]-p22[1]);
  p2[1]=(1-p11[1])/(2-p11[1]-p22[1]);
  
  
  X_tlag[1]=p11[1]*p1[1]+(1-p22[1])*p2[1];  ##predictive prob for state 1, time t =1
  X_tlag[2]=(1-p11[1])*p1[1]+p22[1]*p2[1]; ##predictive prob for state 2, time t =1
  
  t=1
  
  S[t]<-sample(c(1,2),1,prob=c(X_tlag[1],X_tlag[2]))
  if(S[t]==1) y.sim[t]<-sim.y(par1)  ## y[1]
  if(S[t]==2) y.sim[t]<-sim.y(par2)
  
  
  eta[1]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[1]-mu0)^2/sigma2);
  eta[2]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[1]-mu1)^2/sigma2);
  
  logLik[1]=log(eta[1]*X_tlag[1]+eta[2]*X_tlag[2]);
  X_t[1]=eta[1]*X_tlag[1]/exp(logLik[1]);  ##filter prob for state 1, time t =1
  X_t[2]=eta[2]*X_tlag[2]/exp(logLik[1]); ##filter e prob for state 2, time t =1
  
  
  
  d1=as.numeric(0)
  d2=as.numeric(0)
  den=as.numeric(0)
  
  I_star=0;
  
  I[1]=0;
  
  mu_weights=as.numeric(median(yemp))
  sigma_weights=sd(yemp)
  
  for (j in 1:30) {
    
    d1=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(nodes[j]-mu0)^2/sigma2);
    
    d2=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(nodes[j]-mu1)^2/sigma2);
    
    den=d1*(p11[1]*p1[1]+(1-p22[1])*p2[1])+d2*((1-p22[1])*p1[1]+p22[1]*p2[1]);
    
    if(den ==0){
      I_star=0;	
    } else{
      
      I_star=((d1-d2)^2/den)*(2*pi*sigma_weights*sigma_weights)^0.5*exp((nodes[j]-mu_weights)^2)/(2*sigma_weights*sigma_weights);
    }
    
    I[1]=I[1]+I_star*weights[j];
  }
  
  
  
  g= as.numeric(2*N)
  
  
  S1=(eta[1]-eta[2])/exp(logLik[1]);
  g[1]=	 p1[1]*p11[1]*(1-p11[1]);
  g[2]=	-p2[1]*p22[1]*(1-p22[1]);
  
  
  g_mod=(g[1]*g[1]+g[2]*g[2])^0.5;
  
  if(g_mod==0) {
    g[1]=0;
    g[2]=0;	
  } else{
    
    g[1]=g[1]/g_mod;
    g[2]=g[2]/g_mod;	
    
  }
  
  
  
  if(I[1]==0){
    S_star=0;	
  } else{
    S_star=S1/I[1]^0.5;	
  }
  
  score_SCAL[1]= S_star*g[1];  
  score_SCAL[2]= S_star*g[2];
  
  f1[2]=omega1+ A1*score_SCAL[1]+ B1*(f1[1]-omega1);
  f2[2]=omega2+ A2*score_SCAL[2]+ B2*(f2[1]-omega2);
  
  
  p11[2]=1e-10+(1-2*1e-10)/(1+exp(-f1[2]));
  p22[2]=1e-10+(1-2*1e-10)/(1+exp(-f2[2]));
  
  
  ##loop and updating begin 

    for (t in 1:(N-2)) {
      

      
      X_tlag[t*2+1]        =p11[t+1]        *X_t[(t-1)*2+1]+    (1-p22[t+1])    *X_t[(t-1)*2+2]; ##predictive prob for state 1, time t =t
      X_tlag[(t*2)+2]    =(1-p11[t+1])    *X_t[(t-1)*2+1]+    p22[t+1]        *X_t[(t-1)*2+2]; ##predictive prob for state 2, time t =t
      
      S[t+1]<-sample(c(1,2),1,prob=c( X_tlag[t*2+1],X_tlag[(t*2)+2]))  ##
      if(S[t+1]==1) y.sim[t+1]<-sim.y(par1) ##y[t+1]
      if(S[t+1]==2) y.sim[t+1]<-sim.y(par2)
      
      
      eta[t*2+1]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[t+1]-mu0)^2/sigma2);
      eta[t*2+2]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[t+1]-mu1)^2/sigma2);
      
      logLik[t+1]=log(eta[t*2+1]*X_tlag[t*2+1]+eta[(t*2)+2]*X_tlag[(t*2)+2]);
      X_t[t*2+1]=eta[t*2+1]*X_tlag[t*2+1]/exp(logLik[t+1]);
      X_t[t*2+2]=eta[t*2+2]*X_tlag[t*2+2]/exp(logLik[t+1]);
      
      
      
      logLik[t+1]=log(eta[t*2+1]*X_tlag[t*2+1]+eta[(t*2)+2]*X_tlag[(t*2)+2]);
      X_t[t*2+2]=eta[t*2+2]*X_tlag[t*2+2]/exp(logLik[t+1]);
      X_t[t*2+1]=1-X_t[t*2+2];
      
      
      I_star=0;
      I[t+1]=0;
      
      for (j in 1:30) {
        
        d1=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(nodes[j]-mu0)^2/sigma2);
        
        d2=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(nodes[j]-mu1)^2/sigma2);
        
        den=d1*(p11[t+1]*X_t[(t-1)*2+1]+(1-p22[t+1])*X_t[(t-1)*2+2])+d2*((1-p22[t+1])*X_t[(t-1)*2+1]+p22[t+1]*X_t[(t-1)*2+2]);
        
        if(den ==0){
          I_star=0;	
        } else{
          
          I_star=((d1-d2)^2/den)*(2*pi*sigma_weights*sigma_weights)^0.5*exp((nodes[j]-mu_weights)^2)/(2*sigma_weights*sigma_weights);
        }
        
        I[t+1]=I[t+1]+I_star*weights[j];
      }
      
      
      S=(eta[t*2+1]-eta[t*2+2])/exp(logLik[t+1]);
      g[t*2+1]=	 	X_t[(t-1)*2+1]*p11[t+1]*(1-p11[t+1]);
      g[t*2+2]=	-X_t[(t-1)*2+2]*p22[t+1]*(1-p22[t+1]);
      
      
      g_mod=(g[t*2]*g[t*2]+g[t*2+1]*g[t*2+1])^ 0.5;
      
      if(g_mod==0) {
        g[t*2+1]=0;
        g[t*2+2]=0;	
      } else{
        
        g[t*2+1]=g[t*2+1]/g_mod;
        g[t*2+2]=g[t*2+2]/g_mod;	
        
      }
      
      
      
      if(I[t+1]==0){
        S_star=0;	
      } else{
        S_star=S/I[t+1]^0.5;	
      }
      
      score_SCAL[t*2+1]= S_star*g[t*2+1];   ##update scaled score
      score_SCAL[t*2+2]= S_star*g[t*2+2];
      
        
        f1[t+2]=omega1+ A1*score_SCAL[t*2+1]		+B1*(f1[t+1]-omega1);  #update f
        f2[t+2]=omega2+ A2*score_SCAL[t*2+2]	+B2*(f2[t+1]-omega2);
        
        p11[t+2]=1e-10+(1-2*1e-10)/(1+exp(-f1[t+2])); ##update transition probability
        p22[t+2]=1e-10+(1-2*1e-10)/(1+exp(-f2[t+2]));
      
    }
    
    data[i,][1:N-1]=y.sim
  }

  data
  
}


M=5  ###for the GAS model, estimation take long time, so I just set M=5 
N=1000
dataGASCDsim=dataGASCD(M,N)[,c(200:N-10)] ###dataGASCD(M,N) is simulated M path data, with each datalength N

parGASCD=matrix(0,M,9)

par=c(-1,1,0.5,0.5,0.5,0,0,0.9,0.9) ##this par is the start value in the following nlmb function,

##it can be changed to other initial value, which will not effect the estimation result as long as it does not be so far from real value

for (i in 1:M) {
  
  GAS.est<-nlminb(start= par.trasf.inv_GAS(par),filtering.single.trasf_GAS, y=dataGASCDsim[i,],B=B,control=list(eval.max=1e6,iter.max=1e6,trace=0))
  
  GAS.par.fin<-par.trasf_GAS(GAS.est$par)
  
  parGASCD[i,]=GAS.par.fin
  
}
parGASCD
colMeans(parGASCD)  ##means of 5 path simulation, each path data can give out one parameter estimation. I averge M time estimation


##to check bias, seemslow bias, which means both the datasimulation code and estimation function should be right


c(par.true,0.5,0.5,0,0,0.9,0.9)




#############TVP


dataTVPCD= function (M,N) {
  
  
  par<-c(par.true,0.5,0.5,0,0)  ##this par should be the parameter estimated from empirical data
  
  
  data=matrix(0,M,N)
  
  mu0=par[1];
  mu1=par[2];
  sigma2=par[3];
  
  omega1_LR=log(par[4]/(1-par[4]));
  omega2_LR=log(par[5]/(1-par[5]));
  
  
  A1=par[6];
  A2=par[7];
  
  omega1=omega1_LR*(1-A1);
  omega2=omega2_LR*(1-A2);
  
  
  sim.y<-function(par)
  {
    mu1<-par[1]
    sigma2<-par[2]
    mean.y<-mu1
    rnorm(1, mean.y,sqrt(sigma2))}
  
  par1=c(mu0,sigma2)
  par2=c(mu1,sigma2)
  
  
  
  for (i in 1:M) {
    
    eta =as.numeric(2*N)
    logLik = as.numeric(N)
    
    f1= as.numeric(N)
    f2= as.numeric(N)
    
    p11= as.numeric(N)
    p22= as.numeric(N)
    p1= as.numeric(N)
    p2= as.numeric(N)
    S= as.numeric(N)
    
    y.sim=as.numeric(N)
    X_tlag= as.numeric(2*N)
    X_t= as.numeric(2*N)
    
    f1[1]=omega1_LR;
    f2[1]=omega2_LR;
    
    
    p11[1]=1/(1+exp(-f1[1]));
    p22[1]=1/(1+exp(-f2[1]));
    
    p1[1]=(1-p22[1])/(2-p11[1]-p22[1]);
    p2[1]=(1-p11[1])/(2-p11[1]-p22[1]);
    
    
    X_tlag[1]=p11[1]*p1[1]+(1-p22[1])*p2[1];
    X_tlag[2]=(1-p11[1])*p1[1]+p22[1]*p2[1];
    
    t=1
    
    S[t]<-sample(c(1,2),1,prob=c( X_tlag[1], X_tlag[2]))
    if(S[t]==1) y.sim[t]<-sim.y(par1)
    if(S[t]==2) y.sim[t]<-sim.y(par2)
    
    
    eta[1]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[1]-mu0)^2/sigma2);
    eta[2]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[1]-mu1)^2/sigma2);
    
    logLik[1]=log(eta[1]*X_tlag[1]+eta[2]*X_tlag[2]);
    X_t[1]=eta[1]*X_tlag[1]/exp(logLik[1]);
    X_t[2]=eta[2]*X_tlag[2]/exp(logLik[1]);
    

    for(t in 1:(N-1)) {
      
      f1[t+1]=omega1+ A1*y.sim[t]; ##update f in the beginning here
      f2[t+1]=omega2+ A2*y.sim[t];
      
      p11[t+1]=1/(1+exp(-f1[t+1]));
      p22[t+1]=1/(1+exp(-f2[t+1]));
      
      X_tlag[t*2+1]        =p11[t+1]        *X_t[(t-1)*2+1]+    (1-p22[t+1])    *X_t[(t-1)*2+2];
      X_tlag[(t*2)+2]    =(1-p11[t+1])    *X_t[(t-1)*2+1]+    p22[t+1]        *X_t[(t-1)*2+2];
      
      
      S[t+1]<-sample(c(1,2),1,prob=c( X_tlag[t*2+1],X_tlag[(t*2)+2]))
      
      if(S[t+1]==1) y.sim[t+1]<-sim.y(par1)
      if(S[t+1]==2) y.sim[t+1]<-sim.y(par2)
      
      eta[t*2+1]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[t+1]-mu0)^2/sigma2);
      eta[t*2+2]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[t+1]-mu1)^2/sigma2);
      
      logLik[t+1]=log(eta[t*2+1]*X_tlag[t*2+1]+eta[(t*2)+2]*X_tlag[(t*2)+2]);
      X_t[t*2+1]=eta[t*2+1]*X_tlag[t*2+1]/exp(logLik[t+1]);
      X_t[t*2+2]=eta[t*2+2]*X_tlag[t*2+2]/exp(logLik[t+1]);
      
    }
    
    data[i,]=y.sim
  }
  
  data
}

M=10
N=1000
dataTVPCDsim=dataTVPCD(M,N)[,c(200:N-10)]

partvpCD=matrix(0,M,7)


par<-c(-1,1,0.5,0.5,0.5,0,0) ##again, this par is just the initial start value in follloiwng nlminb function
###you can try par setting as other values as well, as long as not far away from real value, will not effect result

for (i in 1:M){
  TVP.est <- nlminb(start = par.trasf.inv_TVP(par), filtering.single.trasf_TVP, 
                    y = dataTVPCDsim[i,], B = B, control = list(eval.max = 1e+06, iter.max = 1e+06, 
                                                              trace = 0))
  
  TVP.par.fin <- par.trasf_TVP(TVP.est$par)
  partvpCD[i,]=TVP.par.fin
}


partvpCD
colMeans(partvpCD)  ##means of 10 times simulation to check bias, seems  low bias, which means both the datasimulation code and estimation function should be right
c(par.true,0.5,0.5,0,0)







######XExo


dataexoCD= function (M,N) {
  
  
  par<-c(par.true,0.5,0.5,0,0)  ##this par is the parameter estimated from empirical data
  
  data=matrix(0,M,N)
  
  
  
  mu0=par[1];
  mu1=par[2];
  sigma2=par[3];
  
  sim.y<-function(par)
  {
    mu1<-par[1]
    sigma2<-par[2]
    mean.y<-mu1
    rnorm(1, mean.y,sqrt(sigma2))}
  
  par1=c(mu0,sigma2)
  par2=c(mu1,sigma2)
  
  omega1_LR=log(par[4]/(1-par[4]));
  omega2_LR=log(par[5]/(1-par[5]));
  
  A1=par[6];
  A2=par[7];
  
  omega1=omega1_LR*(1-A1);
  omega2=omega2_LR*(1-A2);
  
  
  for (i in 1:M) {
    eta =as.numeric(2*N)
    logLik = as.numeric(N)
  f1= as.numeric(N)
  f2= as.numeric(N)
  p11= as.numeric(N)
  p22= as.numeric(N)
  p1= as.numeric(N)
  p2= as.numeric(N)
  X_tlag= as.numeric(2*N)
  X_t= as.numeric(2*N)
  S= as.numeric(N)
  
  y.sim=as.numeric(N)
  
  f1[1]=omega1+A1*X_Exo[1]
  f2[1]=omega2+A2*X_Exo[1]
  
  
  p11[1]=1/(1+exp(-f1[1]));
  p22[1]=1/(1+exp(-f2[1]));
  
  p1[1]=(1-p22[1])/(2-p11[1]-p22[1]);
  p2[1]=(1-p11[1])/(2-p11[1]-p22[1]);
  
  
  X_tlag[1]=p11[1]*p1[1]+(1-p22[1])*p2[1];
  X_tlag[2]=(1-p11[1])*p1[1]+p22[1]*p2[1];
  

  t=1
  
  
  S[t]<-sample(c(1,2),1,prob=c( X_tlag[1], X_tlag[2]))
  if(S[t]==1) y.sim[t]<-sim.y(par1)
  if(S[t]==2) y.sim[t]<-sim.y(par2)
  
  
  
  eta[1]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[1]-mu0)^2/sigma2);
  eta[2]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[1]-mu1)^2/sigma2);
  
  logLik[1]=log(eta[1]*X_tlag[1]+eta[2]*X_tlag[2]);
  X_t[1]=eta[1]*X_tlag[1]/exp(logLik[1]);
  X_t[2]=eta[2]*X_tlag[2]/exp(logLik[1]);
  
  f1[2]=omega1+A1*X_Exo[2];
  f2[2]=omega2+A2*X_Exo[2];
  
  
  p11[2]=1/(1+exp(-f1[2]));
  p22[2]=1/(1+exp(-f2[2]));
  
  
  for (t in 1:(N-2)) {
    
    
    X_tlag[t*2+1]        =p11[t+1]        *X_t[(t-1)*2+1]+    (1-p22[t+1])    *X_t[(t-1)*2+2];
    X_tlag[(t*2)+2]    =(1-p11[t+1])    *X_t[(t-1)*2+1]+    p22[t+1]        *X_t[(t-1)*2+2];
    
    S[t+1]<-sample(c(1,2),1,prob=c( X_tlag[t*2+1],X_tlag[(t*2)+2]))  ##
    if(S[t+1]==1) y.sim[t+1]<-sim.y(par1) ##y[t+1]
    if(S[t+1]==2) y.sim[t+1]<-sim.y(par2)
    
    
    eta[t*2+1]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[t+1]-mu0)^2/sigma2);
    eta[t*2+2]=exp(-0.5*log(2*pi)-0.5*log(sigma2)-0.5*(y.sim[t+1]-mu1)^2/sigma2);
    
    logLik[t+1]=log(eta[t*2+1]*X_tlag[t*2+1]+eta[(t*2)+2]*X_tlag[(t*2)+2]);
    X_t[t*2+1]=eta[t*2+1]*X_tlag[t*2+1]/exp(logLik[t+1]);
    X_t[t*2+2]=eta[t*2+2]*X_tlag[t*2+2]/exp(logLik[t+1]);
    

      
      f1[t+2]=omega1+ A1*X_Exo[t+2];
      f2[t+2]=omega2+ A2*X_Exo[t+2];
      
      p11[t+2]=1/(1+exp(-f1[t+2]));
      p22[t+2]=1/(1+exp(-f2[t+2]));
      
    
  }
  
  data[i,][1:N-1]=y.sim
  }
  
  data
  
}


M=10
N=1000

X_Exo=rnorm(N)

dataexoCDsim=dataexoCD(M,N)[,c(200:N-10)]


par<-c(-0.8,0.8,0.4,0.6,0.6,0,0) ##fexample, here I change the fouth and fifth start parameter as 0.6, but it will not effect result
parexoCD=matrix(0,M,7)

for (i in 1:M) {

  TVP.XExo.est<-nlminb(start= par.trasf.inv_TVPXExo(par),filtering.single.trasf_TVPXExo, y=dataexoCDsim[i,],B=B, X_Exo =X_Exo,control=list(eval.max=1e6,iter.max=1e6,trace=0))
  
  TVP.XExo.par.fin<-par.trasf_TVPXExo(TVP.XExo.est$par)[1:7]
  
  parexoCD[i,]=TVP.XExo.par.fin
  
}

parexoCD
colMeans(parexoCD)

c(par.true,0.5,0.5,0,0)

