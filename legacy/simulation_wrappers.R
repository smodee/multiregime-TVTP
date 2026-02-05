library(HMMGAS)
library(statmod)


# Source required helper functions from new implementation
source("helpers/utility_functions.R")
source("helpers/transition_helpers.R") 
source("helpers/parameter_transforms.R")



##########GAS



dataGASCD_original <- function(M, N, mu, sigma2, init_trans, A, B, weights, nodes, mu_weights, sigma_weights) {
  
  data=matrix(0,M,N)
  
  # Convert to original internal format:
  mu0 = mu[1]; mu1 = mu[2]
  p11 = 1 - init_trans[1]; p22 = 1 - init_trans[2]  
  A1 = A[1]; A2 = A[2]
  B1 = B[1]; B2 = B[2]  # Key addition: B parameters for persistence
  
  # GAS model uses different omega calculation than TVP/Exogenous
  omega1_LR = log(p11/(1-p11))
  omega2_LR = log(p22/(1-p22))
  
  # For GAS model: omega = omega_LR (no A multiplication like in TVP/Exo)
  omega1 = omega1_LR
  omega2 = omega2_LR
  
  sim.y<-function(par)
  {
    mu1<-par[1]
    sigma2<-par[2]
    mean.y<-mu1
    rnorm(1, mean.y,sqrt(sigma2))}
  
  par1=c(mu0,sigma2)
  par2=c(mu1,sigma2)
  
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


run_original_GAS_estimation <- function(M, N, mu, sigma2, init_trans, A, B, B_burnin = 0, C_cutoff = 0) {
  
  # STEP 1: Set up initial quadrature (use regime parameters for representative data)
  # Create representative sample from regime parameters
  regime_sample <- c(rnorm(500, mu[1], sqrt(sigma2)), rnorm(500, mu[2], sqrt(sigma2)))
  sample_median <- median(regime_sample)
  sample_sd <- sd(regime_sample)
  GQ <- gauss.quad.prob(30, "normal", mu = sample_median, sigma = sample_sd)
  assign("weights", GQ$weights, envir = .GlobalEnv)
  assign("nodes", GQ$nodes, envir = .GlobalEnv)
  
  # STEP 2: Now generate data with proper quadrature setup
  dataGASCDsim <- dataGASCD_original(M, N, mu, sigma2, init_trans, A, B, 
                                     weights, nodes, sample_median, sample_sd)
  
  # STEP 3: Apply burn-in and cutoff if needed
  if (B_burnin > 0 || C_cutoff > 0) {
    start_idx <- B_burnin + 1
    end_idx <- N - C_cutoff
    dataGASCDsim <- dataGASCDsim[, start_idx:end_idx]
  }
  
  # Results matrix (9 parameters for GAS model)
  parGASCD <- matrix(0, M, 9)
  
  # Estimation loop
  for (i in 1:M) {
    
    # Generate starting point using the new implementation's method
    y_for_start <- dataGASCDsim[i,]
    starting_points <- generate_starting_points(y_for_start, K = 2, model_type = "gas", n_starts = 1, seed = i)
    par_start_new_format <- starting_points[[1]]
    
    # Convert from new implementation format to original format
    # New format: [mu1, mu2, sigma2_1, sigma2_2, init_trans, A, B]  
    # Original format: [mu1, mu2, sigma2, p11, p22, A1, A2, B1, B2]
    mu_start <- par_start_new_format[1:2]
    sigma2_start <- par_start_new_format[3]  # Use first variance
    init_trans_start <- par_start_new_format[5:6]  # off-diagonal transitions
    A_start <- par_start_new_format[7:8]
    B_start <- par_start_new_format[9:10]  # B parameters
    
    # Convert transitions: new format -> original format
    p11_start <- 1 - init_trans_start[1]
    p22_start <- 1 - init_trans_start[2]
    
    par_start_original <- c(mu_start[1], mu_start[2], sigma2_start, p11_start, p22_start, 
                            A_start[1], A_start[2], B_start[1], B_start[2])
    
    GAS.est <- nlminb(
      start = par.trasf.inv_GAS(par_start_original), 
      filtering.single.trasf_GAS, 
      y = dataGASCDsim[i,], 
      B = B_burnin,
      control = list(eval.max = 1e+06, iter.max = 1e+06, trace = 0)
    )
    
    GAS.par.fin <- par.trasf_GAS(GAS.est$par)
    parGASCD[i,] <- GAS.par.fin
  }
  
  return(parGASCD)
}




#############TVP


dataTVPCD_original <- function(M, N, mu, sigma2, init_trans, A) {
  
  data=matrix(0,M,N)
  
  # Convert to original internal format:
  mu0 = mu[1]; mu1 = mu[2]
  p11 = 1 - init_trans[1]; p22 = 1 - init_trans[2]  
  A1 = A[1]; A2 = A[2]
  
  omega1_LR = log(p11/(1-p11))
  omega2_LR = log(p22/(1-p22))
  
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

run_original_TVP_estimation <- function(M, N, mu, sigma2, init_trans, A, B_burnin = 0, C_cutoff = 0) {
  
  # Generate data using our new parameterized function (full series)
  dataTVPCDsim <- dataTVPCD_original(M, N, mu, sigma2, init_trans, A)
  
  # Apply burn-in and cutoff if needed
  if (B_burnin > 0 || C_cutoff > 0) {
    start_idx <- B_burnin + 1
    end_idx <- N - C_cutoff
    dataTVPCDsim <- dataTVPCDsim[, start_idx:end_idx]
  }
  
  # Results matrix
  partvpCD <- matrix(0, M, 7)
  
  # Estimation loop
  for (i in 1:M) {
    
    # Generate starting point using the new implementation's method
    y_for_start <- dataTVPCDsim[i,]  # Use this data series for starting point generation
    starting_points <- generate_starting_points(y_for_start, K = 2, model_type = "tvp", n_starts = 1, seed = i)
    par_start_new_format <- starting_points[[1]]
    
    # Convert from new implementation format to original format
    # New format: [mu1, mu2, sigma2_1, sigma2_2, init_trans, A]  
    # Original format: [mu1, mu2, sigma2, p11, p22, A1, A2]
    mu_start <- par_start_new_format[1:2]
    sigma2_start <- par_start_new_format[3]  # Use first variance (they should be equal anyway)
    init_trans_start <- par_start_new_format[5:6]  # off-diagonal transitions
    A_start <- par_start_new_format[7:8]
    
    # Convert transitions: new format -> original format
    p11_start <- 1 - init_trans_start[1]
    p22_start <- 1 - init_trans_start[2]
    
    par_start_original <- c(mu_start[1], mu_start[2], sigma2_start, p11_start, p22_start, A_start[1], A_start[2])
    
    TVP.est <- nlminb(
      start = par.trasf.inv_TVP(par_start_original), 
      filtering.single.trasf_TVP, 
      y = dataTVPCDsim[i,], 
      B = B_burnin, 
      control = list(eval.max = 1e+06, iter.max = 1e+06, trace = 0)
    )
    
    TVP.par.fin <- par.trasf_TVP(TVP.est$par)
    partvpCD[i,] <- TVP.par.fin
  }
  
  return(partvpCD)
}






######XExo


dataexoCD_original <- function(M, N, mu, sigma2, init_trans, A, X_Exo) {
  
  data=matrix(0,M,N)
  
  # Convert to original internal format:
  mu0 = mu[1]; mu1 = mu[2]
  p11 = 1 - init_trans[1]; p22 = 1 - init_trans[2]  
  A1 = A[1]; A2 = A[2]
  
  # Key difference from TVP: exogenous model uses different omega calculation
  omega1_LR = log(p11/(1-p11))
  omega2_LR = log(p22/(1-p22))
  
  # For exogenous model: omega = omega_LR * (1 - A), just like TVP
  omega1 = omega1_LR*(1-A1)
  omega2 = omega2_LR*(1-A2)
  
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


run_original_EXO_estimation <- function(M, N, mu, sigma2, init_trans, A, X_Exo, B_burnin = 0, C_cutoff = 0) {
  
  # Generate data using our new parameterized function (full series)
  dataexoCDsim <- dataexoCD_original(M, N, mu, sigma2, init_trans, A, X_Exo)
  
  # Apply burn-in and cutoff if needed
  if (B_burnin > 0 || C_cutoff > 0) {
    start_idx <- B_burnin + 1
    end_idx <- N - C_cutoff
    dataexoCDsim <- dataexoCDsim[, start_idx:end_idx]
    # Also trim X_Exo to match the trimmed data
    X_Exo_trimmed <- X_Exo[start_idx:end_idx]
  } else {
    X_Exo_trimmed <- X_Exo
  }
  
  # Results matrix
  parexoCD <- matrix(0, M, 7)
  
  # Estimation loop
  for (i in 1:M) {
    
    # Generate starting point using the new implementation's method
    y_for_start <- dataexoCDsim[i,]
    starting_points <- generate_starting_points(y_for_start, K = 2, model_type = "exogenous", n_starts = 1, seed = i)
    par_start_new_format <- starting_points[[1]]
    
    # Convert from new implementation format to original format
    # New format: [mu1, mu2, sigma2_1, sigma2_2, init_trans, A]  
    # Original format: [mu1, mu2, sigma2, p11, p22, A1, A2]
    mu_start <- par_start_new_format[1:2]
    sigma2_start <- par_start_new_format[3]  # Use first variance
    init_trans_start <- par_start_new_format[5:6]  # off-diagonal transitions
    A_start <- par_start_new_format[7:8]
    
    # Convert transitions: new format -> original format
    p11_start <- 1 - init_trans_start[1]
    p22_start <- 1 - init_trans_start[2]
    
    par_start_original <- c(mu_start[1], mu_start[2], sigma2_start, p11_start, p22_start, A_start[1], A_start[2])
    
    TVP.XExo.est <- nlminb(
      start = par.trasf.inv_TVPXExo(par_start_original), 
      filtering.single.trasf_TVPXExo, 
      y = dataexoCDsim[i,], 
      B = B_burnin, 
      X_Exo = X_Exo_trimmed,  # KEY: Pass the exogenous variable
      control = list(eval.max = 1e+06, iter.max = 1e+06, trace = 0)
    )
    
    TVP.XExo.par.fin <- par.trasf_TVPXExo(TVP.XExo.est$par)[1:7]  # Take first 7 elements like original
    parexoCD[i,] <- TVP.XExo.par.fin
  }
  
  return(parexoCD)
}


