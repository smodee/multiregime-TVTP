/*
 * C implementations of filtering functions for multiregimeTVTP
 * Uses plain R C API (.Call interface) - no Rcpp dependency
 *
 * Implements: Cfiltering_Const, Cfiltering_TVP, Cfiltering_Exo, Cfiltering_GAS
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define EPS_MACHINE 2.220446e-16
#define EPS_PROB    1e-10
#define LOG_EPS     -36.04365

/* =========================================================================
   Scalar helpers
   ========================================================================= */

static double logit_s(double x) {
  if (x <= EPS_MACHINE) x = EPS_MACHINE;
  if (x >= 1.0 - EPS_MACHINE) x = 1.0 - EPS_MACHINE;
  return log(x / (1.0 - x));
}

static double logistic_s(double x) {
  return 1.0 / (1.0 + exp(-x));
}

static double logistic_clamped_s(double x) {
  return EPS_PROB + (1.0 - 2.0 * EPS_PROB) / (1.0 + exp(-x));
}

/* =========================================================================
   Build K x K transition matrix (row-major flat array)
   ========================================================================= */

/* Diagonal parameterization: diag_probs[K] -> P[K*K] */
static void build_tmat_diag(const double *dp, int K, double *P) {
  for (int i = 0; i < K; i++) {
    double off = (1.0 - dp[i]) / (K - 1);
    for (int j = 0; j < K; j++)
      P[i*K + j] = (i == j) ? dp[i] : off;
  }
}

/* Off-diagonal parameterization: off_diag[K*(K-1)] -> P[K*K] */
static void build_tmat_offdiag(const double *od, int K, double *P) {
  int idx = 0;
  for (int i = 0; i < K; i++) {
    double rs = 0.0;
    for (int j = 0; j < K; j++) {
      if (i != j) { P[i*K + j] = od[idx]; rs += od[idx]; idx++; }
    }
    P[i*K + i] = 1.0 - rs;
  }
}

static void build_tmat(const double *tp, int n_trans, int diag_probs, int K, double *P) {
  if (diag_probs) build_tmat_diag(tp, K, P);
  else build_tmat_offdiag(tp, K, P);
}

/* =========================================================================
   Clamp off-diagonal row sums
   ========================================================================= */
static void clamp_offdiag(double *p, int K) {
  int npr = K - 1;
  for (int i = 0; i < K; i++) {
    double s = 0.0;
    for (int j = 0; j < npr; j++) s += p[i*npr + j];
    if (s > 1.0 - 1e-6) {
      double sc = (1.0 - 1e-6) / s;
      for (int j = 0; j < npr; j++) p[i*npr + j] *= sc;
    }
  }
}

/* =========================================================================
   mat^T * vec: result = P^T * x (P row-major K x K)
   ========================================================================= */
static void matT_vec(const double *P, int K, const double *x, double *result) {
  memset(result, 0, K * sizeof(double));
  for (int i = 0; i < K; i++)
    for (int j = 0; j < K; j++)
      result[j] += P[i*K + j] * x[i];
}

/* =========================================================================
   Stationary distribution via power iteration
   ========================================================================= */
static void stat_dist(const double *P, int K, double *pi) {
  for (int k = 0; k < K; k++) pi[k] = 1.0 / K;
  double *pi_new = (double *)R_alloc(K, sizeof(double));

  for (int iter = 0; iter < 1000; iter++) {
    matT_vec(P, K, pi, pi_new);
    double s = 0.0;
    for (int k = 0; k < K; k++) s += pi_new[k];
    if (s > EPS_MACHINE) for (int k = 0; k < K; k++) pi_new[k] /= s;
    double md = 0.0;
    for (int k = 0; k < K; k++) {
      double d = fabs(pi_new[k] - pi[k]);
      if (d > md) md = d;
    }
    memcpy(pi, pi_new, K * sizeof(double));
    if (md < 1e-10) break;
  }
  double s = 0.0;
  for (int k = 0; k < K; k++) { if (pi[k] < EPS_MACHINE) pi[k] = EPS_MACHINE; s += pi[k]; }
  for (int k = 0; k < K; k++) pi[k] /= s;
}

/* =========================================================================
   Compute predicted probs: X_pred = P^T * X_prev, clamp, renormalize
   ========================================================================= */
static void predict_probs(const double *P, int K, const double *X_prev, double *X_pred) {
  matT_vec(P, K, X_prev, X_pred);
  double s = 0.0;
  for (int k = 0; k < K; k++) { if (X_pred[k] < EPS_MACHINE) X_pred[k] = EPS_MACHINE; s += X_pred[k]; }
  for (int k = 0; k < K; k++) X_pred[k] /= s;
}

/* =========================================================================
   f -> trans_prob via logistic (plain or clamped) + optional row clamp
   ========================================================================= */
static void f_to_tp(const double *f, int n, int diag_probs, int clamped, int K, double *p) {
  for (int i = 0; i < n; i++)
    p[i] = clamped ? logistic_clamped_s(f[i]) : logistic_s(f[i]);
  if (!diag_probs) clamp_offdiag(p, K);
}

/* =========================================================================
   Neg log-likelihood
   ========================================================================= */
static double neg_loglik(const double *tot_lik, int M, int n_burnin, int n_cutoff) {
  double nll = 0.0;
  for (int t = n_burnin; t < M - n_cutoff; t++) {
    double ll = log(tot_lik[t]);
    if (!R_FINITE(ll)) ll = LOG_EPS;
    nll -= ll;
  }
  return nll;
}

/* =========================================================================
   Pre-compute emission likelihoods eta[k*M + t] = dnorm(y[t], mu[k], sd[k])
   ========================================================================= */
static double *compute_eta(const double *mu, const double *sigma2, const double *y, int K, int M) {
  double *eta = (double *)R_alloc(K * M, sizeof(double));
  for (int k = 0; k < K; k++) {
    double sd_k = sqrt(sigma2[k]);
    for (int t = 0; t < M; t++)
      eta[k*M + t] = dnorm(y[t], mu[k], sd_k, 0);
  }
  return eta;
}

/* =========================================================================
   CONSTANT MODEL
   ========================================================================= */
static double do_filtering_const(const double *mu, const double *sigma2,
                                  const double *tp, int n_trans, const double *y,
                                  int M, int n_burnin, int n_cutoff, int diag_probs, int K) {
  double *P = (double *)R_alloc(K*K, sizeof(double));
  build_tmat(tp, n_trans, diag_probs, K, P);

  double *eta = compute_eta(mu, sigma2, y, K, M);
  double *sd_pi = (double *)R_alloc(K, sizeof(double));
  stat_dist(P, K, sd_pi);

  double *tot_lik = (double *)R_alloc(M, sizeof(double));
  double *X_t = (double *)R_alloc(K, sizeof(double));
  double *X_tlag = (double *)R_alloc(K, sizeof(double));

  /* t = 0 */
  matT_vec(P, K, sd_pi, X_tlag);
  double tl = 0.0;
  for (int k = 0; k < K; k++) tl += eta[k*M] * X_tlag[k];
  if (tl <= 0.0 || !R_FINITE(tl)) {
    tl = EPS_MACHINE;
    for (int k = 0; k < K; k++) X_t[k] = 1.0/K;
  } else {
    for (int k = 0; k < K; k++) X_t[k] = (eta[k*M] * X_tlag[k]) / tl;
  }
  tot_lik[0] = tl;

  /* Main loop */
  for (int t = 1; t < M; t++) {
    matT_vec(P, K, X_t, X_tlag);
    tl = 0.0;
    for (int k = 0; k < K; k++) tl += eta[k*M + t] * X_tlag[k];
    if (tl <= 0.0 || !R_FINITE(tl)) {
      tl = EPS_MACHINE;
      for (int k = 0; k < K; k++) X_t[k] = 1.0/K;
    } else {
      for (int k = 0; k < K; k++) X_t[k] = (eta[k*M + t] * X_tlag[k]) / tl;
    }
    tot_lik[t] = tl;
  }

  return neg_loglik(tot_lik, M, n_burnin, n_cutoff);
}

SEXP C_filtering_Const(SEXP mu_, SEXP sigma2_, SEXP tp_, SEXP y_,
                        SEXP n_burnin_, SEXP n_cutoff_, SEXP diag_probs_) {
  int K = LENGTH(mu_);
  int M = LENGTH(y_);
  double result = do_filtering_const(REAL(mu_), REAL(sigma2_), REAL(tp_),
                                      LENGTH(tp_), REAL(y_), M,
                                      INTEGER(n_burnin_)[0], INTEGER(n_cutoff_)[0],
                                      LOGICAL(diag_probs_)[0], K);
  return ScalarReal(result);
}

/* =========================================================================
   TVP MODEL
   ========================================================================= */
SEXP C_filtering_TVP(SEXP mu_, SEXP sigma2_, SEXP init_trans_, SEXP A_,
                      SEXP y_, SEXP n_burnin_, SEXP n_cutoff_, SEXP diag_probs_) {
  int K = LENGTH(mu_);
  int M = LENGTH(y_);
  int n_trans = LENGTH(init_trans_);
  int dp = LOGICAL(diag_probs_)[0];
  const double *mu = REAL(mu_), *sigma2 = REAL(sigma2_);
  const double *it = REAL(init_trans_), *A = REAL(A_), *y = REAL(y_);
  int n_burnin = INTEGER(n_burnin_)[0], n_cutoff = INTEGER(n_cutoff_)[0];

  double *eta = compute_eta(mu, sigma2, y, K, M);
  double *omega_LR = (double *)R_alloc(n_trans, sizeof(double));
  double *omega = (double *)R_alloc(n_trans, sizeof(double));
  for (int i = 0; i < n_trans; i++) {
    omega_LR[i] = logit_s(it[i]);
    omega[i] = omega_LR[i] * (1.0 - A[i]);
  }

  double *tot_lik = (double *)R_alloc(M, sizeof(double));
  double *X_t = (double *)R_alloc(K, sizeof(double));
  double *X_tlag = (double *)R_alloc(K, sizeof(double));
  double *f = (double *)R_alloc(n_trans, sizeof(double));
  double *pt = (double *)R_alloc(n_trans, sizeof(double));
  double *P = (double *)R_alloc(K*K, sizeof(double));
  double *sd_pi = (double *)R_alloc(K, sizeof(double));

  /* t=0: f = omega_LR */
  memcpy(f, omega_LR, n_trans * sizeof(double));
  f_to_tp(f, n_trans, dp, 0, K, pt);
  build_tmat(pt, n_trans, dp, K, P);
  stat_dist(P, K, sd_pi);
  matT_vec(P, K, sd_pi, X_tlag);

  double tl = 0.0;
  for (int k = 0; k < K; k++) tl += eta[k*M] * X_tlag[k];
  if (tl <= 0.0 || !R_FINITE(tl)) { tl = EPS_MACHINE; for (int k=0;k<K;k++) X_t[k]=1.0/K; }
  else { for (int k=0;k<K;k++) X_t[k]=(eta[k*M]*X_tlag[k])/tl; }
  tot_lik[0] = tl;

  /* CRITICAL: f[:,1] = omega */
  if (M >= 2) { memcpy(f, omega, n_trans*sizeof(double)); f_to_tp(f,n_trans,dp,0,K,pt); }

  /* Main loop: t=1..M-2 */
  for (int t = 1; t < M-1; t++) {
    build_tmat(pt, n_trans, dp, K, P);
    predict_probs(P, K, X_t, X_tlag);
    tl = 0.0;
    for (int k=0;k<K;k++) tl += eta[k*M+t]*X_tlag[k];
    if (tl<=0.0||!R_FINITE(tl)) { tl=EPS_MACHINE; for(int k=0;k<K;k++) X_t[k]=1.0/K; }
    else { for(int k=0;k<K;k++) X_t[k]=(eta[k*M+t]*X_tlag[k])/tl; }
    tot_lik[t] = tl;
    if (t < M-2) {
      for (int i=0;i<n_trans;i++) f[i] = omega[i] + A[i]*y[t];
      f_to_tp(f,n_trans,dp,0,K,pt);
    }
  }

  /* Last time point */
  if (M >= 2) {
    if (M > 2) { for(int i=0;i<n_trans;i++) f[i]=omega[i]+A[i]*y[M-2]; f_to_tp(f,n_trans,dp,0,K,pt); }
    build_tmat(pt,n_trans,dp,K,P);
    predict_probs(P,K,X_t,X_tlag);
    tl=0.0; for(int k=0;k<K;k++) tl+=eta[k*M+M-1]*X_tlag[k];
    if(tl<=0.0||!R_FINITE(tl)) tl=EPS_MACHINE;
    tot_lik[M-1]=tl;
  }

  return ScalarReal(neg_loglik(tot_lik, M, n_burnin, n_cutoff));
}

/* =========================================================================
   EXOGENOUS MODEL
   ========================================================================= */
SEXP C_filtering_Exo(SEXP mu_, SEXP sigma2_, SEXP init_trans_, SEXP A_,
                      SEXP y_, SEXP X_Exo_, SEXP n_burnin_, SEXP n_cutoff_, SEXP diag_probs_) {
  int K = LENGTH(mu_);
  int M = LENGTH(y_);
  int n_trans = LENGTH(init_trans_);
  int dp = LOGICAL(diag_probs_)[0];
  const double *mu=REAL(mu_), *sigma2=REAL(sigma2_);
  const double *it=REAL(init_trans_), *A=REAL(A_), *y=REAL(y_), *xexo=REAL(X_Exo_);
  int n_burnin=INTEGER(n_burnin_)[0], n_cutoff=INTEGER(n_cutoff_)[0];

  double *eta = compute_eta(mu,sigma2,y,K,M);
  double *omega_LR=(double*)R_alloc(n_trans,sizeof(double));
  double *omega=(double*)R_alloc(n_trans,sizeof(double));
  for(int i=0;i<n_trans;i++) { omega_LR[i]=logit_s(it[i]); omega[i]=omega_LR[i]*(1.0-A[i]); }

  double *tot_lik=(double*)R_alloc(M,sizeof(double));
  double *X_t=(double*)R_alloc(K,sizeof(double));
  double *X_tlag=(double*)R_alloc(K,sizeof(double));
  double *f=(double*)R_alloc(n_trans,sizeof(double));
  double *pt=(double*)R_alloc(n_trans,sizeof(double));
  double *P=(double*)R_alloc(K*K,sizeof(double));
  double *sd_pi=(double*)R_alloc(K,sizeof(double));

  /* t=0: f = omega + A*X_Exo[0] */
  for(int i=0;i<n_trans;i++) f[i]=omega[i]+A[i]*xexo[0];
  f_to_tp(f,n_trans,dp,0,K,pt);
  build_tmat(pt,n_trans,dp,K,P); stat_dist(P,K,sd_pi); matT_vec(P,K,sd_pi,X_tlag);
  double tl=0.0; for(int k=0;k<K;k++) tl+=eta[k*M]*X_tlag[k];
  if(tl<=0.0||!R_FINITE(tl)){tl=EPS_MACHINE;for(int k=0;k<K;k++)X_t[k]=1.0/K;}
  else{for(int k=0;k<K;k++)X_t[k]=(eta[k*M]*X_tlag[k])/tl;}
  tot_lik[0]=tl;

  if(M>=2){for(int i=0;i<n_trans;i++)f[i]=omega[i]+A[i]*xexo[1];f_to_tp(f,n_trans,dp,0,K,pt);}

  for(int t=1;t<M-1;t++){
    build_tmat(pt,n_trans,dp,K,P); predict_probs(P,K,X_t,X_tlag);
    tl=0.0; for(int k=0;k<K;k++)tl+=eta[k*M+t]*X_tlag[k];
    if(tl<=0.0||!R_FINITE(tl)){tl=EPS_MACHINE;for(int k=0;k<K;k++)X_t[k]=1.0/K;}
    else{for(int k=0;k<K;k++)X_t[k]=(eta[k*M+t]*X_tlag[k])/tl;}
    tot_lik[t]=tl;
    if(t<M-2){for(int i=0;i<n_trans;i++)f[i]=omega[i]+A[i]*xexo[t+1];f_to_tp(f,n_trans,dp,0,K,pt);}
  }

  if(M>=2){
    if(M>2){for(int i=0;i<n_trans;i++)f[i]=omega[i]+A[i]*xexo[M-1];f_to_tp(f,n_trans,dp,0,K,pt);}
    build_tmat(pt,n_trans,dp,K,P); predict_probs(P,K,X_t,X_tlag);
    tl=0.0; for(int k=0;k<K;k++)tl+=eta[k*M+M-1]*X_tlag[k];
    if(tl<=0.0||!R_FINITE(tl))tl=EPS_MACHINE;
    tot_lik[M-1]=tl;
  }

  return ScalarReal(neg_loglik(tot_lik,M,n_burnin,n_cutoff));
}

/* =========================================================================
   GAS MODEL HELPERS
   ========================================================================= */

/* Fisher info - diagonal parameterization */
static double fisher_diag(const double *mu, const double *sigma2, const double *X_prev,
                           const double *p_diag, const double *nodes, const double *weights,
                           int n_nodes, double mu_q, double sigma_q, int K) {
  double fi = 0.0;
  double *P = (double *)R_alloc(K*K, sizeof(double));
  double *X_pred = (double *)R_alloc(K, sizeof(double));
  build_tmat_diag(p_diag, K, P);
  matT_vec(P, K, X_prev, X_pred);

  for (int j = 0; j < n_nodes; j++) {
    double yn = nodes[j];
    double d0 = dnorm(yn, mu[0], sqrt(sigma2[0]), 0);
    double dK = dnorm(yn, mu[K-1], sqrt(sigma2[K-1]), 0);
    double td = 0.0;
    for (int k = 0; k < K; k++) {
      double dk = (k==0)?d0:(k==K-1)?dK:dnorm(yn,mu[k],sqrt(sigma2[k]),0);
      td += dk * X_pred[k];
    }
    if (td > EPS_MACHINE) {
      double ld = d0 - dK;
      double wc = sqrt(2.0*M_PI*sigma_q*sigma_q)*exp((yn-mu_q)*(yn-mu_q)/(2.0*sigma_q*sigma_q));
      fi += (ld*ld)/td * wc * weights[j];
    }
  }
  if (!R_FINITE(fi) || fi <= 0.0) fi = 1.0;
  if (fi < EPS_MACHINE) fi = EPS_MACHINE;
  return fi;
}

/* Fisher info - off-diagonal parameterization */
static double fisher_offdiag(const double *mu, const double *sigma2, const double *X_prev,
                              const double *p_trans, int n_trans, const double *nodes,
                              const double *weights, int n_nodes, double mu_q, double sigma_q, int K) {
  double fi = 0.0;
  /* Clamp probs */
  double *pv = (double *)R_alloc(n_trans, sizeof(double));
  int npr = K-1;
  memcpy(pv, p_trans, n_trans*sizeof(double));
  for (int i=0;i<K;i++){
    double rs=0;
    for(int j=0;j<npr;j++){int idx=i*npr+j;if(pv[idx]<0.001)pv[idx]=0.001;if(pv[idx]>0.99)pv[idx]=0.99;rs+=pv[idx];}
    if(rs>0.99){double sc=0.99/rs;for(int j=0;j<npr;j++)pv[i*npr+j]*=sc;}
  }
  double *P=(double*)R_alloc(K*K,sizeof(double));
  double *X_pred=(double*)R_alloc(K,sizeof(double));
  build_tmat_offdiag(pv,K,P); matT_vec(P,K,X_prev,X_pred);

  for(int j=0;j<n_nodes;j++){
    double yn=nodes[j];
    double *rd=(double*)R_alloc(K,sizeof(double));
    for(int k=0;k<K;k++) rd[k]=dnorm(yn,mu[k],sqrt(sigma2[k]),0);
    double td=0; for(int k=0;k<K;k++) td+=rd[k]*X_pred[k];
    if(td>EPS_MACHINE){
      double is=0;
      for(int i=0;i<K;i++)for(int l=0;l<K;l++)if(i!=l){double ld=rd[i]-rd[l];is+=X_pred[i]*(ld*ld)/td;}
      double wc=sqrt(2.0*M_PI*sigma_q*sigma_q)*exp((yn-mu_q)*(yn-mu_q)/(2.0*sigma_q*sigma_q));
      fi+=is*wc*weights[j];
    }
  }
  if(!R_FINITE(fi)||fi<0.0)fi=1.0;
  if(fi<EPS_MACHINE)fi=EPS_MACHINE;
  return fi;
}

/* GAS score for one time step */
static void gas_score(double y_obs, const double *mu, const double *sigma2,
                       const double *tp, int n_trans, int diag_probs,
                       const double *X_pred, double tl_ext,
                       const double *nodes, const double *weights, int n_nodes,
                       double mu_q, double sigma_q, int K, double *score_out) {
  memset(score_out, 0, n_trans*sizeof(double));
  if (tl_ext <= EPS_MACHINE || !R_FINITE(tl_ext)) return;

  double *eta_s = (double *)R_alloc(K, sizeof(double));
  for (int k=0;k<K;k++) eta_s[k] = dnorm(y_obs, mu[k], sqrt(sigma2[k]), 0);

  if (diag_probs) {
    /* Factored scaling */
    double fi = fisher_diag(mu,sigma2,X_pred,tp,nodes,weights,n_nodes,mu_q,sigma_q,K);
    double S = (eta_s[0]-eta_s[K-1])/tl_ext;
    double *gv=(double*)R_alloc(K,sizeof(double));
    int sign=1;
    for(int i=0;i<K;i++){gv[i]=sign*X_pred[i]*tp[i]*(1.0-tp[i]);sign=-sign;}
    double gm=0; for(int i=0;i<K;i++) gm+=gv[i]*gv[i]; gm=sqrt(gm);
    double Ss = (fi>EPS_MACHINE)?S/sqrt(fi):0.0;
    if(gm>EPS_MACHINE) for(int i=0;i<K;i++) score_out[i]=Ss*gv[i]/gm;
    for(int i=0;i<K;i++) if(!R_FINITE(score_out[i])){memset(score_out,0,n_trans*sizeof(double));return;}
  } else {
    /* Simple scaling */
    double fi = fisher_offdiag(mu,sigma2,X_pred,tp,n_trans,nodes,weights,n_nodes,mu_q,sigma_q,K);
    double fis=(fi>EPS_MACHINE)?1.0/sqrt(fi):0.0;
    int idx=0;
    for(int i=0;i<K;i++)for(int j=0;j<K;j++)if(i!=j){
      double ld=(eta_s[i]-eta_s[j])/tl_ext;
      double gc=X_pred[i]*tp[idx]*(1.0-tp[idx]);
      score_out[idx]=ld*gc*fis;
      idx++;
    }
    for(int i=0;i<n_trans;i++) if(!R_FINITE(score_out[i])){memset(score_out,0,n_trans*sizeof(double));return;}
  }
}

/* =========================================================================
   GAS MODEL
   ========================================================================= */
SEXP C_filtering_GAS(SEXP mu_, SEXP sigma2_, SEXP init_trans_, SEXP A_, SEXP B_,
                      SEXP y_, SEXP n_burnin_, SEXP n_cutoff_, SEXP diag_probs_,
                      SEXP gh_nodes_, SEXP gh_weights_, SEXP mu_quad_, SEXP sigma_quad_,
                      SEXP use_fallback_, SEXP A_threshold_) {
  int K=LENGTH(mu_), M=LENGTH(y_), n_trans=LENGTH(init_trans_);
  int dp=LOGICAL(diag_probs_)[0];
  const double *mu=REAL(mu_),*sigma2=REAL(sigma2_),*it=REAL(init_trans_);
  const double *A=REAL(A_),*B=REAL(B_),*y=REAL(y_);
  int n_burnin=INTEGER(n_burnin_)[0],n_cutoff=INTEGER(n_cutoff_)[0];
  const double *ghn=REAL(gh_nodes_),*ghw=REAL(gh_weights_);
  double mu_q=REAL(mu_quad_)[0],sigma_q=REAL(sigma_quad_)[0];
  int n_nodes=LENGTH(gh_nodes_);
  int uf=LOGICAL(use_fallback_)[0];
  double at=REAL(A_threshold_)[0];

  /* Fallback */
  if(uf){
    double mA=0; for(int i=0;i<n_trans;i++){double a=fabs(A[i]);if(a>mA)mA=a;}
    if(mA<at) return ScalarReal(do_filtering_const(mu,sigma2,it,n_trans,y,M,n_burnin,n_cutoff,dp,K));
  }

  double *eta=compute_eta(mu,sigma2,y,K,M);
  double *omega=(double*)R_alloc(n_trans,sizeof(double));
  for(int i=0;i<n_trans;i++) omega[i]=logit_s(it[i]);

  double *tot_lik=(double*)R_alloc(M,sizeof(double));
  double *X_t=(double*)R_alloc(K,sizeof(double));
  double *X_t_prev=(double*)R_alloc(K,sizeof(double));
  double *X_tlag=(double*)R_alloc(K,sizeof(double));
  double *f=(double*)R_alloc(n_trans,sizeof(double));
  double *f_prev=(double*)R_alloc(n_trans,sizeof(double));
  double *pt=(double*)R_alloc(n_trans,sizeof(double));
  double *P=(double*)R_alloc(K*K,sizeof(double));
  double *sd_pi=(double*)R_alloc(K,sizeof(double));
  double *sc=(double*)R_alloc(n_trans,sizeof(double));

  /* t=0 */
  memcpy(f,omega,n_trans*sizeof(double));
  f_to_tp(f,n_trans,dp,0,K,pt);
  build_tmat(pt,n_trans,dp,K,P); stat_dist(P,K,sd_pi);
  matT_vec(P,K,sd_pi,X_tlag);
  double tl=0; for(int k=0;k<K;k++) tl+=eta[k*M]*X_tlag[k];
  if(tl<=0.0||!R_FINITE(tl)){tl=EPS_MACHINE;for(int k=0;k<K;k++)X_t[k]=1.0/K;}
  else{for(int k=0;k<K;k++)X_t[k]=(eta[k*M]*X_tlag[k])/tl;
    double s=0;for(int k=0;k<K;k++){if(X_t[k]<EPS_MACHINE)X_t[k]=EPS_MACHINE;s+=X_t[k];}
    for(int k=0;k<K;k++)X_t[k]/=s;}
  tot_lik[0]=tl;

  /* Score at t=0, set f for t=1 */
  gas_score(y[0],mu,sigma2,pt,n_trans,dp,sd_pi,tot_lik[0],ghn,ghw,n_nodes,mu_q,sigma_q,K,sc);
  memcpy(f_prev,f,n_trans*sizeof(double));
  if(M>=2){
    for(int i=0;i<n_trans;i++) f[i]=omega[i]+A[i]*sc[i]+B[i]*(f_prev[i]-omega[i]);
    f_to_tp(f,n_trans,dp,1,K,pt);
  }

  /* Main loop */
  for(int t=1;t<M;t++){
    memcpy(X_t_prev,X_t,K*sizeof(double));
    memcpy(f_prev,f,n_trans*sizeof(double));

    build_tmat(pt,n_trans,dp,K,P); predict_probs(P,K,X_t_prev,X_tlag);
    tl=0; for(int k=0;k<K;k++) tl+=eta[k*M+t]*X_tlag[k];
    if(tl<=0.0||!R_FINITE(tl)){tl=EPS_MACHINE;for(int k=0;k<K;k++)X_t[k]=1.0/K;}
    else{for(int k=0;k<K;k++)X_t[k]=(eta[k*M+t]*X_tlag[k])/tl;
      double s=0;for(int k=0;k<K;k++){if(X_t[k]<EPS_MACHINE)X_t[k]=EPS_MACHINE;s+=X_t[k];}
      for(int k=0;k<K;k++)X_t[k]/=s;}
    tot_lik[t]=tl;

    if(t<M-1){
      gas_score(y[t],mu,sigma2,pt,n_trans,dp,X_t_prev,tot_lik[t],ghn,ghw,n_nodes,mu_q,sigma_q,K,sc);
      for(int i=0;i<n_trans;i++) f[i]=omega[i]+A[i]*sc[i]+B[i]*(f_prev[i]-omega[i]);
      f_to_tp(f,n_trans,dp,1,K,pt);
    }
  }
  if(tot_lik[M-1]<=0.0||!R_FINITE(tot_lik[M-1]))tot_lik[M-1]=EPS_MACHINE;

  return ScalarReal(neg_loglik(tot_lik,M,n_burnin,n_cutoff));
}

/* Debug test function */
SEXP C_test_call(SEXP x_) {
  double result = 0.0;
  int n = LENGTH(x_);
  const double *x = REAL(x_);
  for (int i = 0; i < n; i++) result += x[i];
  return ScalarReal(result);
}

/* =========================================================================
   Registration
   ========================================================================= */
static const R_CallMethodDef CallEntries[] = {
  {"C_test_call",       (DL_FUNC) &C_test_call,       1},
  {"C_filtering_Const", (DL_FUNC) &C_filtering_Const, 7},
  {"C_filtering_TVP",   (DL_FUNC) &C_filtering_TVP,   8},
  {"C_filtering_Exo",   (DL_FUNC) &C_filtering_Exo,   9},
  {"C_filtering_GAS",   (DL_FUNC) &C_filtering_GAS,   15},
  {NULL, NULL, 0}
};

void R_init_multiregimeTVTP(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
