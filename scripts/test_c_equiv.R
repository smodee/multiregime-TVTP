#!/usr/bin/env Rscript
# Test numerical equivalence between R and C filtering implementations
library(multiregimeTVTP)
cat("cpp_available:", cpp_available(), "\n\n")

# === CONSTANT MODEL K=2 ===
set.seed(42)
par_c <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90), K=2, model_type="constant",
                                   diag_probs=TRUE, equal_variances=TRUE)
data_c <- dataConstCD(200, 1, par_c)
nll_r <- Rfiltering_Const(par_c, data_c[,1], 5, 0, use_cpp=FALSE)
nll_c <- Rfiltering_Const(par_c, data_c[,1], 5, 0, use_cpp=TRUE)
cat(sprintf("Const K=2 diag eq_var: R=%.6f C=%.6f diff=%.2e %s\n",
            nll_r, nll_c, abs(nll_r-nll_c), ifelse(abs(nll_r-nll_c)<1e-8, "PASS", "FAIL")))

# Constant K=2, separate variances
par_c2 <- set_parameter_attributes(c(-1, 1, 0.8, 1.2, 0.95, 0.90), K=2, model_type="constant",
                                    diag_probs=TRUE, equal_variances=FALSE)
nll_r2 <- Rfiltering_Const(par_c2, data_c[,1], 5, 0, use_cpp=FALSE)
nll_c2 <- Rfiltering_Const(par_c2, data_c[,1], 5, 0, use_cpp=TRUE)
cat(sprintf("Const K=2 diag sep_var: R=%.6f C=%.6f diff=%.2e %s\n",
            nll_r2, nll_c2, abs(nll_r2-nll_c2), ifelse(abs(nll_r2-nll_c2)<1e-8, "PASS", "FAIL")))

# Constant K=3
set.seed(42)
par_c3 <- set_parameter_attributes(c(-2, 0, 2, 1, 1, 1, 0.90, 0.85, 0.80), K=3, model_type="constant",
                                    diag_probs=TRUE, equal_variances=FALSE)
data_c3 <- dataConstCD(200, 1, par_c3)
nll_r3 <- Rfiltering_Const(par_c3, data_c3[,1], 5, 0, use_cpp=FALSE)
nll_c3 <- Rfiltering_Const(par_c3, data_c3[,1], 5, 0, use_cpp=TRUE)
cat(sprintf("Const K=3 diag:        R=%.6f C=%.6f diff=%.2e %s\n",
            nll_r3, nll_c3, abs(nll_r3-nll_c3), ifelse(abs(nll_r3-nll_c3)<1e-8, "PASS", "FAIL")))

# === TVP MODEL K=2 ===
set.seed(42)
par_t <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90, 0.05, 0.05), K=2, model_type="tvp",
                                   diag_probs=TRUE, equal_variances=TRUE)
data_t <- dataTVPCD(200, 1, par_t)
nll_rt <- Rfiltering_TVP(par_t, data_t[,1], 5, 0, use_cpp=FALSE)
nll_ct <- Rfiltering_TVP(par_t, data_t[,1], 5, 0, use_cpp=TRUE)
cat(sprintf("TVP K=2 diag:          R=%.6f C=%.6f diff=%.2e %s\n",
            nll_rt, nll_ct, abs(nll_rt-nll_ct), ifelse(abs(nll_rt-nll_ct)<1e-8, "PASS", "FAIL")))

# === EXOGENOUS MODEL K=2 ===
set.seed(42)
par_e <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90, 0.05, 0.05), K=2, model_type="exogenous",
                                   diag_probs=TRUE, equal_variances=TRUE)
X_Exo <- rnorm(200)
data_e <- dataTVPXExoCD(200, 1, par_e, X_Exo=X_Exo)
nll_re <- Rfiltering_TVPXExo(par_e, X_Exo, data_e[,1], 5, 0, use_cpp=FALSE)
nll_ce <- Rfiltering_TVPXExo(par_e, X_Exo, data_e[,1], 5, 0, use_cpp=TRUE)
cat(sprintf("Exo K=2 diag:          R=%.6f C=%.6f diff=%.2e %s\n",
            nll_re, nll_ce, abs(nll_re-nll_ce), ifelse(abs(nll_re-nll_ce)<1e-8, "PASS", "FAIL")))

# === GAS MODEL K=2 ===
set.seed(42)
par_g <- set_parameter_attributes(c(-1, 1, 1, 0.95, 0.90, 0.3, 0.3, 0.5, 0.5), K=2, model_type="gas",
                                   diag_probs=TRUE, equal_variances=TRUE)
data_g <- dataGASCD(200, 1, par_g)
# Same seed for GH quadrature sampling
set.seed(123)
nll_rg <- Rfiltering_GAS(par_g, data_g[,1], 5, 0, use_cpp=FALSE)
set.seed(123)
nll_cg <- Rfiltering_GAS(par_g, data_g[,1], 5, 0, use_cpp=TRUE)
cat(sprintf("GAS K=2 diag:          R=%.6f C=%.6f diff=%.2e %s\n",
            nll_rg, nll_cg, abs(nll_rg-nll_cg), ifelse(abs(nll_rg-nll_cg)<1e-6, "PASS", "FAIL")))

cat("\nAll tests complete!\n")
