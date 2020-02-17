rm(list=ls())
library(penalized)
library(optparse)
library(dplyr)
library(flashr)

setwd('C:/Users/aclou/Grad School/2020 Autumn/Thesis/ts_eQTLs-master/ts_eQTLs-master/')
source("sn_spMF/readIn.R")
source("sn_spMF/Update_FL.R")


run_nnf = function(X_mtx, W_mtx, num_facs, max_iters = 100, penalty_L = 0.1, penalty_F = 0.1){
  # Write to table first in order to comply with "readIn" function format
  write.table(X_mtx, 'C:/Users/aclou/Grad School/2020 Autumn/Thesis/ts_eQTLs-master/ts_eQTLs-master/X_sim.txt', row.names = FALSE)
  write.table(W_mtx, 'C:/Users/aclou/Grad School/2020 Autumn/Thesis/ts_eQTLs-master/ts_eQTLs-master/W_sim.txt', row.names = FALSE)

  # Define hyperparams for test
  K = num_facs
  alpha1 = penalty_L
  lambda1 = penalty_F
  inputdir='C:/Users/aclou/Grad School/2020 Autumn/Thesis/ts_eQTLs-master/ts_eQTLs-master/'
  outputdir = inputdir
  xfn=paste0(inputdir,'X_sim.txt')
  wfn=paste0(inputdir,'W_sim.txt')
  Data = readIn(K, alpha1, lambda1, xfn, wfn) 
  X = Data[['X']];
  W = Data[['W']];
  option = Data[['option']];
  option[['iter']]  = max_iters;
  Fn_basename = Data[['Fn_basename']];
  
  ## run MF to learn the factors  
  print(paste0('K=', (K), '; alpha1=', (alpha1),'; lambda1=', (lambda1)));
  Run_iter = Update_FL(X, W, option);
  FactorM = Run_iter[[1]]
  LoadingM = Run_iter[[2]]
  factor_corr = norm(cor(FactorM), 'F')
  L_sparsity = Run_iter[[3]]
  F_sparsity = Run_iter[[4]]
  print(paste0('Sparsity in Factor matrix =', (F_sparsity),'; Sparsity in L =', (L_sparsity), '; '))
  print(paste0((Run_iter[[5]]), ' factors remain; ; correlation between factors = ', (factor_corr)));
  return(list(FactorM = FactorM, LoadingM = LoadingM))
}

# Simulate different matrices and compare empirical results between methods -- Frobenius norm
results_df = matrix(0, nrow = 0, ncol = 13) %>% as.data.frame()

set.seed(9999)
for(n_fac in c(2, 4, 10)){
  for (n_col in c(10, 25, 50)){
    for (n_obs in c(100, 500, 1000, 5000)){
      for (sparsity_lvl in c(0.25, 0.5, 0.9)){
        seed_num = sample(seq(1, 1000000, 1), 1)
        set.seed(seed_num)
        sparseF_sim = rbinom(n = n_fac*n_col, 1, prob = sparsity_lvl)
        F_sim = ifelse(sparseF_sim == 1, 0, rnorm(n = n_fac*n_col, 2, 1)) %>% matrix(nrow = n_col) # replicate F+ (non-neg)  
        # dim(F_sim)
        sparseL_sim = rbinom(n = n_fac*n_obs, 1, prob = sparsity_lvl)
        L_sim = ifelse(sparseL_sim == 1, 0, rnorm(n = n_fac*n_obs, 0, 1)) %>% matrix(nrow = n_obs)
        # dim(L_sim)
        E_sim = rnorm(n = n_obs*n_col, 0, 0.5) %>% matrix(nrow = n_obs)
        X_sim = L_sim%*%t(F_sim) + E_sim
        W_sim = matrix(1, nrow = n_obs, ncol = n_col)
        penal_L = 0.5 # vary this too
        penal_F = 0.5 # vary this too
        t_start = Sys.time()
        nnf_out = run_nnf(X_mtx = X_sim, W_mtx = W_sim, num_facs = n_fac, max_iters = 50, penalty_L = penal_L, penalty_F = penal_F)
        t_nnf = as.numeric(Sys.time() - t_start)
        norm_nnf = sum((nnf_out$FactorM - F_sim)^2)
        norm_nnf_X = sum((nnf_out$LoadingM %*% t(nnf_out$FactorM) - X_sim)^2)
        t_start = Sys.time()
        flash_X = flash(X_sim, backfit = TRUE, greedy = TRUE)
        t_ebmf = Sys.time() - t_start
        ldf_X = flash_X$ldf
        norm_ebmf = ifelse(ncol(ldf_X$f) == ncol(F_sim), 
                           sum((ldf_X$f - F_sim)^2), 
                           -1) # -1 if num predicted factors don't match
        norm_ebmf_X = ifelse(length(ldf_X$d) > 1,
                             sum((ldf_X$l %*% diag(ldf_X$d) %*% t(ldf_X$f) - X_sim)^2),
                             sum((ldf_X$l %*% ldf_X$d %*% t(ldf_X$f) - X_sim)^2)) # if 1 factor predicted, no 'diag'
        results_df = rbind(results_df,
                           c(seed_num, n_fac, n_col, n_obs, 
                             sparsity_lvl, penal_L, penal_F, 
                             norm_nnf, norm_ebmf, norm_nnf_X, norm_ebmf_X, t_nnf, t_ebmf))
      }
    }
  }
}

colnames(results_df) = c("seed", "n_facs", "n_col", "n_obs", 
                         "sparsity", "penal_L", "penal_F", 
                         "norm_nnf", "norm_ebmf", "norm_nnf_X", "norm_ebmf_X", "t_nnf", "t_ebmf")
results_df






