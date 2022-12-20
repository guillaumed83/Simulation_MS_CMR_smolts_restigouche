
## running M1 - M4

setwd("C:/Users/DauphinGU/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Simulation")

library(gtools)
library(coda)
#library(rjags)
library(R2jags) #more flexible to calculate WAIC and deviance
#load.module("glm")
library(dclone) #for multicore computing
library(boot)
library(tidyverse)
library(MCMCvis)
library(mcmcplots)


rm(list=ls())

source("Final_data_simulation_replicates_LogN.R")
#source("C:/Users/DauphinGU/Desktop/smolts_restigouche/Model_CMR/2019-06-27_whole_season_data2019/_data_simulation.R")


set.seed(20221214)

n_burnin <- 100000
n_iter <- 3100000
n_thin <- 600



## n.sims = (n.iter - n.burnin) / n.thin
(n_iter - n_burnin)/n_thin # *n_chain

monitor_=c('mu_theta',
           'theta',
           'sigma_theta',
           'Nm','Nm_tot','Nm_rest','p_smolt_prod',
           'mu_lambda','beta_lambda','alpha_lambda',
           'delta','eta_alphaN'
)

## M5
jags_DM_simul<-list()
samples_DM_simul<-list()
waic_DM <- list()

for(i in 1:n_rep){
  jags_DM_simul[[i]] <- jags(model.file="C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Models/M5_simul_concentration_factor_unif_Nmtot_LogN.bug",
                                     data=data_pooled[[i]],
                                     inits=inits,
                                     parameters.to.save = monitor_,
                                     n.chains=2,
                                     n.burnin = n_burnin,
                                     n.iter = n_iter,
                                     n.thin= n_thin
  )
  
  #mcmcplot(jags_indpt_no_hier_simul )
  jags_DM_simul[[i]]
  
  samples_DM_simul[[i]] <-jags.samples(jags_DM_simul[[i]]$model, c("WAIC", "deviance"),type="mean",
                                               n.burnin = n_burnin,
                                               n.iter = n_iter,
                                               n.thin=n_thin)
  
  samples_DM_simul[[i]]$p_waic <- samples_DM_simul[[i]]$WAIC
  samples_DM_simul[[i]]$waic <- samples_DM_simul[[i]]$deviance + samples_DM_simul[[i]]$p_waic
  tmp <- sapply(samples_DM_simul[[i]], sum)
  waic_DM[[i]] <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
  waic_DM[[i]]
  
  
  
  
}

var_names_parameters_DM_simul <- varnames(jags_DM_simul)


save.image("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_final_update_2022-12-14_M5_LogN.RData")





