
## running M3 - M4

setwd("C:/Users/DauphinGU/Desktop/smolts_restigouche/Model_CMR/2023-07-30_jfb_model/Simulation")

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


set.seed(2021062)

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



## M3
##### M3_simul_indpt_hier.bug #####

jags_indpt_hier_simul<-list()
samples_indpt_hier_simul<-list()
waic_indpt_hier <- list()

for(i in 1:n_rep){
  jags_indpt_hier_simul[[i]] <- jags(model.file="C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2023-07-30_jfb_model/Models/M3_simul_indpt_hier_LogN.bug",
                                     data=data_pooled[[i]],
                                     inits=inits,
                                     parameters.to.save = monitor_,
                                     n.chains=2,
                                     n.burnin = n_burnin,
                                     n.iter = n_iter,
                                     n.thin= n_thin
  )
  
  #mcmcplot(jags_indpt_no_hier_simul )
  jags_indpt_hier_simul[[i]]
  
  samples_indpt_hier_simul[[i]] <-jags.samples(jags_indpt_hier_simul[[i]]$model, c("WAIC", "deviance"),type="mean",
                                               n.burnin = n_burnin,
                                               n.iter = n_iter,
                                               n.thin=n_thin)
  
  samples_indpt_hier_simul[[i]]$p_waic <- samples_indpt_hier_simul[[i]]$WAIC
  samples_indpt_hier_simul[[i]]$waic <- samples_indpt_hier_simul[[i]]$deviance + samples_indpt_hier_simul[[i]]$p_waic
  tmp <- sapply(samples_indpt_hier_simul[[i]], sum)
  waic_indpt_hier[[i]] <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
  waic_indpt_hier[[i]]
  
  
  
  
}

var_names_parameters_indpt_hier_simul <- varnames(jags_indpt_hier_simul)


## M4
##### M4_simul_indpt_hier_split.bug #####


monitor_=c('mu_theta',
           'theta',
           'sigma_theta',
           'Nm','Nm_tot','Nm_rest','p_smolt_prod',
           'mu_lambda','beta_lambda','alpha_lambda',
           'delta','eta_alphaN',
           'mu_split','tau_split', 'sigma_split','split'
)


jags_indpt_hier_split_simul<-list()
samples_indpt_hier_split_simul<-list()
waic_indpt_hier_split <- list()

for(i in 1:n_rep){
  jags_indpt_hier_split_simul[[i]] <- jags(model.file="C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2023-07-30_jfb_model/Models/M4_simul_indpt_hier_split_LogN.bug",
                                           data=data_pooled[[i]],
                                           inits=inits,
                                           parameters.to.save = monitor_,
                                           n.chains=2,
                                           n.burnin = n_burnin,
                                           n.iter = n_iter,
                                           n.thin= n_thin
  )
  
  #mcmcplot(jags_indpt_no_hier_simul )
  jags_indpt_hier_split_simul[[i]]
  
  samples_indpt_hier_split_simul[[i]] <-jags.samples(jags_indpt_hier_split_simul[[i]]$model, c("WAIC", "deviance"),type="mean",
                                                     n.burnin = n_burnin,
                                                     n.iter = n_iter,
                                                     n.thin=n_thin)
  
  samples_indpt_hier_split_simul[[i]]$p_waic <- samples_indpt_hier_split_simul[[i]]$WAIC
  samples_indpt_hier_split_simul[[i]]$waic <- samples_indpt_hier_split_simul[[i]]$deviance + samples_indpt_hier_split_simul[[i]]$p_waic
  tmp <- sapply(samples_indpt_hier_split_simul[[i]], sum)
  waic_indpt_hier_split[[i]] <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
  waic_indpt_hier_split[[i]]
  
  
  
  
}

var_names_parameters_indpt_hier_split_simul <- varnames(jags_indpt_hier_split_simul)

save.image("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2023-07-30_jfb_model/Simulation/outputs/_simulation_coda_R2jags_final_update_2023-07-30_M3-M4_LogN.RData")


