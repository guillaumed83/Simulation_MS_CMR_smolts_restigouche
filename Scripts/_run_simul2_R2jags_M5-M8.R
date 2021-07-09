## We run the different simulations/models in different R sessions
## To save some time 
## Run_simul_M1-M4.R
## Run_simul_M5-M8.R
## Run_simul_M9.R
## Run_simul_M10.R
## Run_simul_M11.R
##
## note that seed is te same for the data generation process: 
## _data_simulation_replicates.R
##
## and we set seed before running models




## running M5 - M8

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

source("_data_simulation_replicates2.R")
#source("C:/Users/DauphinGU/Desktop/smolts_restigouche/Model_CMR/2019-06-27_whole_season_data2019/_data_simulation.R")


set.seed(2021062)



n_burnin <- 100000
n_iter <- 3100000
n_thin <- 600




## n.sims = (n.iter - n.burnin) / n.thin
(n_iter - n_burnin)/n_thin # *n_chain




## M5
##### simul_model_CMR_indpt_hier_split_pooled.bug #####

monitor_=c('mu_theta',
           'theta',
           'sigma_theta',
           'Nm','Nm_tot','Nm_rest','p_smolt_prod',
           'mu_lambda','beta_lambda','alpha_lambda',
           'delta','eta_alphaN','split'
)


jags_indpt_hier_split_pooled_simul<-list()
samples_indpt_hier_split_pooled_simul<-list()
waic_indpt_hier_split_pooled_simul <- list()



for(i in 1:n_rep){
  
  jags_indpt_hier_split_pooled_simul[[i]]<-jags(model.file="C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Models/simul_model_CMR_indpt_hier_split.bug",
                                                data=data_pooled[[i]],
                                                inits=inits,
                                                parameters.to.save = monitor_,
                                                n.chains=2,
                                                n.burnin = n_burnin,
                                                n.iter = n_iter,
                                                n.thin= n_thin
  )
  
  #mcmcplot(jags_indpt_hier_simul )
  jags_indpt_hier_split_pooled_simul[[i]] 
  
  samples_indpt_hier_split_pooled_simul[[i]] <-jags.samples(jags_indpt_hier_split_pooled_simul[[i]]$model, c("WAIC", "deviance"),type="mean",
                                                            n.burnin = n_burnin,
                                                            n.iter = n_iter,
                                                            n.thin=n_thin)
  
  samples_indpt_hier_split_pooled_simul[[i]]$p_waic <- samples_indpt_hier_split_pooled_simul[[i]]$WAIC
  samples_indpt_hier_split_pooled_simul[[i]]$waic <- samples_indpt_hier_split_pooled_simul[[i]]$deviance + samples_indpt_hier_split_pooled_simul[[i]]$p_waic
  tmp <- sapply(samples_indpt_hier_split_pooled_simul[[i]], sum)
  waic_indpt_hier_split_pooled_simul[[i]] <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
  waic_indpt_hier_split_pooled_simul[[i]]
  
}
var_names_parameters_indpt_hier_split_pooled_simul <- varnames(jags_indpt_hier_split_pooled_simul)






## M6
##### simul_model_CMR_indpt_hier_split_other_wheels.bug #####

monitor_=c('mu_theta',
           'theta',
           'sigma_theta',
           'Nm','Nm_tot','Nm_rest','p_smolt_prod',
           'mu_lambda','beta_lambda','alpha_lambda',
           'delta','eta_alphaN','split'
)


jags_indpt_hier_split_other_wheels_simul<-list()
samples_indpt_hier_split_other_wheels_simul<-list()
waic_indpt_hier_split_other_wheels <- list()


for(i in 1:n_rep){
  jags_indpt_hier_split_other_wheels_simul[[i]]<-jags(model.file="C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Models/simul_model_CMR_indpt_hier_split_other_wheels.bug",
                                                      data=data_oth_wheels_R[[i]],
                                                      inits=inits,
                                                      parameters.to.save = monitor_,
                                                      n.chains=2,
                                                      n.burnin = n_burnin,
                                                      n.iter = n_iter,
                                                      n.thin= n_thin
  )
  
  #mcmcplot(jags_indpt_hier_other_wheels_simul[[i]] )
  jags_indpt_hier_split_other_wheels_simul[[i]] 
  
  samples_indpt_hier_split_other_wheels_simul[[i]] <-jags.samples(jags_indpt_hier_split_other_wheels_simul[[i]]$model, c("WAIC", "deviance"),type="mean",
                                                                  n.burnin = n_burnin,
                                                                  n.iter = n_iter,
                                                                  n.thin=n_thin)
  
  samples_indpt_hier_split_other_wheels_simul[[i]]$p_waic <- samples_indpt_hier_split_other_wheels_simul[[i]]$WAIC
  samples_indpt_hier_split_other_wheels_simul[[i]]$waic <- samples_indpt_hier_split_other_wheels_simul[[i]]$deviance + samples_indpt_hier_split_other_wheels_simul[[i]]$p_waic
  tmp <- sapply(samples_indpt_hier_split_other_wheels_simul[[i]], sum)
  waic_indpt_hier_split_other_wheels[[i]] <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
  waic_indpt_hier_split_other_wheels[[i]]
}

var_names_parameters_indpt_hier_split_other_wheels_simul <- varnames(jags_indpt_hier_split_other_wheels_simul)

## M7
##### simul_model_CMR_indpt_hier_no_split_other_wheels.bug #####

monitor_=c('mu_theta',
           'theta',
           'sigma_theta',
           'Nm','Nm_tot','Nm_rest','p_smolt_prod',
           'mu_lambda','beta_lambda','alpha_lambda',
           'delta','eta_alphaN','split'
)


jags_indpt_hier_no_split_other_wheels_simul<-list()
samples_indpt_hier_no_split_other_wheels_simul<-list()
waic_indpt_hier_no_split_other_wheels <- list()


for(i in 1:n_rep){
  jags_indpt_hier_no_split_other_wheels_simul[[i]]<-jags(model.file="C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Models/simul_model_CMR_indpt_hier_other_wheels.bug",
                                                         data=data_oth_wheels_R[[i]],
                                                         inits=inits,
                                                         parameters.to.save = monitor_,
                                                         n.chains=2,
                                                         n.burnin = n_burnin,
                                                         n.iter = n_iter,
                                                         n.thin= n_thin
  )
  
  #mcmcplot(jags_indpt_hier_other_wheels_simul[[i]] )
  jags_indpt_hier_no_split_other_wheels_simul[[i]] 
  
  samples_indpt_hier_no_split_other_wheels_simul[[i]] <-jags.samples(jags_indpt_hier_no_split_other_wheels_simul[[i]]$model, c("WAIC", "deviance"),type="mean",
                                                                     n.burnin = n_burnin,
                                                                     n.iter = n_iter,
                                                                     n.thin=n_thin)
  
  samples_indpt_hier_no_split_other_wheels_simul[[i]]$p_waic <- samples_indpt_hier_no_split_other_wheels_simul[[i]]$WAIC
  samples_indpt_hier_no_split_other_wheels_simul[[i]]$waic <- samples_indpt_hier_no_split_other_wheels_simul[[i]]$deviance + samples_indpt_hier_no_split_other_wheels_simul[[i]]$p_waic
  tmp <- sapply(samples_indpt_hier_no_split_other_wheels_simul[[i]], sum)
  waic_indpt_hier_no_split_other_wheels[[i]] <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
  waic_indpt_hier_no_split_other_wheels[[i]]
}

var_names_parameters_indpt_hier_no_split_other_wheels_simul <- varnames(jags_indpt_hier_no_split_other_wheels_simul)



## M8
##### simul_model_CMR_indpt_hier_merge.bug #####

monitor_=c('mu_theta',
           'theta',
           'sigma_theta',
           'Nm','Nm_tot','Nm_rest','p_smolt_prod',
           'mu_lambda','beta_lambda','alpha_lambda',
           'delta','eta_alphaN','split'
)


jags_indpt_hier_merge_simul<-list()
samples_indpt_hier_merge_simul<-list()
waic_indpt_hier_merge <- list()


for(i in 1:n_rep){
  jags_indpt_hier_merge_simul[[i]]<-jags(model.file="C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Models/simul_model_CMR_indpt_hier_merge.bug",
                                         data=data_merge[[i]],
                                         inits=inits_merge,
                                         parameters.to.save = monitor_,
                                         n.chains=2,
                                         n.burnin = n_burnin,
                                         n.iter = n_iter,
                                         n.thin= n_thin
  )
  
  #mcmcplot(jags_indpt_hier_other_wheels_simul[[i]] )
  jags_indpt_hier_merge_simul[[i]] 
  
  samples_indpt_hier_merge_simul[[i]] <-jags.samples(jags_indpt_hier_merge_simul[[i]]$model, c("WAIC", "deviance"),type="mean",
                                                     n.burnin = n_burnin,
                                                     n.iter = n_iter,
                                                     n.thin=n_thin)
  
  samples_indpt_hier_merge_simul[[i]]$p_waic <- samples_indpt_hier_merge_simul[[i]]$WAIC
  samples_indpt_hier_merge_simul[[i]]$waic <- samples_indpt_hier_merge_simul[[i]]$deviance + samples_indpt_hier_merge_simul[[i]]$p_waic
  tmp <- sapply(samples_indpt_hier_merge_simul[[i]], sum)
  waic_indpt_hier_merge[[i]] <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
  waic_indpt_hier_merge[[i]]
}

var_names_parameters_indpt_hier_merge_simul <- varnames(jags_indpt_hier_merge_simul)




var.to.keep <- c("jags_indpt_hier_split_pooled_simul",
                 "samples_indpt_hier_split_pooled_simul",
                 "waic_indpt_hier_split_pooled_simul",
                 
                 "jags_indpt_hier_split_other_wheels_simul",
                 "samples_indpt_hier_split_other_wheels_simul",
                 "waic_indpt_hier_split_other_wheels",
                 
                 "jags_indpt_hier_no_split_other_wheels_simul",
                 "samples_indpt_hier_no_split_other_wheels_simul",
                 "waic_indpt_hier_no_split_other_wheels", 
                 
                 "jags_indpt_hier_merge_simul",
                 "samples_indpt_hier_merge_simul",
                 "waic_indpt_hier_merge")


rm(list=setdiff(ls(), var.to.keep))





save.image("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_M5-M8.RData")








