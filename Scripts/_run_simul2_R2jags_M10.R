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




## running M10

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


set.seed(2021064)


n_burnin <- 500000
n_iter <- 5500000
n_thin <- 1000

## n.sims = (n.iter - n.burnin) / n.thin
(n_iter - n_burnin)/n_thin # *n_chain
# 

####  simul_model_CMR_concentration_factor_unif_Nmtot.bug #### 
####  using Dirichlet estimating the concentration dactor
monitor_=c('mu_theta',
           'theta',
           'sigma_theta',
           'Nm','Nm_tot','Nm_rest','p_smolt_prod',
           'mu_lambda','beta_lambda','alpha_lambda',
           'delta','eta_alphaN'
)

jags_concentration_factor_unif_Nmtot_simul<-list()
samples_concentration_factor_unif_Nmtot_simul<-list()
waic_concentration_factor_unif_Nmtot <- list()


for(i in 1:n_rep){
  jags_concentration_factor_unif_Nmtot_simul[[i]]<-jags(model.file="C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Models/simul_model_CMR_concentration_factor_unif_Nmtot.bug",
                                                        data=data_indpt[[i]],
                                                        inits=inits,
                                                        parameters.to.save = monitor_,
                                                        n.chains=2,
                                                        n.burnin = n_burnin,
                                                        n.iter = n_iter,
                                                        n.thin= n_thin
  )
  
  #mcmcplot(jags_concentration_factor_simul[[i]])
  jags_concentration_factor_unif_Nmtot_simul[[i]] 
  
  samples_concentration_factor_unif_Nmtot_simul[[i]]<-jags.samples(jags_concentration_factor_unif_Nmtot_simul[[i]]$model, c("WAIC", "deviance"),type="mean",
                                                                   n.burnin = n_burnin,
                                                                   n.iter = n_iter,
                                                                   n.thin=n_thin)
  
  samples_concentration_factor_unif_Nmtot_simul[[i]]$p_waic <- samples_concentration_factor_unif_Nmtot_simul[[i]]$WAIC
  samples_concentration_factor_unif_Nmtot_simul[[i]]$waic <- samples_concentration_factor_unif_Nmtot_simul[[i]]$deviance + samples_concentration_factor_unif_Nmtot_simul[[i]]$p_waic
  tmp <- sapply(samples_concentration_factor_unif_Nmtot_simul[[i]], sum)
  waic_concentration_factor_unif_Nmtot[[i]] <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
  waic_concentration_factor_unif_Nmtot[[i]]
}

var_names_parameters_concentration_factor_unif_Nmtot_simul <- varnames(jags_concentration_factor_unif_Nmtot_simul)





var.to.keep <- c("jags_concentration_factor_unif_Nmtot_simul",
                 "samples_concentration_factor_unif_Nmtot_simul",
                 "waic_concentration_factor_unif_Nmtot")


rm(list=setdiff(ls(), var.to.keep))

save.image("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_M10.RData")
