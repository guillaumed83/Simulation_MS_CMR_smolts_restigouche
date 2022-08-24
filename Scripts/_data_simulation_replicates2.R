## Generating simulation dataset ----  
## various datasets are generated to resemble the REstigouche dataset
## Using similar structure than the case study: 3 upstream RSTs and 2 downstream RSTs/wheels 
## (RST and wheel are used interchangeably)
## 
##
## Author: G. Dauphin
## last update: 24-08-2022
##
##

## Libraries 

library(gtools)
library(coda)
library(boot)
library(tidyverse)

## clear environment
rm(list=ls())


set.seed(8502)

##---------------------------------------------------------------------##
## Parameters/variable : Notation                                      ##
## total and wheel abundance: Tot_N, Tot_w                             ##  
## proportion of each subwatershed: alpha_dir                          ##
## concentration factor: eta                                           ##  
## index of wheel running: y_wheels                                    ##
## Wheel catchability: theta_wheels                                    ##
## CMR                                                                 ##
## Capture: C_wheels                                                   ##
## Mark: M_wheel (a proportion p_mark of the captured fish are marked) ##
## Recapture: R_wheels                                                 ##
##---------------------------------------------------------------------##



## Using similar proportion of each tributat to the REstigouche
## values used in the simulation are historic proportion of wetted area
## these values have been updated and in the case study the actual values used are:
## K = 0.105 ; U = 0.160 ; Ma = 0.205 ; rest =  0.530
p <- c(0.131,0.199,0.191,0.479)

## number of years in each replicate
n_years <- 15

## upper bound assuming 32 millions sqm of habitat  x 0.5 0+/m2 * 0.3 survival (Aprahamian et al 2003)
Tot_N <- runif(n_years, 100000, 1000000 ) 
plot(1:15,Tot_N)

## Generating the proportions of smolts in each subwatershed
## eta = concentration factor
eta <- 15
alpha_dir <- rdirichlet(n_years, p * eta)
## total non-marked abundance
Tot_w<-round(Tot_N*alpha_dir)

### matrix showing which wheels are active 
### proportions based on case study data set
y_wheels <- as.matrix(array(0,dim=c(n_years,5)))

## upstream RSTs
## 1st RST (Kedgwick) works every year 
y_wheels[,1] <- 1
## 2nd RST (Upsalquitch)
y_wheels[,2]<-rbinom(n_years, 1,.67)
## 3rd RST (Matapedia)
y_wheels[,3]<-rbinom(n_years, 1,.17)

## Downstream RSTs
## 4th and 5th RSTs (Butters and Moses)
y_wheels[,4]<-rbinom(n_years, 1,.89)
y_wheels[,5]<-rbinom(n_years, 1,.83)

## Proportion of fish marked = .85
p_mark <- 0.85

## We generate 10 dataset to do a bit of replication
## Ideally more but computing time is an issue 
## (low probability of capture = slow convergence)
n_rep <- 10

## annual wheels catchability (based on initial model estimates)
## theta in logit scale

theta_wheels_list <- list() 

for (i in 1:n_rep){
  theta_wheels <- as.matrix(array(0,dim=c(n_years,5)))
  
  theta_wheels[,1] <- rnorm(n_years,logit(0.02),0.3)
  theta_wheels[,2] <- rnorm(n_years,logit(0.025),0.6)
  theta_wheels[,3] <- rnorm(n_years,logit(0.02),0.2)
  theta_wheels[,4] <- rnorm(n_years,logit(0.0025),0.55)
  theta_wheels[,5] <- rnorm(n_years,logit(0.0025),0.35)
  
  theta_wheels <- round(inv.logit(theta_wheels),6)

## keeping only the years of interest
  theta_wheels <- theta_wheels * y_wheels 
  theta_wheels_list[[i]] <- theta_wheels
}

## this is used in the model that incorporates the spliting process upstream of the 2 downstream RSTs
## probability to go throughout one of the two downstream (ds) wheels

p_ds <- array(NA,dim=c(n_years,5))

p_ds[,4] <- runif(n_years,0.4,0.6)
p_ds[,5] <- 1-p_ds[,4]


## Captures at all wheels
## storing everything in lists

L_C_wheels <- list()
L_M_wheels <- list()
L_M_pooled_wheels <- list()
L_R_wheels <- list()
L_R_pooled_wheels <- list()
L_R_other_wheels_4 <- list()
L_R_other_wheels_5 <- list()

for (i in 1:n_rep){

  C_wheels <- as.matrix(array(NA,dim=c(n_years,5)))
  M_wheels <- as.matrix(array(NA,dim=c(n_years,5)))
  M_pooled_wheels <- as.matrix(array(NA,dim=c(n_years,5)))
  R_wheels <- as.matrix(array(NA,dim=c(n_years,5)))
  R_wheels_pooled <- as.matrix(array(NA,dim=c(n_years,5)))
  
  ## only for the downstream wheels
  ## when operating the downstream RSTs can recapture smolts marked at the upstream RSTs
  R_other_wheels_4 <- as.matrix(array(NA,dim=c(n_years,5)))
  R_other_wheels_5 <- as.matrix(array(NA,dim=c(n_years,5)))
  
  temp_4<-c(1,2,3,5)
  temp_5<-c(1,2,3,4)
  
  for(t in 1:n_years ){
    for(k in 1:3){
      if(y_wheels[t,k]>0){
        C_wheels[t,k] <- rbinom(1,Tot_w[t,k],theta_wheels_list[[i]][t,k])
        M_wheels[t,k] <- round(p_mark * C_wheels[t,k])
        M_pooled_wheels[t,k] <- round(p_mark * C_wheels[t,k])
        R_wheels[t,k] <- rbinom(1,C_wheels[t,k],theta_wheels_list[[i]][t,k])      
      }
    }
    
    ## Wheels 4 and 5 are    
    
    for(k in 4:5){
      if(y_wheels[t,k]>0){
        C_wheels[t,k] <- rbinom(1,round(p_ds[t,k] * Tot_N[t]),theta_wheels_list[[i]][t,k])
        M_wheels[t,k] <- round(p_mark * C_wheels[t,k])
        ## REcaptures from the same wheel      
        R_wheels[t,k] <- rbinom(1,round(p_ds[t,k] * C_wheels[t,k]),theta_wheels_list[[i]][t,k]) 
        
      }  
    }
    
## pooled data  
## This is for the models where the downstream RST recaptures of smolts marked at the upstream 
## are assumed to be the same process 
    
    for(k in 4:5){
      if(y_wheels[t,k]>0){
        ##C_wheels[t,k] <- rbinom(1,Tot_N[t],theta_wheels[t,k])
        ## REcaptures from the same wheel
        M_pooled_wheels[t,k] <- sum(M_wheels[t,1:5],na.rm=T)
        R_wheels_pooled[t,k] <- rbinom(1,round(p_ds[t,k] * M_pooled_wheels[t,k]),theta_wheels_list[[i]][t,k]) 
        
      }  
    }

## REcapture other wheels 
## This is for the models where the downstream RST recaptures of smolts marked at the upstream 
## are assumed to be a different process than the fish marked at a given RST and recpture at the same RST     
    if(y_wheels[t,4]>0){
      for(j in 1:4){
        if(y_wheels[t,temp_4[j]]>0){
          R_other_wheels_4[t,temp_4[j]] <- rbinom(1,round(p_ds[t,4] * M_wheels[t,temp_4[j]]),theta_wheels_list[[i]][t,4]) 
        }
      }
    }
    
    if(y_wheels[t,5]>0){
      for(j in 1:4){
        if(y_wheels[t,temp_5[j]]>0){
          R_other_wheels_5[t,temp_5[j]] <- rbinom(1,round(p_ds[t,5] * M_wheels[t,temp_5[j]]),theta_wheels_list[[i]][t,5]) 
        }
      }
    }

  }  
  L_C_wheels[[i]] <- C_wheels
  L_M_wheels[[i]] <- M_wheels
  L_M_pooled_wheels[[i]] <- M_pooled_wheels
  L_R_wheels[[i]] <- R_wheels
  L_R_pooled_wheels[[i]] <- R_wheels_pooled
  L_R_other_wheels_4[[i]] <- R_other_wheels_4
  L_R_other_wheels_5[[i]] <- R_other_wheels_5

}

## Various transforation to make the dataframes suitable for data and init files in jags ####

temp_p_covered <-as.data.frame( t(t(y_wheels[,1:3])*p[1:3]))
temp_p_covered <-temp_p_covered %>% mutate(rest=1-rowSums(.))
temp_p_covered2 <- temp_p_covered
temp_p_covered2[temp_p_covered2== 0] <- NA

temp_p_covered<-as.matrix(temp_p_covered)
dimnames(temp_p_covered) = NULL

p_covered <-as.data.frame(array(as.numeric(),dim=c(15,4)))

temp_I_p <-as.data.frame( t(t(y_wheels[,1:3])*c(1,2,3)))
temp_I_p <-temp_I_p %>% mutate(last=4)
temp_I_p[temp_I_p == 0] <- NA

I_p <-as.data.frame(array(as.numeric(),dim=c(15,4)))

n_dir <- rep(0,15)

for(i in 1:n_years){
  temp<-temp_p_covered2[i,][which(!is.na(temp_p_covered2[i,]))]
  p_covered[i,1:length(temp)]<-temp
  
  n_dir[i]<-length(which(!is.na(p_covered[i,])))
  
  temp<-temp_I_p[i,][which(!is.na(temp_I_p[i,]))]
  I_p[i,1:length(temp)]<-temp
}

I_p<- as.matrix(I_p)
dimnames(I_p) = NULL

p_covered<- as.matrix(p_covered)
dimnames(p_covered) = NULL


theta_wheels_inits1 <- theta_wheels_list[[2]] + rnorm(15*5,0,0.02) 
theta_wheels_inits1 <- theta_wheels_inits1 * y_wheels 

theta_wheels_inits1[theta_wheels_inits1==0]<-NA
theta_wheels_inits1<- abs(theta_wheels_inits1)
L_theta_wheels_inits1 <- logit(theta_wheels_inits1)


theta_wheels_inits2 <- theta_wheels_list[[7]] + rnorm(15*5,0,0.02) 
theta_wheels_inits2 <- theta_wheels_inits1 * y_wheels 

theta_wheels_inits2[theta_wheels_inits2==0]<-NA
theta_wheels_inits2<- abs(theta_wheels_inits2)
L_theta_wheels_inits2 <- logit(theta_wheels_inits2)


p_smolt_prod <-as.matrix(array(as.numeric(0),dim=c(15,4)))

to_match<-c(1,2,3,4)

for(i in 1:n_years){
  
  temp <- na.omit(match(to_match,I_p[i,]))
  for(k in 1:length(temp)){
    p_smolt_prod[i,I_p[i,k]]<-NA
  }
}

#p_smolt_prod[is.na(p_smolt_prod)]<-1

# y_wheels2 <- as.data.frame(y_wheels)
# colnames(y_wheels2) <- c("a","b","c","d","e")
# I_tot <- y_wheels2 %>% mutate(tot=max(d,e) )

##_______________________________________________________________
## Creating dataset list that are used by Jags ####

## data used for models M1, M2, M3, M9, M10
data_indpt <- list()

for (i in 1:n_rep){

    data_indpt[[i]] =list(
    I_w = y_wheels,
    Y=n_years,
    n_p_dir = n_dir,
    I_p=I_p,
    p_smolt_prod=p_smolt_prod,
    alpha = temp_p_covered ,
    temp4=c(1,2,3,5),
    temp5=c(1,2,3,4),
    
    C= L_C_wheels[[i]], 
    M= L_M_wheels[[i]],
    R= L_R_wheels[[i]]

  )

}



## Merging the two downstream RSTs ####
## data_merge used for models M8 (data_merge2 was used in an alternate model, not presented in the MS)

I_merge <- apply(y_wheels[,4:5],1,sum)
I_merge[I_merge>0]<-1


data_merge <- list()
data_merge2 <- list()

temp <- L_R_wheels[[i]][,5]
temp[is.na(temp)]<-0


for (i in 1:n_rep){
  
  temp1 <- L_R_wheels[[i]][,4]
  temp1[is.na(temp1)]<-0

  temp2 <- L_R_wheels[[i]][,5]
  temp2[is.na(temp2)]<-0

  data_merge[[i]] =list(
    I_w = y_wheels,
    Y=n_years,
    n_p_dir = n_dir,
    I_p=I_p,
    p_smolt_prod=p_smolt_prod,
    alpha = temp_p_covered ,
    temp4=c(1,2,3,5),
    temp5=c(1,2,3,4),
    
    C=L_C_wheels[[i]], R= L_R_wheels[[i]],M= L_M_wheels[[i]],
    I_merge = I_merge,
    C_merge = apply(L_C_wheels[[i]][,4:5],1,sum,na.rm=T),
    M_merge = apply(L_M_wheels[[i]][,1:5],1,sum,na.rm=T),
    R_merge = temp1+
              temp2+ 
              apply(L_R_other_wheels_4[[i]],1,sum,na.rm=T) + apply(L_R_other_wheels_5[[i]],1,sum,na.rm=T)
    
    
  )
  
  data_merge2[[i]] =list(
    I_w = y_wheels,
    Y=n_years,
    n_p_dir = n_dir,
    I_p=I_p,
    p_smolt_prod=p_smolt_prod,
    alpha = temp_p_covered ,
    temp4=c(1,2,3,5),
    temp5=c(1,2,3,4),
    
    C=L_C_wheels[[i]], R= L_R_wheels[[i]],M= L_M_wheels[[i]],
    I_merge = I_merge,
    C_merge = apply(L_C_wheels[[i]][,4:5],1,sum,na.rm=T),
    M_merge = apply(L_M_wheels[[i]][,4:5],1,sum,na.rm=T),
    R_merge = temp1+
              temp2
    
    
  )

}

## Recaptures from other RSTs at the downstream RSTs treated as a different process ####
## data used for models M6, M11
data_oth_wheels_R <- list()


for (i in 1:n_rep){
  data_oth_wheels_R[[i]]=list(

    I_w = y_wheels,
    Y=n_years,
    n_p_dir = n_dir,
    I_p=I_p,
    p_smolt_prod=p_smolt_prod,
    alpha = temp_p_covered ,
    C=L_C_wheels[[i]], 
    M= L_M_wheels[[i]],
    R= L_R_wheels[[i]],
    R_other_4=L_R_other_wheels_4[[i]],
    R_other_5=L_R_other_wheels_5[[i]],
    
    temp4=c(1,2,3,5),
    temp5=c(1,2,3,4)
    
    )
}


## All Recaptures at each downstream RSTs pooled together #### 
## data used for models M4, M5
data_pooled <- list()


for (i in 1:n_rep){
  L_R_pooled_wheels[[i]][,1:3] <- L_R_wheels[[i]][,1:3]
  
  data_pooled[[i]]=list(
    
    I_w = y_wheels,
    Y=n_years,
    n_p_dir = n_dir,
    I_p=I_p,
    p_smolt_prod=p_smolt_prod,
    alpha = temp_p_covered ,
    C=L_C_wheels[[i]], 
    M= L_M_pooled_wheels[[i]],
    R= L_R_pooled_wheels[[i]],
    
    temp4=c(1,2,3,5),
    temp5=c(1,2,3,4)
    
  )
}


##________________________________________________________
#### Generating the Inits for 2 chains                ####
## The merge dataset requires slighlty modified inits

delta_inits <- p_covered
#delta_inits[is.na(delta_inits)]<-0.5

inits =list( 
  list(
    mu_theta = rep(-0.2,5),
    sigma_theta = rep(0.1,5),
    mu_lambda = 400000,
    beta_lambda = 0.001,
    lambda_Nm_tot = rep(400000,n_years),
    
    eta_alphaN=10,
    
    delta = delta_inits,#I_p,
    
    L_theta = L_theta_wheels_inits1
  ),
  list(
    mu_theta = rep(5,5),
    sigma_theta = rep(5,5),
    mu_lambda = 100000,
    beta_lambda = 0.5,
    lambda_Nm_tot = rep(100000,n_years),
    
    eta_alphaN=1,
    
    delta = delta_inits,#I_p,
    
    L_theta = L_theta_wheels_inits2
  )  
 
  
)

inits_merge =list( 
  list(
    mu_theta = rep(-0.2,5),
    sigma_theta = rep(0.1,5),
    mu_lambda = 400000,
    beta_lambda = 0.001,
    lambda_Nm_tot = rep(400000,n_years),
    
    eta_alphaN=10,
    
    delta = delta_inits,#I_p,
    
    L_theta = L_theta_wheels_inits1[,1:4]
  ),
  list(
    mu_theta = rep(5,5),
    sigma_theta = rep(5,5),
    mu_lambda = 100000,
    beta_lambda = 0.5,
    lambda_Nm_tot = rep(100000,n_years),
    
    eta_alphaN=1,
    
    delta = delta_inits,#I_p,
    
    L_theta = L_theta_wheels_inits2[,1:4]
  )  
  
  
  
)









