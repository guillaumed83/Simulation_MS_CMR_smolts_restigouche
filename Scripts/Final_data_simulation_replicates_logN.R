## Generating simulation dataset ----  
## various datasets are generated to resemble the REstigouche dataset
## Using similar structure than the case study: 3 upstream RSTs and 2 downstream RSTs/wheels 
## (RST and wheel are used interchangeably)
## 
##
## Author: G. Dauphin
## last update: 30-07-2021
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
## Parameters/variable : Notation            ----
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

## upper bound assuming 32 millions sqm of habitat  x 0.5 0+/m2 * 0.3 survival (Aprahamian)

EF_data <- read.csv("C:/Users/DauphinGu/Desktop/DFO projects/CSAS/A_salmon/2023/Update of indicators/Output/SFA15-16_EF_annual_d.csv")

temp <- EF_data %>% filter(Catchment=="Restigouche",Lifestage!="fry") 

temp1 <- temp %>% select(Year,Subcathment,Lifestage,Median) %>%
                  group_by(Year,Subcathment)%>% pivot_wider(names_from = Lifestage,values_from = Median) %>%
                  mutate(tot_parr =sum(c(`small parr`,`large parr`),na.rm=T) ) %>%
                  ungroup() %>%filter(Year>1971)%>% group_by(Year) %>% summarise(d_tot_parr =mean(tot_parr), n_tot_parr = mean(tot_parr) * 33.2*10^6 ) 

temp_K <- temp %>% select(Year,Subcathment,Lifestage,Median) %>%
  group_by(Year,Subcathment)%>% pivot_wider(names_from = Lifestage,values_from = Median) %>%
  mutate(tot_parr =sum(c(`small parr`,`large parr`),na.rm=T) ) %>%
  ungroup() %>%filter(Year>1971,Subcathment=="KED")%>% group_by(Year) %>% summarise(d_tot_parr =mean(tot_parr), n_tot_parr = mean(tot_parr) * 3.5*10^6 ) 

temp_U <- temp %>% select(Year,Subcathment,Lifestage,Median) %>%
  group_by(Year,Subcathment)%>% pivot_wider(names_from = Lifestage,values_from = Median) %>%
  mutate(tot_parr =sum(c(`small parr`,`large parr`),na.rm=T) ) %>%
  ungroup() %>%filter(Year>1971,Subcathment=="UPS")%>% group_by(Year) %>% summarise(d_tot_parr =mean(tot_parr), n_tot_parr = mean(tot_parr) * 3.5*10^6 ) 



mean(temp_K$n_tot_parr)
sd(temp_K$n_tot_parr)
L_mu_parr_K <- mean(log(temp_K$n_tot_parr))
L_sd_parr_K<- sd(log(temp_K$n_tot_parr))

precision_K <- 1/(L_sd_parr_K*L_sd_parr_K)

mean(temp_U$n_tot_parr)
sd(temp_U$n_tot_parr)
L_mu_parr_U <- mean(log(temp_U$n_tot_parr))
L_sd_parr_U<- sd(log(temp_U$n_tot_parr))

precision_U <- 1/(L_sd_parr_U*L_sd_parr_U)


mean(temp1$n_tot_parr*0.3)
sd(temp1$n_tot_parr*0.3)
L_mu_parr <- mean(log(temp1$n_tot_parr*0.3))
L_sd_parr<- sd(log(temp1$n_tot_parr*0.3))

precision <- 1/(L_sd_parr*L_sd_parr)
## wetted areas  
## Ked = 3.50 M
## Ups = 5.31 M
## Mat = 6.81 M
## rest =17.58 M


Tot_N <- rlnorm(n_years,L_mu_parr,L_sd_parr) #runif(n_years, 100000, 1000000 ) 
plot(1:15,Tot_N)

## Generating the proportions of smolts in each subwatershed
## eta = concentration factor
eta <- 15
alpha_dir <- rdirichlet(n_years, p * eta)
## total non-marked abundance
Tot_w<-round(Tot_N*alpha_dir)

Tot_w2 <- as.data.frame(Tot_w)
names(Tot_w2) <- c("N1","N2","N3","N4")


p_all <- Tot_w2 %>% mutate(p1 = N1 / (N1+N2+N3+N4),
                           p2 = N2 / (N1+N2+N3+N4),
                           p3 = N3 / (N1+N2+N3+N4),
                           p4 = N4 / (N1+N2+N3+N4) ) %>%
  summarise(mean_p1= mean(p1),
            mean_p2= mean(p2),
            mean_p3= mean(p3),
            mean_p4= mean(p4),
            sd_p1= sd(p1),
            sd_p2= sd(p2),
            sd_p3= sd(p3),
            sd_p4= sd(p4)
  )



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

## this is used in the model that incorporate the spliting process upstream of the 2 downstream RSTs
## probability to go throughout one of the two ds wheels

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
        ## Recaptures from the same wheel
        R_wheels[t,k] <- rbinom(1,round(p_ds[t,k] * C_wheels[t,k]),theta_wheels_list[[i]][t,k])
      }
    }
    for(k in 4:5){
      ## Recaptures from all wheels together
      M_pooled_wheels[t,k] <- sum(M_wheels[t,1:5],na.rm=T)
      R_wheels_pooled[t,k] <- rbinom(1,round(p_ds[t,k] * M_pooled_wheels[t,k]),theta_wheels_list[[i]][t,k])
    }
  }  
  L_C_wheels[[i]] <- C_wheels
  L_M_wheels[[i]] <- M_wheels
  L_M_pooled_wheels[[i]] <- M_pooled_wheels
  L_R_wheels[[i]] <- R_wheels
  L_R_pooled_wheels[[i]] <- cbind(R_wheels[,1:3],R_wheels_pooled[,4:5])
  # L_R_other_wheels_4[[i]] <- R_other_wheels_4
  # L_R_other_wheels_5[[i]] <- R_other_wheels_5

}

## Various transformation to make the dataframes suitable for data and init files in jags ####

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



## All Recaptures at each downstream RSTs pooled together #### 
## data used for models M2,M3,M4,M5


Nm_inits = 
  matrix(
         c(NA,NA,NA,
           NA,NA,NA,
           NA,NA,0,
           NA,0,0,
           NA,NA,0,
           NA,0,0,
           NA,0,NA,
           NA,NA,0,
           NA,NA,NA,
           NA,NA,NA,
           NA,NA,0,
           NA,0,0,
           NA,NA,0,
           NA,0,NA,
           NA,NA,NA),nrow=15,byrow=T
           ) 


data_pooled <- list()


for (i in 1:n_rep){
  #L_R_pooled_wheels[[i]][,1:3] <- L_R_wheels[[i]][,1:3]
  
  
  data_pooled[[i]]=list(
    Nm = Nm_inits,
    I_w = y_wheels,
    Y=n_years,
    n_p_dir = n_dir,
    I_p=I_p,
    p_smolt_prod=p_smolt_prod,
    alpha = temp_p_covered ,
    C=L_C_wheels[[i]], 
    M= L_M_pooled_wheels[[i]],
    R= L_R_pooled_wheels[[i]]
    
  )
}

## Merging the two downstream RSTs ####
## data_merge used for models M8 (data_merge2 was used in an alternate model, not presented in the MS)

I_merge <- apply(y_wheels[,4:5],1,sum)
I_merge[I_merge>0]<-1

data_merge <- list()


for (i in 1:n_rep){
  
  data_merge[[i]] =list(
    I_w = y_wheels,
    Y=n_years,
    n_p_dir = n_dir,
    I_p=I_p,
    p_smolt_prod=p_smolt_prod,
    alpha = temp_p_covered ,
    temp4=c(1,2,3,5),
    temp5=c(1,2,3,4),
    
    C=L_C_wheels[[i]], 
    M= L_M_pooled_wheels[[i]],
    R= L_R_pooled_wheels[[i]],
    I_merge = I_merge,
    C_merge = apply(L_C_wheels[[i]][,4:5],1,sum,na.rm=T),
    M_merge = apply(L_M_pooled_wheels[[i]][,1:4],1,sum,na.rm=T),
    R_merge = apply(L_R_pooled_wheels[[i]][,4:5],1,sum,na.rm=T)

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









