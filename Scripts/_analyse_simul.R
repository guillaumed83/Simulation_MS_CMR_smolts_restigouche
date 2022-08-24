
setwd("C:/Users/DauphinGU/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Simulation")



rm(list=ls())


library(data.table) ## to be able to collapse lists

library(tidyverse) ## data manipulation
library(MCMCvis)
library(mcmcplots)
library(rlist)
library(DirichletReg)
library(ggridges)
library(ggpubr)
library(scales)
theme_set(theme_minimal())

#load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_split.RData")


load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_M1-M4.RData")
load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_M5-M8.RData")
load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_M8b.RData")
load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_M9.RData")
load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_M10.RData")
load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_M11.RData")


K <- runif(10000,0.01,10000)

temp1 <- K %*% t(p)
par(mfrow=c(2,2))

for (i in 1:4){
  plot(density(temp1[,i]))
}




temp<- rdirichlet(10000,K %*% t(p) )
par(mfrow=c(2,2))

for (i in 1:4){
  plot(density(temp[,i]))
}

apply(temp,2, mean)



par(mfrow=c(2,4))
temp<- rdirichlet(10000,0.01*p)

for (i in 1:4){
  plot(density(temp[,i]))
}
temp<- rdirichlet(10000,15*p)

for (i in 1:4){
  plot(density(temp[,i]))
}



temp<- rdirichlet(10000,0.1*p)
par(mfrow=c(2,2))
for (i in 1:4){
  plot(density(temp[,i]))
}
for (i in 1:4){
plot(1:10000,temp[,i])
}


for (i in 1:4){
  plot(1:10000,temp[,i])
}

temp<- rdirichlet(10000,2*p)
par(mfrow=c(2,2))
for (i in 1:4){
  plot(density(temp[,i]))
}

temp<- rdirichlet(10000,5*p)
par(mfrow=c(2,2))
for (i in 1:4){
  plot(density(temp[,i]))
}
temp<- rdirichlet(10000,10*p)
par(mfrow=c(2,2))
for (i in 1:4){
  plot(density(temp[,i]))
}

temp<- rdirichlet(10000,15*p)
par(mfrow=c(2,2))
for (i in 1:4){
  plot(density(temp[,i]))
}
for (i in 1:4){
  plot(1:10000,temp[,i])
}

temp<- rdirichlet(10000,50*p)
par(mfrow=c(2,2))
for (i in 1:4){
  plot(density(temp[,i]))
}



waic_1 <- as.data.frame(list.rbind(waic_indpt_no_hier)) %>% mutate(ID=1,replicate=rownames(.))
waic_2 <- as.data.frame(list.rbind(waic_indpt_hier)) %>% mutate(ID=2,replicate=rownames(.))
waic_3 <- as.data.frame(list.rbind(waic_indpt_hier_split)) %>% mutate(ID=3,replicate=rownames(.))
waic_4 <- as.data.frame(list.rbind(waic_indpt_hier_pooled)) %>% mutate(ID=4,replicate=rownames(.))
waic_5 <- as.data.frame(list.rbind(waic_indpt_hier_split_pooled_simul)) %>% mutate(ID=5,replicate=rownames(.))
waic_6 <- as.data.frame(list.rbind(waic_indpt_hier_split_other_wheels)) %>% mutate(ID=6,replicate=rownames(.))
waic_7 <- as.data.frame(list.rbind(waic_indpt_hier_no_split_other_wheels)) %>% mutate(ID=7,replicate=rownames(.))
waic_8 <- as.data.frame(list.rbind(waic_indpt_hier_merge)) %>% mutate(ID=8,replicate=rownames(.))
waic_8b <- as.data.frame(list.rbind(waic_indpt_hier_merge2)) %>% mutate(ID=9,replicate=rownames(.))
waic_9 <- as.data.frame(list.rbind(waic_informative_dir_unif_Nmtot)) %>% mutate(ID=10,replicate=rownames(.))
waic_10 <- as.data.frame(list.rbind(waic_concentration_factor_unif_Nmtot)) %>% mutate(ID=11,replicate=rownames(.))
waic_11 <- as.data.frame(list.rbind(waic_concentration_factor_other_R_unif_Nmtot)) %>% mutate(ID=12,replicate=rownames(.))


all_waic <- rbind(waic_1,waic_2,waic_3, waic_4,waic_5,waic_6, waic_7,waic_8,waic_8b,waic_9,waic_10,waic_11)
all_waic %>% group_by(ID) %>%summarise(MEAN=mean(waic),SD=sd(waic),MEAN_p = mean(p_waic),SD_p = sd(p_waic))  

par(mfrow=c(1,1))
pretty_y <- pretty(all_waic$waic)

plot(0,0,type="n",
     xlab="Replicates",ylab="waic",
     xlim=c(0,10),ylim=c(min(pretty_y),max(pretty_y)))#c(600,800)) #c(min(pretty_y),max(pretty_y)))#

col_ <- c("darkorange4","darkorange3","darkorange2","darkorange",
          "goldenrod1","goldenrod2","goldenrod3",
          "firebrick3","firebrick4"  ,
          "skyblue1","skyblue3","skyblue4" ) #"goldenrod3",


for (i in 1:length(col_)){
  temp <- all_waic %>% filter(ID == i)
  points(temp$replicate,sort(temp$waic),type="o",pch=16,col=col_[i])
  
  abline(h=mean(temp$waic),lty=2,lwd=2,col=col_[i])
  
}
legend("topleft",legend=c("indpt_no_hier","indpt_hier","indpt_hier_split","indpt_hier_pooled",
                          "indpt_hier_split_pooled","indpt_hier_split_other_wheels", "indpt_hier_no_split_other_wheels",
                          "indpt_hier_merge","indpt_hier_merge2",
                          "informative_dir","concentration_factor","concentration_factor_other_R"),
       pch=rep(16,11),col=col_)


mean(waic_1[,1])  
mean(waic_2[,1])
mean(waic_3[,1])
mean(waic_4[,1])
mean(waic_5[,1])
mean(waic_6[,1])
mean(waic_7[,1])
mean(waic_8[,1])
mean(waic_8b[,1])
mean(waic_9[,1])
mean(waic_10[,1])
mean(waic_11[,1])


source("simul_helper.R")


## We drop model 8b (merge of the 2 downstream RSTs using only the downstream recaptures)
## 


### Keeping the simulations we are interested in
simul_names <- c("jags_indpt_no_hier_simul",
                 "jags_indpt_hier_simul",
                 "jags_indpt_hier_split_simul",
                 "jags_indpt_hier_simul_pooled",
                 "jags_indpt_hier_split_pooled_simul",
                 "jags_indpt_hier_split_other_wheels_simul",
                 "jags_indpt_hier_no_split_other_wheels_simul",
                 "jags_indpt_hier_merge_simul",
                 #"jags_indpt_hier_merge2_simul",
                 "jags_informative_dir_unif_Nmtot_simul",
                 "jags_concentration_factor_unif_Nmtot_simul",
                 "jags_concentration_factor_other_R_unif_Nmtot_simul"
                 
                 
                 )#<-  ls()[grepl("jags_", ls())] 


temp_names <-character(0)
for(i in 1:15){
  temp_names[i]<-paste("Nm_tot[",i,"]",sep="")
}



png("Figures/Simulation_total_abundance.png",height=1300,width=1800,pointsize = 32)
#par(mfrow=c(2,4),mar=c(2,2,3,2),oma=c(2,4,2,1))
par(mfrow=c(3,4),mar=c(2,2,3,2),oma=c(2,4,2,1))

## M1
## all wheels independant no hierarchy 
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M1",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_no_hier_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M2
## all wheels independant annual hierarchy
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M2",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_hier_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M3
## all wheels independant annual hierarchy + split
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M3",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_hier_split_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M4
## all wheels independant annual hierarchy + pooled
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M4",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_hier_simul_pooled,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M5
## all wheels independant annual hierarchy + split +pooled
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M5",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_hier_split_pooled_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)


## M6
## all wheels independant annual hierarchy + split + others
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M6",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_hier_split_other_wheels_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M7
## all wheels independant annual hierarchy + no split + others
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M7",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_hier_no_split_other_wheels_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M8
## all wheels independant annual hierarchy - bottom wheels merged using all recaptures including upstream RST
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M8",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_hier_merge_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M8b
## all wheels independant annual hierarchy - bottom wheels merged only using downstream recaptures
# plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M8b",
#      xlab="",ylab="",axes=F)
# axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
# axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
# abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
# box()
# 
# for(k in 1:3){
#   for (i in 1:10){
#     plot_sim(jags_indpt_hier_merge2_simul,i,k)
#   }
# }
# points(1:15,Tot_N,col="red",type="b",pch=16)




#######################################
## M9
## Informative dirichlet Dirichlet 
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M9",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_informative_dir_unif_Nmtot_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M10
## Concentration factor Dirichlet

plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M10",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_concentration_factor_unif_Nmtot_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)


## M11
## Concentration factor Dirichlet - downstream tags

plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,5000000),pch=16,type="n",main="M11",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_concentration_factor_other_R_unif_Nmtot_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)
# 
# 
# ## Concentration factor Dirichlet - uniform on Nm_to
# 
# plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,3500000),pch=16,type="n",main="M6",
#      xlab="",ylab="",axes=F)
# axis(2,at=seq(0,3500000,500000),labels=seq(0,3500,500),las=2)
# axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
# abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
# box()
# 
# for(k in 1:3){
#   for (i in 1:10){
#     plot_sim(jags_concentration_factor_unif_Nmtot_simul,i,k)
#   }
# }
# points(1:15,Tot_N,col="red",type="b",pch=16)
# 


# plot.new()
# legend("center",legend=c("True value","Estimated median","Estimated 25th-75th percentiles","Estimated 2.5th-97.5th percentiles" ),
#        cex=0.9,bty="n",pt.cex=1.3,
#        pch=c(16,16,15,15),col=c("red",rgb(50,50,50,alpha=125,maxColorValue = 255),rgb(90,180,172,alpha=125,maxColorValue = 255),rgb(216,179,101,alpha=125,maxColorValue = 255))
# )


mtext( "Total abundance (x1000)",2,out=T,line=2)
mtext( "Years",1,out=T,line=1)


#points(1:15,ifelse(y_wheels[,4]==0|y_wheels[,5]==0,100000,NA),pch="*",cex=1.2)
dev.off()


### CVs
## We drop model 8b (merge of the 2 downstream RSTs using only the downstream recaptures)
## 



CV_Nm_tot <- list()


all_simul <- list(jags_indpt_no_hier_simul,
                  jags_indpt_hier_simul,
                  jags_indpt_hier_split_simul,
                  jags_indpt_hier_simul_pooled,
                  jags_indpt_hier_split_pooled_simul,
                  jags_indpt_hier_split_other_wheels_simul,
                  jags_indpt_hier_no_split_other_wheels_simul,
                  jags_indpt_hier_merge_simul,
                  #jags_indpt_hier_merge2_simul,
                  jags_informative_dir_unif_Nmtot_simul,
                  jags_concentration_factor_unif_Nmtot_simul,
                  jags_concentration_factor_other_R_unif_Nmtot_simul)

for (n in 1:length(simul_names)){
  temp_array <- array(0,dim=c(n_years,n_rep))
  
  for (i in 1:n_rep){
    temp_sum <- as.data.frame(
      all_simul[[n]][[i]]$BUGSoutput$summary[row.names(all_simul[[n]][[i]]$BUGSoutput$summary) %in% temp_names,])
    
    temp_array[,i] <- temp_sum$sd / temp_sum$mean
    
  }
  CV_Nm_tot[[simul_names[n]]]<-temp_array
}

for (n in 1:length(simul_names)){
  print(max(CV_Nm_tot[[n]]) )
}


png("Figures/Simulation_CV_total_abundance.png",height=1200,width=1800,pointsize = 32)
par(mfrow=c(3,4),mar=c(2,2,3,2),oma=c(2,4,2,1))
for (n in 1:length(simul_names)){
  plot(rep(1:15,10),CV_Nm_tot[[n]],col=rgb(30,130,181,alpha=75, maxColorValue = 255),ylim=c(0,2.5),pch=16,main=paste("M",n,sep=""),yaxt="none",type="n")
  axis(2,at=seq(0,4,0.2),labels=seq(0,4,0.2),las=2)
  abline(h=seq(0,2.5,0.5),lty=2,col="grey85")
  points(rep(1:15,10),CV_Nm_tot[[n]],col=rgb(30,130,181,alpha=75, maxColorValue = 255),ylim=c(0,1),pch=16)
  abline(h=apply(CV_Nm_tot[[n]],2,mean),col=rgb(30,130,181,alpha=170, maxColorValue = 255),lty=2)
  abline(h=mean(CV_Nm_tot[[n]]),col="red",lty=2,lwd=2)
  
  text(11,2.2,paste("Avg. CV = ",round(mean(CV_Nm_tot[[n]]),3)))
  
}  
mtext( "Coefficient of Variation (std. dev. / mean)",2,out=T,line=2)

dev.off()



####  difference obs - pred - Nm_tot ####
png("Figures/Simulation_diff_total_abundance.png",height=1200,width=1800,pointsize = 32)
list_diff_Tot_N <- list()
par(mfrow=c(3,4),mar=c(2,2,3,2),oma=c(2,4,2,1))
for (n in 1:length(simul_names)){
  list_diff_Tot_N[[simul_names[n]]] <- diff_obs_pred(simul_names[n],"Tot_N")
  
  #y_axis <- pretty(as.matrix( rbindlist(list_diff_Tot_N[[simul_names[n]]][1:10])[,1:5] )    )
  y_axis <-seq(-5000000,1000000,1000000)#c(-5000000,-2000000,-1000000,0,1000000)
  plot(1:15,Tot_N,xlim=c(1,15),ylim=c(min(y_axis),max(y_axis)),pch=16,type="n",main=paste("M",n,sep=""),
  xlab="",ylab="",axes=F)
  axis(2,at=seq(-5000000,1000000,500000),labels=seq(-5000,1000,500),las=2)
  axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
  
  abline(h=seq(-5000000,1000000,500000),lty=2,col="grey85")
  
  for(k in 1:3){
    for (i in 1:10){
      plot_sim(list_diff_Tot_N[[simul_names[n]]],i,k)
    }
  }
  abline(h=0,lty=2,lwd=2,col="blue")
  abline(h=mean(dplyr::bind_rows(list_diff_Tot_N[[simul_names[n]]])$mean),lty=2,lwd=2,col="red")
  
  text(10,-4000000, labels = paste("Avg. Diff. = ",round(mean(dplyr::bind_rows(list_diff_Tot_N[[simul_names[n]]])$mean)) )   )
  
  
}

# plot.new()
# legend("center",legend=c("True value","Estimated median","Estimated 25th-75th percentiles","Estimated 2.5th-97.5th percentiles" ),
#        cex=0.9,bty="n",pt.cex=1.3,
#        pch=c(16,16,15,15),col=c("red",rgb(50,50,50,alpha=125,maxColorValue = 255),rgb(90,180,172,alpha=125,maxColorValue = 255),rgb(216,179,101,alpha=125,maxColorValue = 255))
# )

mtext( "Total abundance (x1000)",2,out=T,line=2)
mtext( "Years",1,out=T,line=1)

dev.off()




list_diff_tot_N_plot <- list()
gg_list <- list()

for (n in 1:length(simul_names)){
  
  list_diff_tot_N_plot[[n]] <- diff_obs_pred2(simul_names[n],"Tot_N")
}  



for( i in 1:length(all_simul)){
  
  temp_list <- list()

  for(j in 1:n_rep){
    
  temp <- as.data.frame(list_diff_tot_N_plot[[i]][[j]])
  names(temp)<-seq(1,15,1)
  temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
  temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
  temp_list[[j]] <- temp
  
  }
  
  
  
  
  gg1 <- ggplot(temp_list[[1]], aes(x=value,y=Year)) +
          scale_x_continuous(limits=c(-5000000,1000000))+
          geom_vline(xintercept = mean(temp_list[[1]]$value),col="#0072B250",linetype="dashed",size=1,alpha=0.2)+
          geom_density_ridges(scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  
  for(k in 2:10){
    
    gg1 <- gg1 + 
      geom_vline(xintercept = mean(temp_list[[k]]$value),col="#0072B250",linetype="dashed",size=1,alpha=0.2)+
      geom_density_ridges(data=temp_list[[k]], aes(x=value,y=Year),scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  }

  gg1 <- gg1 + 
    #geom_point(data=dat_tot_N,mapping=aes(x=0, y= Year),col="red",size=1.5) + 
    geom_vline(xintercept = 0,col="grey65",linetype="dashed",size=1.5)+
    geom_vline(xintercept = mean(dplyr::bind_rows(temp_list)$value),col="red",linetype="dashed",size=1.5  ) + 
    annotate("text",x=-Inf,y=14.5, label = paste("Avg. Diff. = ",round(mean(dplyr::bind_rows(temp_list)$value))),  hjust = -0.1) +
    xlab("diff. Total abundance") +
    #xlim(0,2.5)+
    theme(legend.position="none")
  
  gg_list[[i]] <- gg1
  

}

gg_list[[1]]


##### calculation (True value - estimate)/ True value ####


diff_pcent_Nm_tot <- list()

summary_ratio_diff_true <- data.frame(

  mean = rep(NA,length(simul_names)),
  sd = rep(NA,length(simul_names)),
  min = rep(NA,length(simul_names)),
  max = rep(NA,length(simul_names))
)

for (n in 1:length(simul_names)){
  temp_array <- array(0,dim=c(n_years,n_rep))

  for (i in 1:n_rep){
    
   temp <- sweep(-all_simul[[n]][[i]]$BUGSoutput$sims.list$Nm_tot,2, Tot_N,'+')
   
   
   temp <-  sweep(temp,2,Tot_N,'/')

   temp_array[,i]<- apply(temp,2,mean)    
    


  }
  summary_ratio_diff_true[n,] <- c(
    mean(temp_array),
    sd(temp_array),
    min(temp_array),
    max(temp_array)
  )
  
  
  diff_pcent_Nm_tot[[simul_names[n]]]<-temp_array
}

write.csv(round(summary_ratio_diff_true,2),"outputs/ratio_diff_true.csv")


####



tiff("Figures/diff_Nm_tot_ggplot.tiff", units="in", width=20, height=15, res=250)

ggarrange(gg_list[[1]],
          gg_list[[2]],
          gg_list[[3]],
          gg_list[[4]],
          gg_list[[5]],
          gg_list[[6]],
          gg_list[[7]],
          gg_list[[8]],
          gg_list[[9]],
          gg_list[[10]],
          gg_list[[11]],
          labels = c("M1", "M2","M3","M4","M5","M6","M7","M8","M9","M10","M11"),
          ncol = 4, nrow = 3) +
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm")) 

dev.off()



###################
### first wheel ###
###################


temp_names <-character(0)
for(i in 1:15){
  temp_names[i]<-paste("Nm[",i,",1]",sep="")
}


par(mfrow=c(3,4))



for( n in 1:length(all_simul)){
  plot(1:15,Tot_w[,1],xlim=c(1,15),ylim=c(0,400000),pch=16,type="n",main=simul_names[n])
  abline(h=seq(0,400000,100000),lty=2,col="grey85")
  
  for(k in 1:3){
    for (i in 1:10){
      plot_sim(all_simul[[n]],i,k)
    }
  }
  points(1:15,Tot_w[,1],col="red",type="b",pch=16)

}


## difference obs - pred - Nm_tot

list_diff_N_1 <- list()
par(mfrow=c(3,4))
for (n in 1:length(simul_names)){
  list_diff_N_1[[simul_names[n]]] <- diff_obs_pred(simul_names[n],"N_1")
  
  #y_axis <- pretty(as.matrix( rbindlist(list_diff_Tot_N[[simul_names[n]]][1:10])[,1:5] )    )
  y_axis <-c(-100000,100000)
  plot(1:15,Tot_w[,1],xlim=c(1,15),ylim=c(min(y_axis),max(y_axis)),pch=16,type="n",main=simul_names[n])
  for(k in 1:3){
    for (i in 1:10){
      plot_sim(list_diff_N_1[[simul_names[n]]],i,k)
    }
  }
  abline(h=0,lty=2,lwd=2,col="red")
}


list_diff_N1_plot <- list()
gg_N1_list <- list()

for (n in 1:length(simul_names)){
  
  list_diff_N1_plot[[n]] <- diff_obs_pred2(simul_names[n],"N_1")
}  



for( i in 1:length(all_simul)){
  
  temp_list <- list()
  
  for(j in 1:n_rep){
    
    temp <- as.data.frame(list_diff_N1_plot[[i]][[j]])
    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    temp_list[[j]] <- temp
    
  }
  
  
  
  
  gg1_N1 <- ggplot(temp_list[[1]], aes(x=value,y=Year)) +
    scale_x_continuous(limits=c(-100000,100000))+
    geom_vline(xintercept = mean(temp_list[[1]]$value),col="#0072B250",linetype="dashed",size=1,alpha=0.2)+
    geom_density_ridges(scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  
  for(k in 2:10){
    
    gg1_N1 <- gg1_N1 + 
      geom_vline(xintercept = mean(temp_list[[k]]$value),col="#0072B250",linetype="dashed",size=1,alpha=0.2)+
      geom_density_ridges(data=temp_list[[k]], aes(x=value,y=Year),scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  }
  
  gg1_N1 <- gg1_N1 + 
    #geom_point(data=dat_tot_N,mapping=aes(x=0, y= Year),col="red",size=1.5) + 
    geom_vline(xintercept = 0,col="grey65",linetype="dashed",size=1.5)+
    geom_vline(xintercept = mean(dplyr::bind_rows(temp_list)$value),col="red",linetype="dashed",size=1.5  ) + 
    annotate("text",x=-Inf,y=14.5, label = paste("Avg. Diff. = ",round(mean(dplyr::bind_rows(temp_list)$value))),  hjust = -0.1) +
    xlab("diff. Total abundance") +
    #xlim(0,2.5)+
    theme(legend.position="none")
  
  gg_N1_list[[i]] <- gg1_N1
  
  
}

gg_N1_list[[2]]




##### calculation (True value - estimate)/ True value ####


diff_pcent_N1 <- list()

summary_ratio_diff_true_N1 <- data.frame(
  
  mean = rep(NA,length(simul_names)),
  sd = rep(NA,length(simul_names)),
  min = rep(NA,length(simul_names)),
  max = rep(NA,length(simul_names))
)

for (n in 1:length(simul_names)){
  temp_array <- array(0,dim=c(n_years,n_rep))
  
  for (i in 1:n_rep){

    temp <- sweep(-all_simul[[n]][[i]]$BUGSoutput$sims.list$Nm[,,1],2, Tot_w[,1],'+')
    
    
    temp <-  sweep(temp,2,Tot_w[,1],'/')
    
    temp_array[,i]<- apply(temp,2,mean)    
    
    
    
  }
  summary_ratio_diff_true_N1[n,] <- c(
    mean(temp_array),
    sd(temp_array),
    min(temp_array),
    max(temp_array)
  )
  
  
  diff_pcent_N1[[simul_names[n]]]<-temp_array
}

write.csv(round(summary_ratio_diff_true_N1,2),"outputs/ratio_diff_true_N1.csv")


####



tiff("Figures/diff_N1_ggplot.tiff", units="in", width=20, height=15, res=250)

ggarrange(gg_N1_list[[1]],
          gg_N1_list[[2]],
          gg_N1_list[[3]],
          gg_N1_list[[4]],
          gg_N1_list[[5]],
          gg_N1_list[[6]],
          gg_N1_list[[7]],
          gg_N1_list[[8]],
          gg_N1_list[[9]],
          gg_N1_list[[10]],
          gg_N1_list[[11]],
          labels = c("M1", "M2","M3","M4","M5","M6","M7","M8","M9","M10","M11"),
          ncol = 4, nrow = 3) +
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm")) 

dev.off()





###___________________________________________________________________________
### Proportion of each wheel   ####


source("simul_helper.R")

p_w <- as.data.frame(Tot_w)
names(p_w) <- c("W_1","W_2","W_3","W_4")

p_w <-p_w %>% mutate(p_1 = round( W_1 /rowSums(select(.,starts_with("W_"))), 3),
                     p_2 = round( W_2 /rowSums(select(.,starts_with("W_"))), 3),
                     p_3 = round( W_3 /rowSums(select(.,starts_with("W_"))), 3),
                     p_4 = round( W_4 /rowSums(select(.,starts_with("W_"))), 3))

p_w$Year <- seq(1,15,1)


p1 <- as.data.frame(all_simul[[2]][[1]]$BUGSoutput$sims.list$Nm[,,1] / all_simul[[2]][[1]]$BUGSoutput$sims.list$Nm_tot)
p2 <- as.data.frame(all_simul[[2]][[1]]$BUGSoutput$sims.list$Nm[,,2] / all_simul[[2]][[1]]$BUGSoutput$sims.list$Nm_tot)
years2keep <- y_wheels[,2]*1:15
years2keep <- years2keep[years2keep!=0]
yy <- 1:15
yy <- yy[!yy %in% years2keep]

for(i in 1:length(yy)){
  p2[,yy[i]] <- 0
  
}

p3 <- as.data.frame(all_simul[[2]][[1]]$BUGSoutput$sims.list$Nm[,,3] / all_simul[[2]][[1]]$BUGSoutput$sims.list$Nm_tot)
years2keep <- y_wheels[,3]*1:15
years2keep <- years2keep[years2keep!=0]
yy <- 1:15
yy <- yy[!yy %in% years2keep]

for(i in 1:length(yy)){
  p3[,yy[i]] <- 0
  
}


p_all <-  p1+p2+p3
names(p_all)<-seq(1,15,1)
p_all <- p_all %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
p_all$Year<-factor(p_all$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))


 ggplot(p_all, aes(x=value,y=Year)) +
  geom_density_ridges(scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")+
   geom_vline(aes(xintercept = 1),col="red",linetype="dashed")



for( n in 1:length(all_simul)){
  # png(paste("Figures/proportions_",simul_names[i],".png",sep=""),
  #     height=1200,width=1200,pointsize =16)
  # 
  
  tiff(paste("Figures/proportions_",simul_names[n],".tiff",sep=""),
       units="in", width=7, height=7, res=300)
  
  print(plot_p_smolt(all_simul[[n]],n))
  
  dev.off()
}

# plot_p_smolt(jags_concentration_factor_other_R_unif_Nmtot_simul)
# plot_p_smolt(jags_indpt_hier_simul)
# 
# plot_p_smolt(jags_indpt_hier_split_simul)
# plot_p_smolt(jags_indpt_hier_merge_simul)


#### ridge plot Nm tot ####

dat_tot_N <-data.frame(Year=1:15,tot_N=Tot_N)

gg_list <-list()

for( i in 1:length(all_simul)){
  
  list_p_simul <- list()
  
  for(n in 1:n_rep){
    
    temp <- as.data.frame(all_simul[[i]][[n]]$BUGSoutput$sims.list$Nm_tot)
    
    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    list_p_simul[[n]] <- temp
    
  }
  
  gg1 <-  ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+#,aes(height=1))+#..ndensity..))+#+
    scale_x_continuous(breaks=c(0,1000000,2000000,3000000,4000000,5000000),labels=c(0,1000,2000,3000,4000,5000))+
    labs(x="Total abundance (x1000)")
      #labels = unit_format(unit = "M", scale = 1e-6))
  

   
  for(n in 2:10){
    
    gg1 <- gg1 + geom_density_ridges(data=list_p_simul[[n]], aes(x=value,y=Year),rel_min_height = 0.001,
                                     col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)
  }
  
  gg1 <- gg1 + geom_point(data=dat_tot_N,mapping=aes(x=tot_N, y= Year),col="red",size=1.5) + 
    xlab("Total abundance") +
    #xlim(0,2.5)+
    theme(legend.position="none")
  
  gg_list[[i]] <- gg1

  
}

gg_all <- ggarrange(gg_list[[1]],
          gg_list[[2]],
          gg_list[[3]],
          gg_list[[4]],
          gg_list[[5]],
          gg_list[[6]],
          gg_list[[7]],
          gg_list[[8]],
          gg_list[[9]],
          gg_list[[10]],
          gg_list[[11]],
           labels = c("M1", "M2","M3","M4","M5","M6","M7","M8","M9","M10","M11"),
          ncol = 4, nrow = 3) +
        theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm")) 



tiff("Figures/Nm_tot_ggplot.tiff", units="in", width=20, height=15, res=250)

gg_all

dev.off()



#### ridge plot N1 ####

dat_tot_N <-data.frame(Year=1:15,tot_N=Tot_w[,1])

gg_list <-list()

for( i in 1:length(all_simul)){
  
  list_p_simul <- list()
  
  for(n in 1:n_rep){
   
    temp <- as.data.frame(all_simul[[i]][[n]]$BUGSoutput$sims.list$Nm[,,1])
    
    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    list_p_simul[[n]] <- temp
    
  }
  
  gg1 <-  ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+#,aes(height=1))+#..ndensity..))+#+
    scale_x_continuous(breaks=c(0,250000,500000),labels=c(0,250,500))+
    xlim(0,500000)+
    labs(x="Abundance wheel 1 (x1000)")
  #labels = unit_format(unit = "M", scale = 1e-6))
  
  
  
  for(n in 2:10){
    
    gg1 <- gg1 + geom_density_ridges(data=list_p_simul[[n]], aes(x=value,y=Year),rel_min_height = 0.001,
                                     col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)
  }
  
  gg1 <- gg1 + geom_point(data=dat_tot_N,mapping=aes(x=tot_N, y= Year),col="red",size=1.5) + 
    xlab("Abundance wheel 1") +
    #xlim(0,2.5)+
    theme(legend.position="none")
  
  gg_list[[i]] <- gg1
  
  
}

gg_all <- ggarrange(gg_list[[1]],
                    gg_list[[2]],
                    gg_list[[3]],
                    gg_list[[4]],
                    gg_list[[5]],
                    gg_list[[6]],
                    gg_list[[7]],
                    gg_list[[8]],
                    gg_list[[9]],
                    gg_list[[10]],
                    gg_list[[11]],
                    labels = c("M1", "M2","M3","M4","M5","M6","M7","M8","M9","M10","M11"),
                    ncol = 4, nrow = 3) +
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm")) 



tiff("Figures/N1_ggplot.tiff", units="in", width=20, height=15, res=250)

gg_all

dev.off()




















# ggsave("Figures/Nm_tot_ggplot.png",width=10,height=18,units="cm")
# 
# png("Figures/Nm_tot_ggplot.png", width = 1000, height = 1800,units = "px")
# gg_all
# dev.off()




  traceplot(as.matrix(jags_concentration_factor_other_R_unif_Nmtot_simul[[1]]$BUGSoutput$sims.list$eta_alphaN))
     
  
  plot(density(jags_concentration_factor_other_R_unif_Nmtot_simul[[1]]$BUGSoutput$sims.list$eta_alphaN))



mcmcplot(jags_indpt_hier_split_simul[[1]],parms = "split")
mcmcplot(jags_indpt_hier_split_simul[[1]],parms = "mu_split")
