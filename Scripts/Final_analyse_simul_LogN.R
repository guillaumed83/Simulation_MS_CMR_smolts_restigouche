
setwd("C:/Users/DauphinGU/Desktop/smolts_restigouche/Model_CMR/2023-07-30_jfb_model/Simulation")



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


library(grid) ## To combine base plot and ggplot 

theme_set(theme_minimal())

#load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_split.RData")


load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2023-07-30_jfb_model/Simulation/outputs/_simulation_coda_R2jags_final_update_2023-07-30_M1-M2_LogN.RData")
load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2023-07-30_jfb_model/Simulation/outputs/_simulation_coda_R2jags_final_update_2023-07-30_M3-M4_LogN.RData")

#load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2020-07-07_clean_scripts_data_2019/Output/_simulation_coda_R2jags_final_update_2022-12-14_M5_LogN.RData")
load("C:/Users/DauphinGu/Desktop/smolts_restigouche/Model_CMR/2023-07-30_jfb_model/Simulation/outputs/_simulation_coda_R2jags_final_update_2023-07-29_M5b_LogN.RData")



## checking some dirichhlet stuff

# K <- runif(10000,0.01,10000)
# 
# temp1 <- K %*% t(p)
# par(mfrow=c(2,2))
# 
# for (i in 1:4){
#   plot(density(temp1[,i]))
# }
# 
# temp<- rdirichlet(10000,K %*% t(p) )
# par(mfrow=c(2,2))
# 
# for (i in 1:4){
#   plot(density(temp[,i]))
# }
# 
# apply(temp,2, mean)
# 
# par(mfrow=c(2,4))
# temp<- rdirichlet(10000,0.01*p)
# 
# for (i in 1:4){
#   plot(density(temp[,i],na.rm=T))
# }
# temp<- rdirichlet(10000,15*p)
# 
# for (i in 1:4){
#   plot(density(temp[,i]))
# }
# 
# temp<- rdirichlet(10000,0.1*p)
# par(mfrow=c(2,2))
# for (i in 1:4){
#   plot(density(temp[,i],na.rm=T))
# }
# for (i in 1:4){
# plot(1:10000,temp[,i])
# }
# 
# for (i in 1:4){
#   plot(1:10000,temp[,i])
# }
# 
# temp<- rdirichlet(10000,1*p)
# par(mfrow=c(2,2))
# for (i in 1:4){
#   plot(density(temp[,i]))
#   abline(v=p[i],col="red")
# }
# 
# temp<- rdirichlet(10000,2*p)
# par(mfrow=c(2,2))
# for (i in 1:4){
#   plot(density(temp[,i]))
#   abline(v=p[i],col="red")
# }
# 
# temp<- rdirichlet(10000,5*p)
# par(mfrow=c(2,2))
# for (i in 1:4){
#   plot(density(temp[,i]))
#   abline(v=p[i],col="red")
# }
# temp<- rdirichlet(10000,10*p)
# par(mfrow=c(2,2))
# for (i in 1:4){
#   plot(density(temp[,i]))
#   abline(v=p[i],col="red")
# }
# 
# temp<- rdirichlet(10000,15*p)
# par(mfrow=c(2,2))
# for (i in 1:4){
#   plot(density(temp[,i]))
#   abline(v=p[i],col="red")
# }
# for (i in 1:4){
#   plot(1:10000,temp[,i])
# }
# 
# temp<- rdirichlet(10000,50*p)
# par(mfrow=c(2,2))
# for (i in 1:4){
#   plot(density(temp[,i]))
#   abline(v=p[i],col="red")
# }


waic_1 <- as.data.frame(list.rbind(waic_indpt_merge_no_hier )) %>% mutate(ID=1,replicate=rownames(.))
waic_2 <- as.data.frame(list.rbind(waic_indpt_no_hier)) %>% mutate(ID=2,replicate=rownames(.))
waic_3 <- as.data.frame(list.rbind(waic_indpt_hier)) %>% mutate(ID=3,replicate=rownames(.))
waic_4 <- as.data.frame(list.rbind(waic_indpt_hier_split)) %>% mutate(ID=4,replicate=rownames(.))
#waic_5 <- as.data.frame(list.rbind(waic_DM)) %>% mutate(ID=5,replicate=rownames(.))
waic_5 <- as.data.frame(list.rbind(waic_DMb)) %>% mutate(ID=5,replicate=rownames(.))

all_waic <- rbind(waic_1,waic_2,waic_3, waic_4,waic_5)
all_waic %>% group_by(ID) %>%summarise(MEAN=mean(waic),SD=sd(waic),MEAN_p = mean(p_waic),SD_p = sd(p_waic))  

par(mfrow=c(1,1))
pretty_y <- pretty(all_waic$waic)

plot(0,0,type="n",
     xlab="Replicates",ylab="waic",
     xlim=c(0,10),ylim=c(min(pretty_y),max(pretty_y)))#c(600,800)) #c(min(pretty_y),max(pretty_y)))#

col_ <- c("darkorange4",
          "goldenrod1","goldenrod2",
          "firebrick3",
          "skyblue1" ) #"goldenrod3",


for (i in 1:length(col_)){
  temp <- all_waic %>% filter(ID == i)
  points(temp$replicate,sort(temp$waic),type="o",pch=16,col=col_[i])
  
  abline(h=mean(temp$waic),lty=2,lwd=2,col=col_[i])
  
}
legend("topleft",legend=c("indpt_merge_no_hier",
                          "indpt_no_hier",
                          "indpt_hier","indpt_hier_split",
                          "DM"),
       pch=rep(16,11),col=col_)


waic_sum <- data.frame(mean_waic=round(c(mean(waic_1[,1]),
                                    mean(waic_2[,1]),
                                    mean(waic_3[,1]),
                                    mean(waic_4[,1]),
                                    mean(waic_5[,1]))),
                       sd_waic=round(c(sd(waic_1[,1]),
                                 sd(waic_2[,1]),
                                 sd(waic_3[,1]),
                                 sd(waic_4[,1]),
                                 sd(waic_5[,1])),2))

waic_sum$mean_waic <- c(

mean(waic_1[,1]),


mean(waic_2[,1]),
mean(waic_3[,1]),
mean(waic_4[,1]),
mean(waic_5[,1]))


round(waic_sum$mean_waic)

source("simul_helper.R")


## We drop model 8b (merge of the 2 downstream RSTs using only the downstream recaptures)
## 


### Keeping the simulations we are interested in
simul_names <- c("jags_indpt_merge_no_hier_simul",
                 "jags_indpt_no_hier_simul",
                 "jags_indpt_hier_simul",
                 "jags_indpt_hier_split_simul",
                 "jags_DMb_simul"
             
                 )#<-  ls()[grepl("jags_", ls())] 

temp_names <-character(0)
for(i in 1:15){
  temp_names[i]<-paste("Nm_tot[",i,"]",sep="")
}

split_names <-character(0)
for(i in 1:15){
  split_names[i]<-paste("split[",i,"]",sep="")
}

theta4_names <-character(0)
for(i in 1:15){
 theta4_names[i]<-paste("theta[",i,",4]",sep="")
}
theta5_names <-character(0)
for(i in 1:15){
  theta5_names[i]<-paste("theta[",i,",5]",sep="")
}


par(mfrow=c(3,1))

for(i in 1:10){
  temp_split_sum <- as.data.frame(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary[row.names(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary) %in% split_names,] )
  temp_theta4_sum <- as.data.frame(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary[row.names(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary) %in% theta4_names,] )
  temp_theta5_sum <- as.data.frame(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary[row.names(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary) %in% theta5_names,] )
  
  if(i==1){
    plot(temp_split_sum$`50%`[y_wheels[,4]==1],temp_theta4_sum$`50%`,
         xlim=c(0,1),ylim=c(0,1),xlab="split",ylab="theta 45",
         pch=16,col=rgb(75,75,75,alpha=100,maxColorValue = 255))
  }else{
    points(temp_split_sum$`50%`[y_wheels[,4]==1],temp_theta4_sum$`50%`,
           pch=16,col=rgb(75,75,75,alpha=100,maxColorValue = 255))
  }
}

for(i in 1:10){
  temp_split_sum <- as.data.frame(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary[row.names(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary) %in% split_names,] )
  temp_theta4_sum <- as.data.frame(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary[row.names(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary) %in% theta4_names,] )
  temp_theta5_sum <- as.data.frame(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary[row.names(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary) %in% theta5_names,] )
  
  if(i==1){
    plot(temp_split_sum$`50%`[y_wheels[,5]==1],temp_theta5_sum$`50%`,
         xlim=c(0,1),ylim=c(0,1),xlab="split",ylab="theta 5",
         pch=16,col=rgb(75,75,75,alpha=100,maxColorValue = 255))
  }else{
    points(temp_split_sum$`50%`[y_wheels[,5]==1],temp_theta5_sum$`50%`,
           pch=16,col=rgb(75,75,75,alpha=100,maxColorValue = 255))
  }
}

for(i in 1:10){
  temp_split_sum <- as.data.frame(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary[row.names(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary) %in% split_names,] )
  temp_theta4_sum <- as.data.frame(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary[row.names(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary) %in% theta4_names,] )
  temp_theta5_sum <- as.data.frame(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary[row.names(jags_indpt_hier_split_simul[[i]]$BUGSoutput$summary) %in% theta5_names,] )
  
  if(i==1){
    plot(temp_theta4_sum$`50%`[c(2,3,4,7,8,9,10,12,13,14)],temp_theta5_sum$`50%`[c(1,3,4,5,6,7,8,9,10,11)],
         xlim=c(0,1),ylim=c(0,1),xlab="theta 4",ylab="theta 5",
         pch=16,col=rgb(75,75,75,alpha=100,maxColorValue = 255))
  }else{
    points(temp_theta4_sum$`50%`[c(2,3,4,7,8,9,10,12,13,14)],temp_theta5_sum$`50%`[c(1,3,4,5,6,7,8,9,10,11)],
           pch=16,col=rgb(75,75,75,alpha=100,maxColorValue = 255))
  }
}

png("Figures/Simulation_total_abundance.png",height=1300,width=1800,pointsize = 32)
#par(mfrow=c(2,4),mar=c(2,2,3,2),oma=c(2,4,2,1))
par(mfrow=c(2,3),mar=c(2,2,3,2),oma=c(2,4,2,1))

ymax <- 3000000

## M1
## all wheels independant no hierarchy 
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,ymax),pch=16,type="n",main="M1",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,ymax,500000),labels=seq(0,ymax/1000,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_merge_no_hier_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M2
## all wheels independant annual hierarchy
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,ymax),pch=16,type="n",main="M2",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,ymax,500000),labels=seq(0,ymax/1000,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_no_hier_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M3
## all wheels independant annual hierarchy + split
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,ymax),pch=16,type="n",main="M3",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,ymax,500000),labels=seq(0,ymax/1000,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_hier_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

## M4
## all wheels independant annual hierarchy + pooled
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,ymax),pch=16,type="n",main="M4",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,ymax,500000),labels=seq(0,ymax/1000,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_indpt_hier_split_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)

# ## M5
# ## all wheels independant annual hierarchy + split +pooled
# plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,ymax),pch=16,type="n",main="M5",
#      xlab="",ylab="",axes=F)
# axis(2,at=seq(0,ymax,500000),labels=seq(0,ymax/1000,500),las=2)
# axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
# abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
# box()
# 
# for(k in 1:3){
#   for (i in 1:10){
#     plot_sim(jags_DM_simul,i,k)
#   }
# }
# points(1:15,Tot_N,col="red",type="b",pch=16)



## M5
## all wheels independant annual hierarchy + split +pooled
plot(1:15,Tot_N,xlim=c(1,15),ylim=c(50000,ymax),pch=16,type="n",main="M5b",
     xlab="",ylab="",axes=F)
axis(2,at=seq(0,ymax,500000),labels=seq(0,ymax/1000,500),las=2)
axis(1,at=seq(1,15,1),labels=c(1,NA,3,NA,5,NA,7,NA,9,NA,11,NA,13,NA,15))
abline(h=seq(500000,3500000,500000),lty=2,col="grey85")
box()

for(k in 1:3){
  for (i in 1:10){
    plot_sim(jags_DMb_simul,i,k)
  }
}
points(1:15,Tot_N,col="red",type="b",pch=16)



mtext( "Total abundance (x1000)",2,out=T,line=2)
mtext( "Years",1,out=T,line=1)


#points(1:15,ifelse(y_wheels[,4]==0|y_wheels[,5]==0,100000,NA),pch="*",cex=1.2)
dev.off()


### CVs
## We drop model 8b (merge of the 2 downstream RSTs using only the downstream recaptures)
## 

CV_Nm_tot <- list()

all_simul <- list(jags_indpt_merge_no_hier_simul,
                  jags_indpt_no_hier_simul,
                  jags_indpt_hier_simul,
                  jags_indpt_hier_split_simul,
                  jags_DMb_simul
                  )

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

for (n in 1:length(simul_names)){
  print(mean(CV_Nm_tot[[n]]) )
}



png("Figures/Simulation_CV_total_abundance.png",height=1800,width=1300,pointsize = 32)
par(mfrow=c(3,2),mar=c(2,2,3,2),oma=c(2,4,2,1))
for (n in 1:length(simul_names)){
  plot(rep(1:15,10),CV_Nm_tot[[n]],col=rgb(30,130,181,alpha=75, maxColorValue = 255),ylim=c(0,1),pch=16,main=paste("M",n,sep=""),yaxt="none",type="n")
  axis(2,at=seq(0,4,0.2),labels=seq(0,4,0.2),las=2)
  abline(h=seq(0,1,0.2),lty=2,col="grey85")
  points(rep(1:15,10),CV_Nm_tot[[n]],col=rgb(30,130,181,alpha=75, maxColorValue = 255),ylim=c(0,1),pch=16)
  abline(h=apply(CV_Nm_tot[[n]],2,mean),col=rgb(30,130,181,alpha=170, maxColorValue = 255),lty=2)
  abline(h=mean(CV_Nm_tot[[n]]),col="red",lty=2,lwd=2)
  
  text(12.5,1,paste("Avg. CV = ",round(mean(CV_Nm_tot[[n]]),3)))
  
}  
mtext( "Coefficient of Variation (std. dev. / mean)",2,out=T,line=2)

dev.off()



####  difference obs - pred - Nm_tot ####
png("Figures/Simulation_diff_total_abundance.png",height=1200,width=1800,pointsize = 32)
list_diff_Tot_N <- list()
par(mfrow=c(2,3),mar=c(2,2,3,2),oma=c(2,4,2,1))
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





#### Fig3 old ####
list_diff_tot_N_plot <- list()
gg_list <- list()

for (n in 1:length(simul_names)){
  diff_obs_pred2(simul_names[n],"Tot_N")
  list_diff_tot_N_plot[[n]] <- temp_diff_par
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
    geom_vline(xintercept = 0,col="black",linetype="dashed",size=1.5)+
    geom_vline(xintercept = mean(dplyr::bind_rows(temp_list)$value),col="red",linetype="dashed",size=1.5  ) + 
    annotate("text",x=-Inf,y=14.5, label = paste("Avg. Diff. = ",round(mean(dplyr::bind_rows(temp_list)$value))),  hjust = -0.1) +
    xlab("diff. Total abundance") +
    #xlim(0,2.5)+
    theme(legend.position="none")
  
  gg_list[[i]] <- gg1
}

gg_list[[1]]

####

tiff("Figures/diff_Nm_tot_ggplot.tiff", units="in", width=20, height=15, res=250)

ggarrange(gg_list[[1]],
          gg_list[[2]],
          gg_list[[3]],
          gg_list[[4]],
          gg_list[[5]],
          
          labels = c("M1", "M2","M3","M4","M5"),
          ncol = 2, nrow = 3) +
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm")) 

dev.off()

##### calculation (True value - estimate)/ True value ####

#diff_obs_pred_pct

#### Fig3 new ####
list_diff_tot_N__pct_plot <- list()
gg_list <- list()

for (n in 1:length(simul_names)){
  diff_obs_pred_pct(simul_names[n],"Tot_N")
  list_diff_tot_N__pct_plot[[n]] <- temp_diff_par
}  

for( i in 1:length(all_simul)){
  
  temp_list <- list()
  
  for(j in 1:n_rep){
    temp <- as.data.frame(list_diff_tot_N__pct_plot[[i]][[j]])
    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    temp_list[[j]] <- temp
  }

  gg1 <- ggplot(temp_list[[1]], aes(x=value,y=Year)) +
    scale_x_continuous(limits=c(-420,100))+
    #geom_vline(xintercept = mean(temp_list[[1]]$value),col="#0072B250",linetype="dashed",size=1,alpha=0.2)+
    geom_density_ridges(scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  
  for(k in 2:10){
    gg1 <- gg1 + 
      #geom_vline(xintercept = mean(temp_list[[k]]$value),col="#0072B250",linetype="dashed",size=1,alpha=0.2)+
      geom_density_ridges(data=temp_list[[k]], aes(x=value,y=Year),scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  }
  
  gg1 <- gg1 + 
    #geom_point(data=dat_tot_N,mapping=aes(x=0, y= Year),col="red",size=1.5) + 
    geom_vline(xintercept = 0,col="black",linetype="dashed",size=1.5)+
    geom_vline(xintercept = mean(dplyr::bind_rows(temp_list)$value),col="red",linetype="dashed",size=1.5  ) + 
    annotate("text",x=-Inf,y=14.5, label = paste("Avg. Diff. = ",round(mean(dplyr::bind_rows(temp_list)$value),1)),  hjust = -0.1) +
    xlab("(True - Est.) / True") +#xlab("diff. Total abundance") +
    #xlim(0,2.5)+
    theme(legend.position="none",
          axis.title=element_text(size=16))
  
  gg_list[[i]] <- gg1
}

gg_list[[1]]

####

tiff("Figures/diff_Nm_tot_pcent_ggplot.tiff", units="in", width=20, height=15, res=250)

ggarrange(gg_list[[1]],
          gg_list[[2]],
          gg_list[[3]],
          gg_list[[4]],
          gg_list[[5]],
          
          labels = c("M1", "M2","M3","M4","M5"),
          ncol = 2, nrow = 3) +
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm")) 

dev.off()

#####

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

###################
### first wheel ###
###################


temp_names <-character(0)
for(i in 1:15){
  temp_names[i]<-paste("Nm[",i,",1]",sep="")
}

par(mfrow=c(2,3))

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
  diff_obs_pred2(simul_names[n],"N_1")
  list_diff_N1_plot[[n]] <- temp_diff_par
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

tiff("Figures/diff_N1_ggplot.tiff", units="in", width=20, height=15, res=250)

ggarrange(gg_N1_list[[1]],
          gg_N1_list[[2]],
          gg_N1_list[[3]],
          gg_N1_list[[4]],
          gg_N1_list[[5]],
         
          labels = c("M1", "M2","M3","M4","M5"),
          ncol = 2, nrow = 3) +
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm")) 

dev.off()


#######

list_diff_N1_plot <- list()
gg_N1_list <- list()

for (n in 1:length(simul_names)){
  diff_obs_pred_pct(simul_names[n],"N_1")
  list_diff_N1_plot[[n]] <- temp_diff_par
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
    scale_x_continuous(limits=c(-200,100))+
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

tiff("Figures/diff_N1_ggplot_pct.tiff", units="in", width=20, height=15, res=250)

ggarrange(gg_N1_list[[1]],
          gg_N1_list[[2]],
          gg_N1_list[[3]],
          gg_N1_list[[4]],
          gg_N1_list[[5]],
          
          labels = c("M1", "M2","M3","M4","M5"),
          ncol = 2, nrow = 3) +
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm")) 

dev.off()


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
  
  print(plot_p_smolt(all_simul[[n]],n,FALSE))
  
  dev.off()
}
for( n in 1:length(all_simul)){
  # png(paste("Figures/proportions_",simul_names[i],".png",sep=""),
  #     height=1200,width=1200,pointsize =16)
  # 
  
  tiff(paste("Figures/proportions_all_",simul_names[n],".tiff",sep=""),
       units="in", width=7, height=7, res=300)
  
  print(plot_p_smolt(all_simul[[n]],n,TRUE))
  
  dev.off()
}




# 
# 
# ## Return proportions of smolt abundances
# n=1 ## model indpt no hier
# 
# return_p_smolt(all_simul[[n]],n)
# 
# 
# p_all_indpt_no_hier_merge <- list_pall_simul
# 
# mean_p_all <- function(x){
#   temp <- as.data.frame(x) %>% mutate(p_greater_1 = ifelse(value>1,1,0)) %>%
#     group_by(Year) %>% summarise(mean_p_greater_1_pcent = mean(p_greater_1)*100,
#                                  mean_p=mean(value),
#                                  q_75=quantile(value,probs=0.75),
#                                  q_90=quantile(value,probs=0.90),
#                                  q_95=quantile(value,probs=0.95)
#     )
#   
# }
# 
# sum_p_all_indpt_no_hier_merge <- p_all_indpt_no_hier_merge |> lmap(mean_p_all)
# mean_p_greater_indpt_no_hier_merge_all_reps <- rbindlist(sum_p_all_indpt_no_hier_merge) %>% group_by(Year) %>% summarise(mean_p_greater1 = mean(mean_p_greater_1_pcent))
# 
# ## Years with 3 upstream wheel
# years_n_rst <- rowSums(y_wheels[,1:3])
# col_rst <- c("#e66101","#fdb863","#5e3c99")
# 
# ## Return proportions of smolt abundances
# return_p_smolt(all_simul[[n]],n)
# 
# 
# tiff("Figures/2022-12-15/LogN/probability_proportions_greater_than_1_indpt_no_hier_merge.tiff", units="in", width=12, height=8, res=250)
# 
# par(mfrow=c(1,2))
# 
# 
# for(i in 1:n_rep){
#   
#   if(i==1){
#     plot(1:15,sum_p_all_indpt_no_hier_merge[[i]]$mean_p_greater_1_pcent,pch=16, col=col_rst[years_n_rst],ylim=c(0,100),
#          xlab="Years",
#          ylab="Probability that upstream proportions > 1 (%)",
#          main=simul_names[n])
#   }else{
#     points(1:15,sum_p_all_indpt_no_hier_merge[[i]]$mean_p_greater_1_pcent,pch=16, col=col_rst[years_n_rst] )
#   }
#   abline(h=5,lty=2)
#   legend("topleft",legend = c("Number of RSTs active","1","2","3"), pch=c(NA,16,16,16),col=c(NA,col_rst), )
# }
# 
# 
# ## Plotting the distribution of the proportion of all upstream RSTS
# gg1 <-  ggplot(p_all_indpt_no_hier_merge[[1]], aes(x=value,y=Year)) +
#   geom_density_ridges(rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+#,aes(height=1))+#..ndensity..))+#+
#   geom_point(data=sum_p_all_indpt_no_hier_merge[[1]],mapping=aes(x=mean_p,y=Year),col=col_rst[years_n_rst])+
#   xlim(0,3)+
#   annotate("rect",xmin=1,xmax=3,ymin=0,ymax=17,fill="red",alpha=0.15) +
#   geom_vline(aes(xintercept = 1),col="red",linetype="dashed") +
#   #scale_x_continuous(breaks=c(0,0.25,0.50,0.75,1,1.25),labels=c(0,0.25,0.50,0.75,1,1.25))+
#   labs(x="sum of upstream RSTs proportions")
# #labels = unit_format(unit = "M", scale = 1e-6))
# 
# 
# 
# for(k in 2:10){
#   gg1 <- gg1 + geom_density_ridges(data=p_all_indpt_no_hier_merge[[k]], aes(x=value,y=Year),rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+
#     geom_point(data=sum_p_all_indpt_no_hier_merge[[k]],mapping=aes(x=mean_p,y=Year),col=col_rst[years_n_rst])
# }
# gg_no_hier_merge <- gg1 +  annotate("text", x =2.8, y =(1:16)+0.2, label = c(round(mean_p_greater_indpt_no_hier_merge_all_reps$mean_p_greater1,1),"Avg. p > 1"))
# 
# vp.Right <- viewport(height=unit(0.85, "npc"), width=unit(0.5, "npc"), 
#                      just=c("left","top"), 
#                      y=0.92, x=0.5)
# 
# 
# print(gg_no_hier_merge,vp=vp.Right)
# mtext(3,simul_names[n])
# 
# dev.off()
# 
# 
# 
# ## Return proportions of smolt abundances
# n=2 ## model indpt no hier
# return_p_smolt(all_simul[[n]],n)
# p_all_indpt_no_hier <- list_pall_simul
# sum_p_all_indpt_no_hier <- p_all_indpt_no_hier |> lmap(mean_p_all)
# mean_p_greater_indpt_no_hier_all_reps <- rbindlist(sum_p_all_indpt_no_hier) %>% group_by(Year) %>% summarise(mean_p_greater1 = mean(mean_p_greater_1_pcent))
# 
# ## Years with 3 upstream wheel
# years_n_rst <- rowSums(y_wheels[,1:3])
# col_rst <- c("#e66101","#fdb863","#5e3c99")
# 
# tiff("Figures/2022-12-15/LogN/probability_proportions_greater_than_1_indpt_no_hier.tiff", units="in", width=12, height=8, res=250)
# 
# par(mfrow=c(1,2))
# for(i in 1:n_rep){
#   if(i==1){
#     plot(1:15,sum_p_all_indpt_no_hier[[i]]$mean_p_greater_1_pcent,pch=16, col=col_rst[years_n_rst],ylim=c(0,100),
#          xlab="Years",
#          ylab="Probability that upstream proportions > 1 (%)",
#          main=simul_names[n])
#   }else{
#     points(1:15,sum_p_all_indpt_no_hier[[i]]$mean_p_greater_1_pcent,pch=16, col=col_rst[years_n_rst] )
#   }
#   abline(h=5,lty=2)
#   legend("topleft",legend = c("Number of RSTs active","1","2","3"), pch=c(NA,16,16,16),col=c(NA,col_rst), )
# }
# 
# 
# ## Plotting the distribution of the proportion of all upstream RSTS
# gg1 <-  ggplot(p_all_indpt_no_hier[[1]], aes(x=value,y=Year)) +
#   geom_density_ridges(rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+#,aes(height=1))+#..ndensity..))+#+
#   geom_point(data=sum_p_all_indpt_no_hier[[1]],mapping=aes(x=mean_p,y=Year),col=col_rst[years_n_rst])+
#   xlim(0,3)+
#   annotate("rect",xmin=1,xmax=3,ymin=0,ymax=17,fill="red",alpha=0.15) +
#   geom_vline(aes(xintercept = 1),col="red",linetype="dashed") +
#   #scale_x_continuous(breaks=c(0,0.25,0.50,0.75,1,1.25),labels=c(0,0.25,0.50,0.75,1,1.25))+
#   labs(x="sum of upstream RSTs proportions")
# #labels = unit_format(unit = "M", scale = 1e-6))
# 
# for(k in 2:10){
#   gg1 <- gg1 + geom_density_ridges(data=p_all_indpt_no_hier[[k]], aes(x=value,y=Year),rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+
#     geom_point(data=sum_p_all_indpt_no_hier[[k]],mapping=aes(x=mean_p,y=Year),col=col_rst[years_n_rst])
# }
# 
# gg_no_hier <- gg1 +  annotate("text", x =2.8, y =(1:16)+0.2, label = c(round(mean_p_greater_indpt_no_hier_all_reps$mean_p_greater1,1),"Avg. p > 1"))
# 
# 
# vp.Right <- viewport(height=unit(0.85, "npc"), width=unit(0.5, "npc"), 
#                            just=c("left","top"), 
#                            y=0.92, x=0.5)
# 
# print(gg_no_hier,vp=vp.Right)
# mtext(3,simul_names[n])
# 
# dev.off()
# 
# 
# 
# 
# n=5 ## model indpt with annual hier
# 
# return_p_smolt(all_simul[[n]],n)
# 
# 
# temp_<-as.data.frame(jags_DMb_simul[[2]]$BUGSoutput$sims.list$Nm[,,3] / jags_DMb_simul[[2]]$BUGSoutput$sims.list$Nm_tot)
# 
# 
# jags_DMb_simul[[2]]$BUGSoutput$sims.list$Nm_tot - 
#   (jags_DMb_simul[[2]]$BUGSoutput$sims.list$Nm[,,1]+
#   jags_DMb_simul[[2]]$BUGSoutput$sims.list$Nm[,,2]+
#   jags_DMb_simul[[2]]$BUGSoutput$sims.list$Nm[,,3]+
#   jags_DMb_simul[[2]]$BUGSoutput$sims.list$Nm_rest)
# 
# 
# jags_DM_simul[[2]]$BUGSoutput$sims.list$p_smolt_prod[,,1]+
#   jags_DM_simul[[2]]$BUGSoutput$sims.list$p_smolt_prod[,,2]+
#   jags_DM_simul[[2]]$BUGSoutput$sims.list$p_smolt_prod[,,3]+
#   jags_DM_simul[[2]]$BUGSoutput$sims.list$p_smolt_prod[,,4]
# 


## Function that helps us generating summary stats on p>1 
mean_p_all <- function(x){
  temp <- as.data.frame(x) %>% mutate(p_greater_1 = ifelse(value>1,1,0)) %>%
    group_by(Year) %>% summarise(mean_p_greater_1_pcent = mean(p_greater_1)*100,
                                 mean_p=mean(value),
                                 q_75=quantile(value,probs=0.75),
                                 q_90=quantile(value,probs=0.90),
                                 q_95=quantile(value,probs=0.95)
    )
  
} 

for( n in 1:length(all_simul)){
  
  ## Return proportions of smolt abundances
  ## in object list_pall_simul
  return_p_smolt(all_simul[[n]],n)

  sum_p_all <- list_pall_simul |> lmap(mean_p_all)
  
  #sum_p_all <- list_p3_simul |> lmap(mean_p_all)
  mean_p_greater_all_reps <- rbindlist(sum_p_all) %>% group_by(Year) %>% summarise(mean_p_greater1 = mean(mean_p_greater_1_pcent))
  
  ## Years with 3 upstream wheel
  years_n_rst <- rowSums(y_wheels[,1:3])
  col_rst <- c("#e66101","#fdb863","#5e3c99")
  
  tiff(paste("Figures/probability_proportions_greater_than_1_",simul_names[n],".tiff",sep=""),
       units="in", width=10, height=8, res=300)
  
  par(mfrow=c(1,2))
  
  
  for(i in 1:n_rep){
    
    if(i==1){
      plot(1:15,sum_p_all[[i]]$mean_p_greater_1_pcent,pch=16, col=col_rst[years_n_rst],ylim=c(0,60),
           xlab="Years",
           ylab="Probability that upstream proportions > 1 (%)",
           main=simul_names[n])
    }else{
      points(1:15,sum_p_all[[i]]$mean_p_greater_1_pcent,pch=16, col=col_rst[years_n_rst] )
    }
    abline(h=5,lty=2)
    legend("topleft",legend = c("Number of RSTs active","1","2","3"), pch=c(NA,16,16,16),col=c(NA,col_rst), )
  }
  
  
  ## Plotting the distribution of the proportion of all upstream RSTS
  gg1 <-  ggplot(list_pall_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+#,aes(height=1))+#..ndensity..))+#+
    geom_point(data=sum_p_all[[1]],mapping=aes(x=mean_p,y=Year),col=col_rst[years_n_rst])+
    xlim(0,3)+
    annotate("rect",xmin=1,xmax=3,ymin=0,ymax=17,fill="red",alpha=0.15) +
    geom_vline(aes(xintercept = 1),col="red",linetype="dashed") +
    #scale_x_continuous(breaks=c(0,0.25,0.50,0.75,1,1.25),labels=c(0,0.25,0.50,0.75,1,1.25))+
    labs(x="sum of upstream RSTs proportions")
  #labels = unit_format(unit = "M", scale = 1e-6))
  
  
  
  for(k in 2:10){
    gg1 <- gg1 + geom_density_ridges(data=list_pall_simul[[k]], aes(x=value,y=Year),rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+
      geom_point(data=sum_p_all[[k]],mapping=aes(x=mean_p,y=Year),col=col_rst[years_n_rst])
  }
  
  gg_all <- gg1 +  annotate("text", x =2.8, y =(1:16)+0.2, label = c(round(mean_p_greater_all_reps$mean_p_greater1,1),"Avg. p > 1"))
  
  
  vp.Right <- viewport(height=unit(0.85, "npc"), width=unit(0.5, "npc"), 
                       just=c("left","top"), 
                       y=0.92, x=0.5)
 
  print(gg_all,vp=vp.Right)
  mtext(3,simul_names[n])

  dev.off()
}

  

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
          
           labels = c("M1", "M2","M3","M4","M5"),
          ncol = 2, nrow = 3) +
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
                    
                    labels = c("M1", "M2","M3","M4","M5"),
                    ncol = 2, nrow = 3) +
  theme(plot.margin = margin(0.3,0.3,0.3,0.3, "cm")) 



tiff("Figures/N1_ggplot.tiff", units="in", width=20, height=15, res=250)

gg_all

dev.off()






# ggsave("Figures/Nm_tot_ggplot.png",width=10,height=18,units="cm")
# 
# png("Figures/Nm_tot_ggplot.png", width = 1000, height = 1800,units = "px")
# gg_all
# dev.off()
sum_fun <- function(x){
  c(mean(x),sd(x),  quantile(x,probs = c(0.025,0.25,0.5,0.75,0.975)))
}

par(mfrow=c(3,4))
sum_eta <- array(0,dim=c(10,7))


for(j in 1:10){
  
  plot(density(jags_DMb_simul[[j]]$BUGSoutput$sims.list$eta_alphaN),main=paste0("replicate",j))
  abline(v=15,col="red")
  abline(v= median(jags_DMb_simul[[j]]$BUGSoutput$sims.list$eta_alphaN),col="black")
  sum_eta[j,]<-sum_fun(jags_DMb_simul[[j]]$BUGSoutput$sims.list$eta_alphaN)
}


mcmcplot(jags_indpt_hier_split_simul[[1]],parms = "split")
mcmcplot(jags_indpt_hier_split_simul[[1]],parms = "mu_split")
