
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


diff_obs_pred <- function(name_simul,par_temp){
  
  temp_simul <-  eval(as.name(name_simul))

  temp_diff_sum <- list()
  if(par_temp=="Tot_N"){
    temp_diff_totN <- list()
  }
  
  for (n in 1:n_rep){
    ## difference observed - predicted 
    ## create one case for each parameter (i'm lazy) 
    
    if(par_temp=="Tot_N"){
      temp_diff <- -sweep(as.matrix(temp_simul[[n]]$BUGSoutput$sims.list$Nm_tot), 2, Tot_N )
    }
    if(par_temp=="N_1"){
      temp_diff <- -sweep(as.matrix(temp_simul[[n]]$BUGSoutput$sims.list$Nm[,,1]), 2, Tot_w[,1] )
    }

    temp_diff_sum[[n]] <- as.data.frame(
                                        t(rbind(
      
                                        apply(temp_diff,2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975)) ,
                                        mean=apply(temp_diff,2,mean),
                                        sd=apply(temp_diff,2,sd)
    ))
    )
  }
  temp_diff_sum <<- temp_diff_sum 
}


###_______________________
diff_obs_pred2 <- function(name_simul,par_temp){
  
  temp_simul <-  eval(as.name(name_simul))
  temp_diff_par <- list()
 
  for (n in 1:n_rep){
    ## difference observed - predicted 
    ## create one case for each parameter (i'm lazy) 
    
    if(par_temp=="Tot_N"){
      temp_diff <- -sweep(as.matrix(temp_simul[[n]]$BUGSoutput$sims.list$Nm_tot), 2, Tot_N )
      temp_diff_par[[n]] <- temp_diff 
    }
    if(par_temp=="N_1"){
      temp_diff <- -sweep(as.matrix(temp_simul[[n]]$BUGSoutput$sims.list$Nm[,,1]), 2, Tot_w[,1] )
      temp_diff_par[[n]] <- temp_diff 
    }
  }
 
  if(par_temp=="Tot_N"){
    temp_diff_par <<-  temp_diff_par
  }
  if(par_temp=="N_1"){
    temp_diff_par <<-  temp_diff_par
  }
}


### plotting differences between true/observed values and predicted/estimated values from a simulation object
### 
diff_obs_pred_pct <- function(name_simul,par_temp){
  
  temp_simul <-  eval(as.name(name_simul))
  temp_diff_par <- list()

  for (n in 1:n_rep){
    ## difference observed - predicted 
    ## create one case for each parameter (i'm lazy) 
    
    if(par_temp=="Tot_N"){
      temp_diff <- -sweep(as.matrix(temp_simul[[n]]$BUGSoutput$sims.list$Nm_tot), 2, Tot_N )
      temp_diff <- sweep(temp_diff ,2,Tot_N,"/")
      temp_diff_par[[n]] <- temp_diff *100
    }
    if(par_temp=="N_1"){
      temp_diff <- -sweep(as.matrix(temp_simul[[n]]$BUGSoutput$sims.list$Nm[,,1]), 2, Tot_w[,1] )
      temp_diff <- sweep(temp_diff ,2,Tot_w[,1],"/")
      temp_diff_par[[n]] <- temp_diff *100
    }
  }
  if(par_temp=="Tot_N"){
    temp_diff_par <<-  temp_diff_par
  }
  if(par_temp=="N_1"){
    temp_diff_par <<-  temp_diff_par
  }
}

## Plotting total abundance from simulation output 
## colored polygons to represent cumulative uncertainty from each replicate

plot_sim <- function(name,i,k){

  if(grepl("all_simul",deparse(substitute(name)) ) | grepl("jags",deparse(substitute(name)) ) ){
    temp_df <- name[[i]]$BUGSoutput$summary[row.names(name[[i]]$BUGSoutput$summary) %in% temp_names,] %>% 
      as.data.frame() %>% select("2.5%","25%","50%","75%","97.5%")
  }else{
    temp_df <- name[[i]]
  }

  if(k==1){
    polygon(x=c(1:15,15:1),y=c(temp_df[,1],rev(temp_df[,2])),col=rgb(216,179,101,alpha=35,maxColorValue = 255),border=NA)
    polygon(x=c(1:15,15:1),y=c(temp_df[,4],rev(temp_df[,5])),col=rgb(216,179,101,alpha=35,maxColorValue = 255),border=NA)
  }
  
  if(k==2){
    polygon(x=c(1:15,15:1),y=c(temp_df[,2],rev(temp_df[,4])),col=rgb(90,180,172,alpha=35,maxColorValue = 255),border=NA)
  }

  if(k==3){
    points(1:15,temp_df[,3],pch=16,type="b",col= rgb(50,50,50,alpha=125,maxColorValue = 255))
  }
}


## Return proportions of abundance associated with each upstream RST (lists of data frames)
## Takes simulation MCMC object and position in the vector of the name of simulation
## (the latter attributes could be obtained within the function)
return_p_smolt <- function(name_simul,j){
  
  list_p1_simul <- list()
  list_p2_simul <- list()
  list_p3_simul <- list()
  list_pall_simul <- list()
  
  for(i in 1:n_rep){
 
      ## Proportion wheel 1
      temp1 <- as.data.frame(name_simul[[i]]$BUGSoutput$sims.list$Nm[,,1] / name_simul[[i]]$BUGSoutput$sims.list$Nm_tot)
      
      ## Proportion wheel 2    
      temp2 <- as.data.frame(name_simul[[i]]$BUGSoutput$sims.list$Nm[,,2] / name_simul[[i]]$BUGSoutput$sims.list$Nm_tot)
      years2keep <- y_wheels[,2]*1:15
      years2keep <- years2keep[years2keep!=0]
      yy <- 1:15
      yy <- yy[!yy %in% years2keep]
      
      ## year when RST not active, filling with 0s
      for(k in 1:length(yy)){ temp2[,yy[k]] <- 0 }
      
      ## Proportion wheel 3
      temp3 <- as.data.frame(name_simul[[i]]$BUGSoutput$sims.list$Nm[,,3] / name_simul[[i]]$BUGSoutput$sims.list$Nm_tot)
      years2keep <- y_wheels[,3]*1:15
      years2keep <- years2keep[years2keep!=0]
      yy <- 1:15
      yy <- yy[!yy %in% years2keep]
      
      ## year when RST not active, filling with 0s
      for(k in 1:length(yy)){ temp3[,yy[k]] <- 0 }
      
      names(temp1)<-seq(1,15,1)
      names(temp2)<-seq(1,15,1)
      names(temp3)<-seq(1,15,1)
      
      temp1 <- temp1 %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
      temp2 <- temp2 %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
      temp3 <- temp3 %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
      
      temp1$Year<-factor(temp1$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
      temp2$Year<-factor(temp2$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
      temp3$Year<-factor(temp3$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))

      temp_all <- data.frame(Year=temp1$Year, 
                             value=temp1$value+temp2$value+temp3$value,
                             anchor=temp1$anchor)
      
      list_p1_simul[[i]] <- temp1
      list_p2_simul[[i]] <- temp2
      list_p3_simul[[i]] <- temp3
      list_pall_simul[[i]] <- temp_all
   # } ## end of conditional loop (independent)
  } ## end of replicate loop
  
  list_p1_simul <<- list_p1_simul
  list_p2_simul <<- list_p2_simul
  list_p3_simul <<- list_p3_simul
  list_pall_simul <<- list_pall_simul
  
} ## end of return p smolt function


## Generates plots of proportions of abundance associated with each upstream RST (lists of data frames)
## Takes simulation MCMC object and position in the vector of the name of simulation
## (the latter attributes could be obtained within the function)
plot_p_smolt <- function(name_simul,j,full_plot){
  
    list_p_simul <- list()
 
  for(i in 1:n_rep){
    #if( grepl("indpt",deparse(substitute(name_simul)) ) ){
    if( grepl("indpt",simul_names[j] ) ){  
      temp <- as.data.frame(name_simul[[i]]$BUGSoutput$sims.list$Nm[,,1] / name_simul[[i]]$BUGSoutput$sims.list$Nm_tot)
    }else{
      temp <-as.data.frame( name_simul[[i]]$BUGSoutput$sims.list$p_smolt_prod[,,1] )
    }

    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    list_p_simul[[i]] <- temp
  }
 
  gg1 <- ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  
  for(i in 2:10){
    gg1 <- gg1 + geom_density_ridges(data=list_p_simul[[i]], aes(x=value,y=Year),scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  }
  
  if( grepl("indpt",simul_names[j] )==FALSE ){
    gg1 <- gg1 + geom_point(data=p_w,mapping=aes(x=p_1, y= Year),col="red",size=1.5) + 
      xlab("Proportion") +
      #xlim(0,2.5)+
      theme(legend.position="none",
            panel.border = element_rect(colour = "black", fill=NA))
  }else{
   
    gg1 <- gg1 + geom_point(data=p_w,mapping=aes(x=p_1, y= Year),col="red",size=1.5) + 
      xlim(0,3)+
      annotate("rect",xmin=1,xmax=3,ymin=0,ymax=17,fill="red",alpha=0.15) +
      geom_vline(aes(xintercept = 1),col="red",linetype="dashed") +
      xlab("Proportion") +
      theme(legend.position="none",
            panel.border = element_rect(colour = "black", fill=NA))
  }

  #______________
  
  for(i in 1:n_rep){
    if( grepl("indpt",simul_names[j] ) ){
      temp <- as.data.frame(name_simul[[i]]$BUGSoutput$sims.list$Nm[,,2] / name_simul[[i]]$BUGSoutput$sims.list$Nm_tot)
    }else{
      temp <-as.data.frame( name_simul[[i]]$BUGSoutput$sims.list$p_smolt_prod[,,2] )
    }
 
    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    
    years2keep <- y_wheels[,2]*1:15
    years2keep <- years2keep[years2keep!=0]
    
    temp <- temp %>% filter(Year %in% years2keep)
    
    yy <- 1:15
    yy <- yy[!yy %in% years2keep]
    
    for(k in 1:length(yy)){
      temp <-rbind(temp,c(yy[k],NA,NA))
     }
    list_p_simul[[i]] <- temp
  }

  gg2 <- ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  
  for(i in 2:10){
    gg2 <- gg2 + geom_density_ridges(data=list_p_simul[[i]], aes(x=value,y=Year),scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  }

  if( grepl("indpt",simul_names[j] )==FALSE ){
    gg2 <- gg2 + geom_point(data=p_w[years2keep,],mapping=aes(x=p_2, y= Year),col="red",size=1.5) +
      xlab("Proportion") +
      #xlim(0,5)+
      theme(legend.position="none",
            panel.border = element_rect(colour = "black", fill=NA))
  }else{
    gg2 <- gg2 + geom_point(data=p_w[years2keep,],mapping=aes(x=p_2, y= Year),col="red",size=1.5) +
      xlim(0,3)+
      annotate("rect",xmin=1,xmax=3,ymin=0,ymax=17,fill="red",alpha=0.15) +
      geom_vline(aes(xintercept = 1),col="red",linetype="dashed") +
      xlab("Proportion") +
      theme(legend.position="none",
            panel.border = element_rect(colour = "black", fill=NA))
  }

  #______________

  for(i in 1:n_rep){
    if( grepl("indpt",simul_names[j] ) ){
      temp <- as.data.frame(name_simul[[i]]$BUGSoutput$sims.list$Nm[,,3] / name_simul[[i]]$BUGSoutput$sims.list$Nm_tot)
    }else{
      temp <-as.data.frame( name_simul[[i]]$BUGSoutput$sims.list$p_smolt_prod[,,3] )
    }

    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    
    years2keep <- y_wheels[,3]*1:15
    years2keep <- years2keep[years2keep!=0]
    
    temp <- temp %>% filter(Year %in% years2keep)
    
    yy <- 1:15
    yy <- yy[!yy %in% years2keep]
    
    for(k in 1:length(yy)){
      temp <-rbind(temp,c(yy[k],NA,NA))
    }
    list_p_simul[[i]] <- temp
  }

  gg3 <- ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  
  for(i in 2:10){
    gg3 <- gg3 + geom_density_ridges(data=list_p_simul[[i]], aes(x=value,y=Year),scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  }
  
  temp_x_max <- max(ggplot_build(gg3)$layout$panel_params[[1]]$x.range )
  print(temp_x_max)
  if( temp_x_max <= 1.2 ){
    gg3 <- gg3 + geom_point(data=p_w[years2keep,],mapping=aes(x=p_3, y= Year),col="red",size=1.5) +
      xlab("Proportion") +
      #xlim(0,5)+
      theme(legend.position="none",
            panel.border = element_rect(colour = "black", fill=NA))
  }else{
    gg3 <- gg3 + geom_point(data=p_w[years2keep,],mapping=aes(x=p_2, y= Year),col="red",size=1.5) +
      xlim(0,3)+
      annotate("rect",xmin=1,xmax=3,ymin=0,ymax=17,fill="red",alpha=0.15) +
      geom_vline(aes(xintercept = 1),col="red",linetype="dashed") +
      xlab("Proportion") +
      theme(legend.position="none",
            panel.border = element_rect(colour = "black", fill=NA))
  }

  #______________
  
  if( grepl("indpt",simul_names[j] )==FALSE ){
    for(i in 1:n_rep){
  
      temp <-as.data.frame( name_simul[[i]]$BUGSoutput$sims.list$p_smolt_prod[,,4] )
      names(temp)<-seq(1,15,1)
      temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
      temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
      list_p_simul[[i]] <- temp
    }
  
    gg4 <- ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
      geom_density_ridges(scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
  
    for(i in 2:10){
       gg4 <- gg4 + geom_density_ridges(data=list_p_simul[[i]], aes(x=value,y=Year),scale=.95,alpha=0.2, rel_min_height = 0.001,col="#0072B250",fill="#0072B250")
    }
  
    p_rest <- as.data.frame(alpha_dir[,1:3]*y_wheels[,1:3])
    p_rest <- p_rest %>% mutate(rest=1-rowSums(.))
    p_rest $Year <- seq(1,15,1)
  
    gg4 <- gg4 + geom_point(data=p_rest,mapping=aes(x=rest, y= Year),col="red",size=1.5) +
      xlab("Proportion") +
      theme(legend.position="none",
            panel.border = element_rect(colour = "black", fill=NA))
  }

  ## Return proportions of smolt abundances
  ## in object list_pall_simul
  return_p_smolt(all_simul[[j]],j)
  
  sum_p_all <- list_pall_simul |> lmap(mean_p_all)
  #print(sum_p_all[[1]])
  #sum_p_all <- list_p3_simul |> lmap(mean_p_all)
  mean_p_greater_all_reps <- rbindlist(sum_p_all) %>% group_by(Year) %>% summarise(mean_p_greater1 = mean(mean_p_greater_1_pcent))
  
  ## Years with 3 upstream wheel
  years_n_rst <- rowSums(y_wheels[,1:3])
  col_rst <- c("#e66101","#fdb863","#5e3c99")
  
  ## Plotting the distribution of the proportion of all upstream RSTS
  gg_all <-  ggplot(list_pall_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+#,aes(height=1))+#..ndensity..))+#+
    geom_point(data=sum_p_all[[1]],mapping=aes(x=mean_p,y=Year),col=col_rst[years_n_rst])+
    xlim(0,3)+
    annotate("rect",xmin=1,xmax=3,ymin=0,ymax=17,fill="red",alpha=0.15) +
    geom_vline(aes(xintercept = 1),col="red",linetype="dashed") +
    #scale_x_continuous(breaks=c(0,0.25,0.50,0.75,1,1.25),labels=c(0,0.25,0.50,0.75,1,1.25))+
    labs(x="sum of upstream RSTs proportions")
  #labels = unit_format(unit = "M", scale = 1e-6))
 
  for(r in 2:10){
    gg_all <- gg_all + geom_density_ridges(data=list_pall_simul[[r]], aes(x=value,y=Year),rel_min_height = 0.001,col="#0072B250",fill="#0072B250",alpha=0.15,scale=.9)+
      geom_point(data=sum_p_all[[r]],mapping=aes(x=mean_p,y=Year),col=col_rst[years_n_rst])
  }
  
  gg_all <- gg_all +  annotate("text", x =2.5, y =(1:16)+0.2, label = c(round(mean_p_greater_all_reps$mean_p_greater1,1),"Avg. p > 1"))

 if( grepl("indpt",simul_names[j] )==FALSE ){
    ggarrange(gg1,gg2,gg3,gg4,
              labels = c("A", "B", "C","D"),
              ncol = 2, nrow = 2)
    }else{
      if(full_plot==T){
        ggarrange(gg1,gg2,gg3,gg_all,
                  labels = c("A", "B", "C","D"),
                  ncol = 2, nrow = 2)
      }else{
        ggarrange(gg1,gg2,gg3,
                  labels = c("A", "B", "C"),
                  ncol = 2, nrow = 2)
      }
    }
}


## Plots total abundance 
ggplot_totNm <- function(name_simul){
  
  list_p_simul <- list()

  for(i in 1:n_rep){
    if( grepl("indpt",deparse(substitute(jags_indpt_hier_no_split_other_wheels_simul)) ) ){
      temp <- as.data.frame(name_simul[[i]]$BUGSoutput$sims.list$Nm[,,1] / name_simul[[i]]$BUGSoutput$sims.list$Nm_tot)
    }else{
      temp <-as.data.frame( name_simul[[i]]$BUGSoutput$sims.list$p_smolt_prod[,,1] )
    }
    
    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    list_p_simul[[i]] <- temp
  }
  
  gg1 <- ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(scale=3,alpha=0.2, rel_min_height = 0.001)
  
  for(i in 2:10){
    gg1 <- gg1 + geom_density_ridges(data=list_p_simul[[i]], aes(x=value,y=Year),scale=3,alpha=0.2, rel_min_height = 0.001)
  }
  
  if( grepl("indpt",deparse(substitute(name_simul)) )==FALSE ){
    gg1 <- gg1 + geom_point(data=p_w,mapping=aes(x=p_1, y= Year),col="red",size=1.5) + 
      xlab("Proportion") +
      #xlim(0,2.5)+
      theme(legend.position="none")
  }else{
    
    #poly_temp <- data.frame(x=c(1,3,3,1),y=c(-0.5,-0.5,15.5,15.5))
    gg1 <- gg1 + geom_point(data=p_w,mapping=aes(x=p_1, y= Year),col="red",size=1.5) + 
      xlim(0,2.5)+
      annotate("rect",xmin=1,xmax=2.5,ymin=0,ymax=16,fill="red",alpha=0.15) +
      geom_vline(aes(xintercept = 1),col="red",linetype="dashed") +
      xlab("Proportion") +
      theme(legend.position="none")
  }

  #______________
  
  for(i in 1:n_rep){
    if( grepl("indpt",deparse(substitute(name_simul)) ) ){
      temp <- as.data.frame(name_simul[[i]]$BUGSoutput$sims.list$Nm[,,2] / name_simul[[i]]$BUGSoutput$sims.list$Nm_tot)
    }else{
      temp <-as.data.frame( name_simul[[i]]$BUGSoutput$sims.list$p_smolt_prod[,,2] )
    }
    
    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    
    years2keep <- y_wheels[,2]*1:15
    years2keep <- years2keep[years2keep!=0]
    
    temp <- temp %>% filter(Year %in% years2keep)
    
    yy <- 1:15
    yy <- yy[!yy %in% years2keep]
    
    for(k in 1:length(yy)){
      temp <-rbind(temp,c(yy[k],NA,NA))
    }
    list_p_simul[[i]] <- temp
  }

  gg2 <- ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(scale=2,alpha=0.2, rel_min_height = 0.001)
  
  for(i in 2:10){
    gg2 <- gg2 + geom_density_ridges(data=list_p_simul[[i]], aes(x=value,y=Year),scale=2,alpha=0.2, rel_min_height = 0.001)
  }

  if( grepl("indpt",deparse(substitute(name_simul)) )==FALSE ){
    gg2 <- gg2 + geom_point(data=p_w[years2keep,],mapping=aes(x=p_2, y= Year),col="red",size=1.5) +
      xlab("Proportion") +
      #xlim(0,5)+
      theme(legend.position="none")
  }else{
    gg2 <- gg2 + geom_point(data=p_w[years2keep,],mapping=aes(x=p_2, y= Year),col="red",size=1.5) +
      xlim(0,5)+
      annotate("rect",xmin=1,xmax=5,ymin=0,ymax=16,fill="red",alpha=0.15) +
      geom_vline(aes(xintercept = 1),col="red",linetype="dashed") +
      xlab("Proportion") +
      theme(legend.position="none")
  }

  #______________
 
  for(i in 1:n_rep){
    if( grepl("indpt",deparse(substitute(name_simul)) ) ){
      temp <- as.data.frame(name_simul[[i]]$BUGSoutput$sims.list$Nm[,,3] / name_simul[[i]]$BUGSoutput$sims.list$Nm_tot)
    }else{
      temp <-as.data.frame( name_simul[[i]]$BUGSoutput$sims.list$p_smolt_prod[,,3] )
    }
    
    names(temp)<-seq(1,15,1)
    temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
    temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
    
    years2keep <- y_wheels[,3]*1:15
    years2keep <- years2keep[years2keep!=0]
    
    temp <- temp %>% filter(Year %in% years2keep)
    
    yy <- 1:15
    yy <- yy[!yy %in% years2keep]
    
    for(k in 1:length(yy)){
      temp <-rbind(temp,c(yy[k],NA,NA))
    }
    list_p_simul[[i]] <- temp
  }

  gg3 <- ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
    geom_density_ridges(scale=2,alpha=0.2, rel_min_height = 0.001)
  
  for(i in 2:10){
    gg3 <- gg3 + geom_density_ridges(data=list_p_simul[[i]], aes(x=value,y=Year),scale=2,alpha=0.2, rel_min_height = 0.001)
  }
  
  gg3 <- gg3 + geom_point(data=p_w[years2keep,],mapping=aes(x=p_3, y= Year),col="red",size=1.5) + 
    xlab("Proportion") +
    theme(legend.position="none")
  
  #______________
  
  if( grepl("indpt",deparse(substitute(name_simul)) )==FALSE ){
    for(i in 1:n_rep){
      
      temp <-as.data.frame( name_simul[[i]]$BUGSoutput$sims.list$p_smolt_prod[,,4] )
      names(temp)<-seq(1,15,1)
      temp <- temp %>% pivot_longer(cols=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'), names_to="Year") %>% mutate(anchor=1)
      temp$Year<-factor(temp$Year,levels=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'))
      list_p_simul[[i]] <- temp
    }
    
    gg4 <- ggplot(list_p_simul[[1]], aes(x=value,y=Year)) +
      geom_density_ridges(scale=2,alpha=0.2, rel_min_height = 0.001)
    
    for(i in 2:10){
      
      gg4 <- gg4 + geom_density_ridges(data=list_p_simul[[i]], aes(x=value,y=Year),scale=2,alpha=0.2, rel_min_height = 0.001)
    }
    
    p_rest <- as.data.frame(alpha_dir[,1:3]*y_wheels[,1:3])
    p_rest <- p_rest %>% mutate(rest=1-rowSums(.))
    p_rest $Year <- seq(1,15,1)
    
    gg4 <- gg4 + geom_point(data=p_rest,mapping=aes(x=rest, y= Year),col="red",size=1.5) +
      xlab("Proportion") +
      theme(legend.position="none")
  }
 
  if( grepl("indpt",deparse(substitute(name_simul)) )==FALSE ){  
    ggarrange(gg1,gg2,gg3,gg4,
              labels = c("A", "B", "C","D"),
              ncol = 2, nrow = 2)
  }else{
    ggarrange(gg1,gg2,gg3,#gg4,
              labels = c("A", "B", "C"),
              ncol = 2, nrow = 2)
    
  }

}











