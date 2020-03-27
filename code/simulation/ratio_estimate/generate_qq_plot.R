theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_line(colour="#f0f0f0"),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold.italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#FFD042","#7fc97f","#004953","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#EF7E3D","chartreuse4","#7fc97f","#004953","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


#plot the ratio estimate distribution
n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.00,0.01,0.03,0.05)
beta_vec <- c(0,0.3,0.5,1)
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
i4 = as.numeric(args[[3]])
setwd("/data/zhangh24/MR_MA/")
load(paste0("./result/simulation/ratio_estimate/ratio_estimate_merged.Rdata"))
#i1 correponding to n
#i2 corresponding to alpha
library(ggplot2)
p <- list()
  alpha_U = 0.1
  beta_U <- 0.1
  sigma_y = 1
  sigma_m = 1
 
       
      temp = 12*(i4-1)+4*(i1-1)+i2
      n <- n_vec[i1]
      alpha_G = alpha_vec[i2]
      beta_M = beta_vec[i4]
      result <- result_final[[temp]]
      Gamma = result[[1]]
      var_Gamma = result[[2]]
      gamma = result[[3]]
      var_gamma = result[[4]]
      ratio_est = result[[5]]
      var_ratio <- result[[6]]
      
      true_distribution = ratio_est
      dot <- seq(0.01,0.99,by=0.002)
      q_true <- quantile(true_distribution,dot)
      #pdot <-ifelse(dot<=0.5,2*dot,2*(1-dot))
      
      
      n.simu = 1e6
      z_Gamma <- rnorm(n.simu,mean =alpha_G*beta_M,sd =sqrt((sigma_y+beta_M^2*sigma_m+(beta_M*alpha_U)^2)/n))
      z_gamma <- rnorm(n.simu,mean = alpha_G,sd = sqrt((alpha_U^2+sigma_m)/n))
      true_distribution2 <- z_Gamma/z_gamma
      q_true2 = quantile(true_distribution2,dot)
      # dot2 = dot
      # for(i in 1:length(dot)){
      #   dot2[i] <- sum(true_distribution2<=q_true[i])/length(true_distribution2)
      #   
      # }
      #pdot2 <- ifelse(dot2<=0.5,2*dot2,2*(1-dot2))
      data <- data.frame(q_true,
                         q_true2)
      colnames(data) <- c("q_true",
                          "q_true2")
     # p[[1]] <- ggplot(data,aes(q_true,q_true2))+
     #    geom_point()+
     #    geom_abline(slope=1)+
     #    xlab(paste0("Empirical distribution quantile"))+
     #    ylab(paste0("Derived distribution quantile"))+
     #    ggtitle(paste0("n=",n,", beta = ",beta_M,", alpha =",alpha_G))+
     #    theme_Publication()
      p[[1]] <- ggplot(data,aes(pdot,pdot2))+
        geom_point()+
        geom_abline(slope=1)+
        xlab(paste0("Empirical distribution"))+
        ylab(paste0("Derived distribution"))+
        ggtitle(paste0("n=",n,", beta = ",beta_M,", alpha =",alpha_G))+
        theme_Publication()+
        scale_y_continuous(limits=c(0,1))
      
      
     n.simu = 1e6
     true_distribution3 = rep(0,n.simu)
     total <- 0
     for(k in 1:1){
       z_Gamma <- rnorm(n.simu,mean =Gamma[k],sd =sqrt(var_Gamma[k]))
       z_gamma <- rnorm(n.simu,mean = gamma[k],sd = sqrt(var_gamma[k]))
       temp_simulation <- z_Gamma/z_gamma  
       true_distribution3[total+(1:n.simu)] <- temp_simulation
       total = total+n.simu
     }
     q_true3 = quantile(true_distribution3,dot)
     data <- data.frame(q_true,
                        q_true3)
     colnames(data) <- c("q_true",
                         "q_true2")
     p[[2]] <- 
      # png(paste0("./result/simulation/ratio_estimate/qqplot/test.png"),width = 16,height = 8,
         #  unit = "in",res = 300)
     ggplot(data,aes(q_true,q_true2))+
       geom_point()+
       geom_abline(slope=1)+
       xlab(paste0("Empirical distribution"))+
       ylab(paste0("Plug in distribution"))+
       ggtitle(paste0("n=",n,", beta = ",beta_M,", alpha =",alpha_G))+
       theme_Publication()
      # dev.off()
     
     
     
     
     
     
     
     
     
     
     
     
     
     
      n.simu = 1e3
      true_distribution3 = rep(0,n.simu*length(Gamma))
      total <- 0
      for(k in 1:length(Gamma)){
        z_Gamma <- rnorm(n.simu,mean =Gamma[k],sd =sqrt(var_Gamma[k]))
        z_gamma <- rnorm(n.simu,mean = gamma[k],sd = sqrt(var_gamma[k]))
        temp_simulation <- z_Gamma/z_gamma  
        true_distribution3[total+(1:n.simu)] <- temp_simulation
        total = total+n.simu
      }
      
      dot3 = dot
      q_true3 = quantile(true_distribution3,dot)
      data <- data.frame(q_true,
                         q_true3)
      colnames(data) <- c("q_true",
                          "q_true2")
      
      p[[2]] <- ggplot(data,aes(pdot,pdot2))+
        geom_point()+
        geom_abline(slope=1)+
        xlab(paste0("Empirical distribution"))+
        ylab(paste0("Plug in distribution"))+
        ggtitle(paste0("n=",n,", beta = ",beta_M,", alpha =",alpha_G))+
        theme_Publication()+
        scale_y_continuous(limits=c(0,1))
      
   
library(gridExtra)
png(paste0("./result/simulation/ratio_estimate/qqplot/qq_",i1,"_",i2,"_",i4,".png"),width = 16,height = 8,
    unit = "in",res = 300)
grid.arrange(p[[1]],p[[2]],nrow=1)
dev.off()

