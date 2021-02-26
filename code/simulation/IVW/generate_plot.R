getwd()
data.IVW = read.csv("./result/simulation/IVW/cover_IVW_table.csv",header=T)
total = nrow(data.IVW)*3*2
beta_vec = rep(0,total)
alpha_vec = rep(0,total)
sample_size_vec = rep(0,total)
method_vec = rep("c",total)
IV_vec = rep(0,total)
coverage_vec = rep(0,total)
temp = 1
sample_sze = c(15000,75000,150000)
for(i in 1:18){
  for(j in 1:3){
    beta_vec[temp] = data.IVW[i,1]
    alpha_vec[temp] = data.IVW[i,2]
    sample_size_vec[temp] = sample_sze[j]
    coverage_vec[temp] = data.IVW[i,j+2]
    method_vec[temp] = "IVW"
    IV_vec[temp] = data.IVW[i,6]
    temp = temp+1
  }
  
}
  
data.IVW = read.csv("./result/simulation/IVW/cover_MRweight_table.csv",header=T)
for(i in 1:18){
  for(j in 1:3){
    beta_vec[temp] = data.IVW[i,1]
    alpha_vec[temp] = data.IVW[i,2]
    sample_size_vec[temp] = sample_sze[j]
    coverage_vec[temp] = data.IVW[i,j+2]
    method_vec[temp] = "MR-Weight"
    IV_vec[temp] = data.IVW[i,6]
    temp = temp+1
  }
  
}



plot.data = data.frame(beta_vec,alpha_vec,sample_size_vec,coverage_vec,method_vec,IV_vec)
for(IV in c(1,5)){
  #generage data for IVW only
  plot.data.sub = plot.data %>% 
    filter(method_vec =="IVW"&
             IV_vec ==IV) %>% 
    mutate(sample_size = factor(sample_size_vec,levels = c("15000","75000","150000")),
           beta = paste0("beta = ",beta_vec),
           alpha = paste0("alpha = ",alpha_vec),
           Method = method_vec)
  p = ggplot(data = plot.data.sub, aes(x=sample_size,y = coverage_vec))+
    geom_bar(aes(fill=Method),
             stat = "identity",
             position="dodge")+
    facet_grid(vars(beta),vars(alpha))+
    theme_Publication()+
    xlab("Sample size")+
    ylab("Coverage")+
    coord_cartesian(ylim=c(0.90,1))+
    scale_fill_Publication()+
    geom_abline(slope = 0,intercept = 0.95,linetype = "dashed",color = "black")+
    theme(legend.position = "none")
  
  png(paste0("./result/simulation/IVW/coverage_IVW_IV_",IV,".png"),height = 6, width = 8,units = "in",res = 300)
  print(p)
  dev.off()
  
  plot.data.sub = plot.data %>% 
    filter(IV_vec ==IV) %>% 
    mutate(sample_size = factor(sample_size_vec,levels = c("15000","75000","150000")),
           beta = paste0("beta = ",beta_vec),
           alpha = paste0("alpha = ",alpha_vec),
           Method = method_vec)
  p = ggplot(data = plot.data.sub, aes(x=sample_size,y = coverage_vec))+
    geom_bar(aes(fill=Method),
             stat = "identity",
             position="dodge")+
    facet_grid(vars(beta),vars(alpha))+
    theme_Publication()+
    xlab("Sample size")+
    ylab("Coverage")+
    coord_cartesian(ylim=c(0.90,1))+
    scale_fill_Publication()+
    geom_abline(slope = 0,intercept = 0.95,linetype = "dashed",color = "black")
  
  png(paste0("./result/simulation/IVW/coverage_all_IV_",IV,".png"),height = 6, width = 8,units = "in",res = 300)
  print(p)
  dev.off()
  
  
}

    
