n.test = 100000
prs_m_mat = matrix(0,n.test,n.rep)
prs_y_mat = matrix(0,n.test,n.rep)

for(i1 in 1:30){
  print(i1)
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/prs_result_twostage_",i1,".rdata"))  
  prs_m_mat = prs_m_mat+prs_result[[1]]
  prs_y_mat = prs_y_mat+prs_result[[2]]
}
prs_result = list(prs_m_mat,prs_y_mat)
save(prs_result,file =paste0("/data/zhangh24/MR_MA/result/simulation/prs/prs_result_combined_twostage.rdata"))

system(paste0("rm /data/zhangh24/MR_MA/result/simulation/prs/prs_result_twostage_*"))
# prs_result <- rep(0,n.rep)
# prs_low_result <- rep(0,n.rep)
# prs_high_result <- rep(0,n.rep)
# for(l in 1:n.rep){
#   model <- lm(prs_y_mat[,l]~prs_m_mat[,l])
#   prs_result[l] <- coefficients(summary(model))[2,1]
#   temp_ci <- confint(model)
#   prs_low_result[l] <- temp_ci[2,1]
#   prs_high_result[l] <- temp_ci[2,2]
# }

