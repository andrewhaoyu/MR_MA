pthres = c(2.895529e-06,1E-05,5E-05,1E-04,5E-04,1E-03)
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
for(i1 in 1:length(pthres)){
  load(paste0(cur.dir,"WMR_result_chr22_i1",i1,".rdata"))
}