library(dplyr)
find_path<-system.file("app", "www", "mvidia/MVIDIA.py", package = "OpDEA")
root_folder<-gsub('MVIDIA.py', '', find_path)
source(paste0(root_folder,'utils.R'))

args <- commandArgs(trailingOnly = TRUE)
#Rscript /home/hui/PycharmProjects/DL_proteomics/benchmark_R_scripts/ensemble_set.R /data/res_DE_files//ensemle_res/CPTAC_DDA_maxquant_hurdle_.csv /data/res_DE_files/CPTAC_DDA/LTQW56/plgem/plgem_maxquant__MinDet_.csv /data/res_DE_files/CPTAC_DDA/LTQW56/DEqMS/DEqMS_maxquant__MLE_center.median.csv /data/res_DE_files/CPTAC_DDA/LTQW56/ProteoMM/ProteoMM_maxquant__bpca_vsn.csv
#Rscript /home/hui/PycharmProjects/DL_proteomics/benchmark_R_scripts/ensemble_set.R /data/res_DE_files//ensemle_res/CPTAC_DDA_LTQW56_maxquant_hurdle_.csv /data/res_DE_files/CPTAC_DDA/LTQW56/plgem/plgem_maxquant__MinDet_.csv /data/res_DE_files/CPTAC_DDA/LTQW56/DEqMS/DEqMS_maxquant__MLE_center.median.csv /data/res_DE_files/CPTAC_DDA/LTQW56/ProteoMM/ProteoMM_maxquant__bpca_vsn.csv
#Rscript /home/hui/PycharmProjects/DL_proteomics/benchmark_R_scripts/ensemble_set.R /data/res_DE_files//ensemle_res/CPTAC_DDA_LTQ86_fragpipe_hurdle_.csv /data/res_DE_files/CPTAC_DDA/LTQ86/plgem/plgem_fragpipe__QRILC_div.mean.csv /data/res_DE_files/CPTAC_DDA/LTQ86/DEP/DEP_fragpipe_MaxLFQ_MinProb_.csv /data/res_DE_files/CPTAC_DDA/LTQ86/ProteoMM/ProteoMM_fragpipe_MaxLFQ_MinDet_max.csv
#Rscript /data/proteomics/data/2022.12.8/Hurdle.R /data/proteomics/data/2022.12.8/temp7/CPTAC_DDA/LTQ86/edgeR/  /data/proteomics/data/2022.12.8/temp7/CPTAC_DDA/LTQ86/edgeR/edgeR_maxquant___.csv /data/proteomics/data/2022.12.8/temp7/CPTAC_DDA/LTQ86/limma/limma_maxquant_LFQ__.csv

# args<-c('/data/res_DE_files//ensemle_res/CPTAC_DDA_LTQ86_fragpipe_hurdle_.csv',
#         '/data/res_DE_files/CPTAC_DDA/LTQ86/plgem/plgem_fragpipe__QRILC_div.mean.csv',
#         '/data/res_DE_files/CPTAC_DDA/LTQ86/DEP/DEP_fragpipe_MaxLFQ_MinProb_.csv',
#         '/data/res_DE_files/CPTAC_DDA/LTQ86/ProteoMM/ProteoMM_fragpipe_MaxLFQ_MinDet_max.csv'
#  )

# args<-c('T', 'HUMAN;YEAST;ECOLI', 'YEAST;ECOLI',
# 'D:/data/benchmark/data/dataset_info/PXD028735_true_fc.txt',
# 'T', 'E:/multi-view_AE_DEA/res_test/temp0/hurdle_real.csv',
# 'E:/multi-view_AE_DEA/res_test/temp21_limma_min_FragPipe/HYE5600735_LFQ/HYE5600735_LFQ_limma_FragPipe_top0__.csv',
# 'E:/multi-view_AE_DEA/res_test/temp21_limma_min_FragPipe/HYE5600735_LFQ/HYE5600735_LFQ_limma_FragPipe_maxlfq__.csv'
#  )




N<-length(args)

save_fold<-args[6]
print_lab<-args[1]
true_organism<-args[2]
DE_organism<-args[3]
true_lgfc<-args[4]
logT<-args[5]
r1<-args[7]


dnum<-function(X){
  n<-nchar(X)
  #if(grepl('_HUMAN', X) | grepl('_YEAST', X) | grepl('_ECOLI', X) | grepl('_UPS', X)){
  if(grepl('_', X)){
    sps<-strsplit(X, '_', fixed = TRUE)
    sps[[1]][length(sps[[1]])] = gsub("\\d","",sps[[1]][length(sps[[1]])])
    outX = sps[[1]][1]
    for (i in 2:length(sps[[1]])) {
      outX = paste0(outX, '_', sps[[1]][i])
    }
  }else{
    outX = X
  }
  #outX = X
  #outX = gsub("\\d","",X)
  # if((substr(X, n, n)!="T") & (substr(X, n, n)!='N') & (substr(X, n, n)!='S') & (substr(X, n, n)!='I')){
  #   outX<-substr(X, 1, n-1)
  # }
  return(outX)
}

res1<-read.table(r1, sep = ',', header = TRUE)#, quote = "")
res1$protein<-apply(as.array(res1$protein), 1, function(x) dnum(x))
base.vars <- c("protein", "contrast")
# Join everything together
b1<-data.frame(cbind(res1$protein, res1$contrast))
colnames(b1)<-c('protein', 'contrast')
res.Base<-b1

ress<-list()

for (i in 8:N) {
  res<-read.table(args[i], sep = ',', header = TRUE)#, quote = "")
  res$protein<-apply(as.array(res$protein), 1, function(x) dnum(x))
  ress[[i]]<-res
  b<-data.frame(cbind(res$protein, res$contrast))
  colnames(b)<-c('protein', 'contrast')
  #base.vars <- c("protein", "UPS", "contrast")
  res.Base<-full_join(res.Base,b,base.vars)
}

###
# plgem pvalue logFC
# ANOVA pvalue logFC
# ProteoMM pvalue logFC
# beta_binomial pvalue logFC
# DEP DEqMS edgeR limma msqrob2 MSStats proDA ROTS SAM siggenes ttest
res11 <- full_join(res.Base, res1)
z11<- -qnorm(as.numeric(res11$pvalue)/2)*sign(as.numeric(res11$logFC))

res.Hurdle <- res.Base
#res.Hurdle$logFC_r1<-res11$logFC
logFCs<-as.numeric(res11$logFC)
dfs<-(!is.na(z11))

all_z<-z11^2
for (i in 8:N) {
  resi <- full_join(res.Base, ress[[i]],base.vars)
  zi<- -qnorm(as.numeric(resi$pvalue)/2)*sign(as.numeric(resi$logFC))
  all_z<-cbind(all_z,zi^2)
  logFCs<-cbind(logFCs, as.numeric(resi$logFC))
  dfs = dfs + (!is.na(zi))
}
logFCs[is.na(logFCs)]<-0
res.Hurdle$chisq <- all_z %>% rowSums(na.rm = TRUE)
res.Hurdle$pensemble <- 1 - pchisq(res.Hurdle$chisq, df = dfs)
#res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) mean(x[which(!is.na(x) & !is.infinite(x))]))
# biggest |logFC|
#res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) sign(x[which(abs(x[which(!is.na(x) & !is.infinite(x))])==max(abs(x[which(!is.na(x) & !is.infinite(x))])))]) * max(abs(x[which(!is.na(x) & !is.infinite(x))])))
res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) sign(x[which(abs(x[which(!is.na(x) & !is.infinite(x))])==max(abs(x[which(!is.na(x) & !is.infinite(x))])))[1]]) * max(abs(x[which(!is.na(x) & !is.infinite(x))]))[1])

res.Hurdle <- res.Hurdle %>% group_by(contrast) %>% mutate(qensemble = pensemble %>%  p.adjust(method = "BH")) %>% ungroup()
res.Hurdle <- res.Hurdle %>% arrange(pensemble)
#res.Hurdle$UPS<-grepl('_HUMAN',res.Hurdle$protein)
colnames(res.Hurdle)<-c('protein', 'contrast', 'chisq', 'pvalue', 'logFC', 'adj.pvalue')#, 'UPS')
res_all<-res.Hurdle
if(print_lab=='T'){
  true_organisms=strsplit(true_organism,';',fixed = T)[[1]]
  DE_organisms=strsplit(DE_organism,';',fixed = T)[[1]]
  labs<-get_labels(res_all, true_organisms, DE_organisms, true_lgfc)
  res_all$organism<-labs$Organism
  res_all$label<-labs$DEP
  res_all$TlogFC<-labs$TlogFC
}

write.table(res_all, save_fold, sep = ',', row.names = TRUE, col.names = TRUE)
