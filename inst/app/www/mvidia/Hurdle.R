library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

N<-length(args)

save_fold<-args[1]
r1<-args[2]


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

for (i in 3:N) {
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
z11<- -qnorm(res11$pvalue/2)*sign(res11$logFC)

res.Hurdle <- res.Base
#res.Hurdle$logFC_r1<-res11$logFC
logFCs<-res11$logFC
dfs<-(!is.na(z11))

all_z<-z11^2
for (i in 3:N) {
  resi <- full_join(res.Base, ress[[i]])
  zi<- -qnorm(resi$pvalue/2)*sign(resi$logFC)
  all_z<-cbind(all_z,zi^2)
  logFCs<-cbind(logFCs, resi$logFC)
  dfs = dfs + (!is.na(zi))
}

res.Hurdle$chisq <- all_z %>% rowSums(na.rm = TRUE)
res.Hurdle$pensemble <- 1 - pchisq(res.Hurdle$chisq, df = dfs)
#res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) mean(x[which(!is.na(x) & !is.infinite(x))]))
# biggest |logFC|
#res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) sign(x[which(abs(x[which(!is.na(x) & !is.infinite(x))])==max(abs(x[which(!is.na(x) & !is.infinite(x))])))]) * max(abs(x[which(!is.na(x) & !is.infinite(x))])))
res.Hurdle$logFC_avg<-apply(logFCs,1,function(x) sign(x[which(abs(x[which(!is.na(x) & !is.infinite(x))])==max(abs(x[which(!is.na(x) & !is.infinite(x))])))[1]]) * max(abs(x[which(!is.na(x) & !is.infinite(x))]))[1])

res.Hurdle <- res.Hurdle %>% group_by(contrast) %>% mutate(qensemble = pensemble %>%  p.adjust(method = "BH")) %>% ungroup()
res.Hurdle <- res.Hurdle %>% arrange(pensemble)
#res.Hurdle$UPS<-grepl('_HUMAN',res.Hurdle$protein)  
colnames(res.Hurdle)<-c('protein', 'contrast', 'chisq', 'pvalue', 'logFC', 'adj.pvalue')
res_all<-res.Hurdle

write.table(res_all, save_fold, sep = ',', row.names = TRUE, col.names = TRUE)