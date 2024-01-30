# coloc 共定位分析

## 加载需要的包
```
library(data.table)
library(dplyr)
library(coloc)
library(MungeSumstats)
```

## 处理eqtl数据，以肺eqtl为例
``` R
exe <- "S:/胸外/twaslc/twaslca/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe"
beqtl <-  "S:/胸外/twaslc/Lung_eqtl/Lung"
out_path <- "S:/胸外/twaslc/Lung_eqtl2"
query=1
code<-paste(exe,"--beqtl-summary",beqtl,"--query",query,"--out",out_path)
system(code)
```

> [!IMPORTANT]
> 注意Windows和Linux的反斜杠处理
### 只有Windows才需要处理
``` R
exe<-gsub("/","\\\\",exe)
beqtl<-gsub("/","\\\\",beqtl)
out_path<-gsub("/","\\\\",out_path)
```



## 处理GWAS数据，以肺癌为例，
> [!IMPORTANT]
>注意版本与eqtl一致，否则需要转换版本
```R
lcgwas<-fread("S:/胸外/twaslc/twaslca/hg37_lcgwas_GCST90043864.tsv")
id.e <- "disease_lc"
lcgwas$FRQ<-ifelse(lcgwas$effect_allele_frequency>=0.5,1-lcgwas$effect_allele_frequency,lcgwas$effect_allele_frequency)

expo<-dplyr::select(lcgwas,SNP=variant_id,CHR=chromosome,BP=base_pair_location,
                    A1=effect_allele,A2=other_allele,BETA=beta,
                    SE=standard_error,P=p_value,EAF=effect_allele_frequency,FRQ,N)
expo <- expo%>%dplyr::mutate(Z=BETA/SE)
expo <- na.omit(expo, cols = "SNP")
expo <- as_tibble(expo)
```

## 继续整理eqtl数据
``` R
Lung_eqtl <- fread("S:/胸外/twaslc/Lung_eqtl/Lung_eqtl.txt")
id.o <- "gene"
outcome<-dplyr::select(Lung_eqtl,Gene=Probe,SNP,CHR=Chr,BP,A1,A2,BETA=b,SE=SE,P=p)
outcome$N <- 450000
outcome <- outcome%>%dplyr::mutate(Z=BETA/SE)
outcome <- na.omit(outcome, cols = "SNP")
genelist<-unique(outcome$Gene)
```

## coloc分析，当p值为0时使用最小浮点数
```R
res <- list()
for (i in  genelist[659:665]) {
  
  
  outcomegene <- outcome[outcome$Gene == i,]
  outcomegene <- as_tibble(outcomegene)
  
  expogene <- expo[expo$SNP %in% outcomegene$SNP,]
  outcomegene <- outcomegene[outcomegene$SNP %in% expogene$SNP ,]
  outcomegene <- outcomegene %>%
    left_join(expogene %>% select(SNP, FRQ), by = "SNP")
  
 
  if (nrow(outcomegene) == 0) {
    next
  }
  
  expogene$P <- ifelse(expogene$P <= 0, .Machine$double.eps, expogene$P)
  outcomegene$P <- ifelse(outcomegene$P <= 0, .Machine$double.eps, outcomegene$P)
  
  
  datalist_names <- c("snp","position","beta","type","MAF","N","s","pvalues")
  eqtl_data_list <- vector("list", length(datalist_names))
  names(eqtl_data_list) <- datalist_names
  eqtl_data_list$snp <- expogene$SNP
  eqtl_data_list$position <- expogene$BP
  eqtl_data_list$beta <- expogene$BETA
  eqtl_data_list$pvalues <- expogene$P
  eqtl_data_list$MAF <- expogene$FRQ
  eqtl_data_list$s <- 0.339
  eqtl_data_list$N <- expogene$N
  eqtl_data_list$type <- "cc"
  # check_dataset(eqtl_data_list)
  
  
  cc_list_names <- c("snp","position","beta","type","sdY","N", "MAF","pvalues")
  gwas_data_list <- vector("list", length(cc_list_names))
  names(gwas_data_list) <- cc_list_names
  gwas_data_list$snp <- outcomegene$SNP
  gwas_data_list$position <- outcomegene$BP
  gwas_data_list$beta <- outcomegene$BETA
  gwas_data_list$pvalues <- outcomegene$P
  gwas_data_list$MAF <- outcomegene$FRQ
  gwas_data_list$sdY <- 1
  #cc 就是 case control 的意思
  #s 是病例占总数的百分比
  gwas_data_list$N <- outcomegene$N
  gwas_data_list$type <- "quant"
  # check_dataset(gwas_data_list)
  # Basic coloc
  
  
  abf_res <- coloc.abf(dataset1=eqtl_data_list,
                       dataset2=gwas_data_list)
  
  res[[i]] <- abf_res
}

library(clusterProfiler)
result <- data.frame()
for (i in names(res)) {
  mid<-t(data.frame(res[[i]]$summary))
  result <- rbind(result,mid)
}
```
## 后续处理,匹配基因及筛选PP.H3+PP.H4>0.8
``` R
rownames(result)<-names(res)
result$ENSEMBL<-rownames(result)
saveRDS(result,"result_coloc_lc.rds")

result_coloc_lc<-readRDS("result_coloc_lc.rds")
result_coloc_lc<-result%>%left_join(bitr(result$ENSEMBL,"ENSEMBL","SYMBOL","org.Hs.eg.db"),by="ENSEMBL")

result_coloc_lc$'PP.H3+pp.H4'<-result_coloc_lc$PP.H3.abf+result2$PP.H4.abf
result.sig<-subset(result_coloc_lc,result_coloc_lc$`PP.H3+pp.H4`>0.8)
```

## 自己保存结果
