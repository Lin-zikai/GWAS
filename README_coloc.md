#coloc 共定位分析

##加载需要的包
`library(data.table)
library(dplyr)
library(coloc)
library(MungeSumstats)`

##处理eqtl数据，以肺eqtl为例
`exe <- "S:/胸外/twaslc/twaslca/smr-1.3.1-win-x86_64/smr-1.3.1-win.exe"
beqtl <-  "S:/胸外/twaslc/Lung_eqtl/Lung"
out_path <- "S:/胸外/twaslc/Lung_eqtl2"
query=1
code<-paste(exe,"--beqtl-summary",beqtl,"--query",query,"--out",out_path)
system(code)`

###注意Windows和Linux的反斜杠处理，只有Windows才需要处理
`exe<-gsub("/","\\\\",exe)
beqtl<-gsub("/","\\\\",beqtl)
out_path<-gsub("/","\\\\",out_path)`






##处理GWAS数据，以肺癌为例，注意版本与eqtl一致，否则需要转换版本

`lcgwas<-fread("S:/胸外/twaslc/twaslca/hg37_lcgwas_GCST90043864.tsv")`
`id.e <- "disease_lc"`
`lcgwas$FRQ<-ifelse(lcgwas$effect_allele_frequency>=0.5,1-lcgwas$effect_allele_frequency,lcgwas$effect_allele_frequency)`

`expo<-dplyr::select(lcgwas,SNP=variant_id,CHR=chromosome,BP=base_pair_location,
                    A1=effect_allele,A2=other_allele,BETA=beta,
                    SE=standard_error,P=p_value,EAF=effect_allele_frequency,FRQ,N)`
`expo <- expo%>%dplyr::mutate(Z=BETA/SE)`
