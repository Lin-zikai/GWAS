# SMR笔记
## smr参考来源：https://yanglab.westlake.edu.cn/software/smr/#Overview

在执行SMR分析时需要以下三个文件包括：

### GWAS-summary 
需要gwas.ma格式文件，一般为GRCH38，如下所示
> [!CAUTION]
> 这个的数据要注意使用的基因组版本
```
SNP    A1  A2  freq    b   se  p   n
rs1001    A   G   0.8493  0.0024  0.0055  0.6653  129850
rs1002    C   G   0.03606 0.0034  0.0115  0.7659  129799
rs1003    A   C   0.5128  0.045   0.038   0.2319  129830
......
```

清洗数据的r脚本如下
``` R
library(data.table)
library(dplyr)
library(MungeSumstats)
colon <- fread("lcgwas_GCST90043864.tsv.gz") #原始文件
head(colon)
data<-dplyr::select(colon,CHR=chromosome,BP=base_pair_location,A1=effect_allele,A2=other_allele,
                    FRQ=effect_allele_frequency,BETA=beta,SE=standard_error,P=p_value,SNP=variant_id,n=N)
dat<-liftover(data,convert_ref_genome = "GRCh38",ref_genome = "GRCh37")  #修改基因组版本

dat <- dat[dat$IMPUTATION_gen_build == "TRUE",]
dat<-dat%>%dplyr::select(SNP,A1,A2,freq=FRQ,
                           b=BETA,se=SE,p=P,n)    #清洗为需要的格式
fwrite(dat,"pap/lungcancer.ma",sep = " ")    #保存数据
```

### 参考人群文件 bed bim 和 fam
千人计划参考基因组
可以从这个下载：http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
> [!CAUTION]
> 这个是使用GRCH37构建的  经过尝试好像用在38没什么影响，可能是由于rsid的原因 如果有大佬懂还请指点一下

这里有GRCH38的：https://www.internationalgenome.org/data-portal/data-collection/30x-grch38

### BESD 格式的 eQTL文件

具体如何进行格式转换可以查看：https://yanglab.westlake.edu.cn/software/smr/#BESDformat



