# This one is mainly about how to do phewas

> [!IMPORTANT]
> It is recommended to run on a server

## 数据来源
我个人比较喜欢的phewas数据库是https://pheweb.org/UKB-TOPMed/ ，我将使用这个进行举例

> [!NOTE]
CKB也有一个https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10435379/，不过我还没拿到权限，老是忘记去问

首先先下载数据，有个大几百份数据的，具体的数据可以看 https://github.com/Lin-zikai/GWAS/blob/main/data/phenotypes%20.tsv 或者自己在网站中下载也行

#现在正式开始进行phewas了
##首先你需要准备一份MR分析的暴露文件

格式如下，我将其命名为expdat.csv,最后一列是蛋白的名字，可以放多个蛋白
``` 
SNP	chr.exposure	pos.exposure	effect_allele.exposure	other_allele.exposure	eaf.exposure	beta.exposure	se.exposure	pval.exposure	exposure	mr_keep.exposure	pval_origin.exposure	id.exposure
rs148574467	2	238511920	A	G	NA	-0.357779	0.10503	0.000658183	exposure	TRUE	reported	HES6
rs6749854	2	239113126	A	G	NA	-0.126647	0.0356435	0.000380635	exposure	TRUE	reported	HES6
rs62194936	2	239174693	C	T	NA	-0.138153	0.0380481	0.000282311	exposure	TRUE	reported	HES6
rs3791454	2	240049647	C	T	NA	0.104537	0.0270646	0.000112236	exposure	TRUE	reported	HES6
rs7672991	4	41910652	T	C	NA	-0.283449	0.0800264	0.000397187	exposure	TRUE	reported	LIMCH1
rs8866	17	65373979	G	C	NA	-0.120051	0.0223952	8.30E-08	exposure	TRUE	reported	PITPNC1
```
不知道怎么制作这份expdat.csv文件话可以使用下面这个代码


现在你应该是已经获得了你需要分析的蛋白名称，你只需要把名字放到genelist中

> [!NOTE]
看不懂下面这个代码的话，可以看https://github.com/Lin-zikai/GWAS/blob/main/local_clump.md

``` R
library(data.table)
library(parallel)
library(TwoSampleMR)
library(ieugwasr)

genelist <- c("HES6", "MIR4435-2HG","LIMCH1","POU2AF2","SMAD9","PITPNC1","GDF15",
              "TGFBR2", "PRKAA1", "COX5B", "COX14", "ATG10", "GRN", "SORT1", "GUCY2C", "GUCA2A") #这里哦

disease <- list.files("./", pattern = "\\.gz$", full.names = TRUE)

ld_clump_local <- function (dat, clump_kb, clump_r2, clump_p, bfile, plink_bin) {
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", 
                  "sh")
  fn <- tempfile()
  write.table(data.frame(SNP = dat[["rsid"]], P = dat[["pval"]]), 
              file = fn, row.names = F, col.names = T, quote = F)
  fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --clump ", shQuote(fn, 
                                                                    type = shell), " --clump-p1 ", clump_p, " --clump-r2 ", 
                 clump_r2, " --clump-kb ", clump_kb, " --out ", shQuote(fn, 
                                                                        type = shell))
  system_output <- try(system(paste(fun2, "2>&1"), intern = TRUE))
  # Check for the specific warning message
  warning_message <- "Warning: No significant --clump results.  Skipping."
  if (warning_message %in% system_output) {
    message("Encountered the specific warning message: ", warning_message)
    return(data.frame()) # return an empty data frame
  } 
  res <- read.table(paste(fn, ".clumped", sep = ""), header = T)
  unlink(paste(fn, "*", sep = ""))
  y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
  if (nrow(y) > 0) {
    message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), 
            " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
}

# gene <- fread("/GPUFS/gyfyy_jxhe_1/User/zyt/coloc/Colon_Transverse/Colon_Transverse/colon.txt")

gene <- fread("/GPUFS/gyfyy_jxhe_1/User/lzk/gtes/Whole_Blood.txt") #这个是你的eqtl文件

exp <- data.table()
for (i in  genelist){
        
genename <- i
gene1 <- gene[gene$Gene == genename,]
gene1 <- gene1[gene1$p<1E-03]

if (nrow(gene1) == 0 ) {
  print(paste0(i," has no snp"))
  next
} 


exp_dat <- format_data(
  dat=gene1,
  type = 'exposure',
  snp_col = "SNP",
  beta_col = "b",
  se_col = "SE",
  eaf_col = "Freq",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  # samplesize_col = "all_meta_sample_N",
  pval_col = "p",
  chr_col = "Chr",
  pos_col = "BP"
)
ld <- ld_clump_local(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure),
                     clump_r2=0.001,clump_kb=10000,clump_p=1, bfile = "/GPUFS/gyfyy_jxhe_1/User/zyt/coloc/1000genomes/EUR", plink_bin = "/GPUFS/gyfyy_jxhe_1/User/zyt/coloc/plink_linux_x86_64_20231018/plink")
exp_dat <- exp_dat[exp_dat$SNP %in% ld$rsid, ]
if (nrow(exp_dat) == 0 ) {
      next
} 

exp_dat$id.exposure <- i
if (nrow(exp_dat) != 0 ) {
      print(paste0(i,"'s snp is ",exp_dat$SNP))
} else {
      print(i,"has no snp")
}
exp <- rbind(exp,exp_dat)
}


write.csv(exp,"expdat.csv")

```













# 开摆过几天写，过年啦！！！
