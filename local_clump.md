# 记录如何本地clump
## 使用的函数可以看这个https://mrcieu.github.io/ieugwasr/reference/ld_clump_local.html

## 加载需要的包
``` R
library(TwoSampleMR)
library(data.table)
library(ieugwasr)
library(plinkbinr)  #这个包是帮你装个plink的，自己装plink也行，就是要注意运行时是否有权限
```

## 读取数据
```
gene_blood_exp_snp<-readRDS("gene_blood_exp_snp.rds") #eQTL数据
snp_epilepsy_out<-readRDS("snp_epilepsy_out.rds")   #你的gwas数据
exp_dat1<-gene_blood_exp_snp
genelist <- unique(exp_dat1$gene)
```

```
head(gene_blood_exp_snp)
 gene hgnc_symbol        SNP chr.exposure pos.exposure effect_allele.exposure other_allele.exposure eaf.exposure z.exposure
1 ENSG00000000938         FGR rs34806307            1     27958339                      T                     C   0.05812083    11.1724
2 ENSG00000000938         FGR  rs4908343            1     27931698                      G                     A   0.18086010     8.3407
3 ENSG00000000938         FGR rs12726763            1     27958245                      G                     A   0.94275336     7.9185
4 ENSG00000000938         FGR rs12749647            1     27895517                      G                     A   0.11428375     6.9687
5 ENSG00000000938         FGR rs72886321            1     27909007                      A                     C   0.01914212     6.5440
6 ENSG00000000938         FGR rs59900503            1     27907201                      C                     T   0.01919834     6.5196
  beta.exposure se.exposure pval.exposure id.exposure exposure         F
1    0.18878560  0.01689750    5.5656e-29       blood exposure 124.82252
2    0.08574590  0.01028042    7.3915e-17       blood exposure  69.56728
3    0.13488910  0.01703468    2.4059e-15       blood exposure  62.70264
4    0.08669946  0.01244127    3.2011e-12       blood exposure  48.56278
5    0.18905513  0.02888984    5.9975e-11       blood exposure  42.82394
6    0.18808058  0.02884848    7.0495e-11       blood exposure  42.50518

head(snp_epilepsy_out)
chr.outcome pos.outcome        SNP effect_allele.outcome other_allele.outcome eaf.outcome z.outcome beta.outcome  se.outcome outcome
1           6   130840091  rs2326918                     A                    G      0.8470    -0.279 -0.002586613 0.009271014 outcome
2           7   145771806  rs6977693                     T                    C      0.8587    -0.874 -0.008373962 0.009581192 outcome
3          11   100009976 rs12364336                     A                    G      0.8772    -1.506 -0.015313714 0.010168469 outcome
4           1   166367755 rs12562373                     A                    G      0.7678     0.760  0.006007180 0.007904185 outcome
5          14    86737556  rs2135099                     A                    G      0.1652    -0.779 -0.007000891 0.008987023 outcome
6           2   201527977 rs57502521                     A                    G      0.9649    -0.088 -0.001595890 0.018135114 outcome
  mr_keep.outcome pval.outcome pval_origin.outcome id.outcome
1            TRUE    0.7802448            inferred     AIkFiR
2            TRUE    0.3821183            inferred     AIkFiR
3            TRUE    0.1320672            inferred     AIkFiR
4            TRUE    0.4472546            inferred     AIkFiR
5            TRUE    0.4359797            inferred     AIkFiR
6            TRUE    0.9298767            inferred     AIkFiR

```
## 为了在后面不报错，我魔改了一下这个函数
但是我忘记我原来对这个函数修改了什么了哈哈哈
> [!IMPORTANT]
> 注意这个是服务器版本的，windows版的fun2要修改一下
>  fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --clump ", gsub('\\\\',"/",fn), " --clump-p1 ", clump_p, " --clump-r2 ", 
                 clump_r2, " --clump-kb ", clump_kb, " --out ", gsub('\\\\',"/",fn))
  system(fun2)

``` R
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
```

## 循环执行函数

``` R
reslist <- list() # 获取基因列表
for (i in genelist) {
  exp_dat <- exp_dat1[exp_dat1$gene == i,]     
  epilepsy <- snp_epilepsy_out[snp_epilepsy_out$SNP %in% exp_dat$SNP, ]
  exp_dat <- exp_dat[exp_dat$SNP %in% epilepsy$SNP, ]

  ld <- ld_clump_local(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure),
                       clump_r2=0.001,clump_kb=10000,clump_p=1, 
                       bfile = "/GPUFS/gyfyy_jxhe_1/User/zyt/mr_gy/0127/1000genomes/EUR", 
                       plink_bin = plinkbinr::get_plink_exe())  #这个是天河服务器版本的，你要改自己的路径哦
  exp_dat <- exp_dat[exp_dat$SNP %in% ld$rsid, ]

if (nrow(exp_dat) == 0) {
  next
}       #这个是为了当没有数据时执行下一个循环

dat <- harmonise_data(exposure_dat=exp_dat, outcome_dat=epilepsy)
if (nrow(dat) == 0) {
  next
}
  res <- mr(dat)     #执行mr
  reslist[[i]] <- res   #把结果写入列表
  
}

result <- data.frame()   
for (i in names(reslist)) {
  if (nrow(reslist[[i]]) == 0) {
    next
  }    #当结果为空时跳到下一个循环
  r1 <- reslist[[i]]
  r1$gene  <-  i
  result <- rbind(result,r1)
}

```

## 我有一计可以多线程运行
数据名和上面的稍微有点区别，其实也就加了个format_data，应该很好理解

``` R
process_gene <- function(i) {
  cough1 <- cough[cough$Gene == i,]    #这个是eqtl文件
  cough1 <- cough1[cough1$p<5E-08]
  if (nrow(cough1) == 0) {
    return(NULL)
  }
  exp_dat <- format_data(
    dat=cough1,
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
  
  asthma2 <- disease[disease$rsid %in% exp_dat$SNP, ]     #这个是你的gwas文件
  if (nrow(asthma2) == 0) {
    return(NULL)
  }
  asthma <- format_data(
    dat=asthma2,
    type = "outcome",
    header = TRUE,
    snps = exp_dat$SNP,
    snp_col = "rsid",
    beta_col = "inv_var_meta_beta",
    se_col = "inv_var_meta_sebeta",
    eaf_col = "all_meta_AF",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    # samplesize_col = "all_meta_sample_N",
    pval_col = "inv_var_meta_p",
    chr_col = "CHR",
    pos_col = "POS"
  )
  
  exp_dat <- exp_dat[exp_dat$SNP %in% asthma$SNP, ]
  
  
  ld <- ld_clump_local(dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure),
                       clump_r2=0.001,clump_kb=10000,clump_p=1, bfile = "./GWAS_data/1000genomic/EUR", plink_bin = "./lin/projects/bin/plink")   #执行本地clump
  exp_dat <- exp_dat[exp_dat$SNP %in% ld$rsid, ]
  
  if (nrow(exp_dat) == 0) {
    return(NULL)
  }
  
  dat <- harmonise_data(exposure_dat=exp_dat, outcome_dat=asthma)
  if (nrow(dat) == 0) {
    return(NULL)
  }
  res <- mr(dat)
  
  return(res)
}

```

运行函数

``` R
library(parallel)
# 您可以根据您的机器的核心数调整mc.cores的值
reslist <- mclapply(genelist, process_gene, mc.cores = 8)

# 清理结果，移除NULL元素
reslist <- reslist[!sapply(reslist, is.null)]
```


