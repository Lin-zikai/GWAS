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
reslist <- list() 
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
}

dat <- harmonise_data(exposure_dat=exp_dat, outcome_dat=epilepsy)
if (nrow(dat) == 0) {
  next
}
  res <- mr(dat)
  reslist[[i]] <- res
  
}

saveRDS(reslist,"reslist.rds")




result <- data.frame()
for (i in names(reslist)) {
  if (nrow(reslist[[i]]) == 0) {
    next
  }
  r1 <- reslist[[i]]
  r1$gene  <-  i
  result <- rbind(result,r1)
}
saveRDS(result,"result1109.rds")

```
累了，明天再写
