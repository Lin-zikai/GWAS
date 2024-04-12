# This one is mainly about how to do phewas

> [!IMPORTANT]
> It is recommended to run on a server

## 数据来源
我个人比较喜欢的phewas数据库是https://pheweb.org/UKB-TOPMed/ ，我将使用这个进行举例

> [!NOTE]
CKB也有一个https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10435379/，不过我还没拿到权限，老是忘记去问

首先先下载数据，有个大几百份数据的，具体的数据可以看 https://github.com/Lin-zikai/GWAS/blob/main/data/phenotypes%20.tsv 或者自己在网站中下载也行


我介绍一下我的大致思路：首先先确定我需要的暴露的snp和数据，根据这些snp，使用shell脚本去获取结局数据，然后合并成一份文件，读进r进行mr分析。我认为这个方案显著压缩了运算时间，如果有更好的方案还请各位指教

#现在正式开始进行phewas了
##1.首先你需要准备两份文件

第一份文件格式如下，我将其命名为expdat.csv,最后一列是蛋白的名字，可以放多个蛋白
``` 
SNP	chr.exposure	pos.exposure	effect_allele.exposure	other_allele.exposure	eaf.exposure	beta.exposure	se.exposure	pval.exposure	exposure	mr_keep.exposure	pval_origin.exposure	id.exposure
rs148574467	2	238511920	A	G	NA	-0.357779	0.10503	0.000658183	exposure	TRUE	reported	HES6
rs6749854	2	239113126	A	G	NA	-0.126647	0.0356435	0.000380635	exposure	TRUE	reported	HES6
rs62194936	2	239174693	C	T	NA	-0.138153	0.0380481	0.000282311	exposure	TRUE	reported	HES6
rs3791454	2	240049647	C	T	NA	0.104537	0.0270646	0.000112236	exposure	TRUE	reported	HES6
rs7672991	4	41910652	T	C	NA	-0.283449	0.0800264	0.000397187	exposure	TRUE	reported	LIMCH1
rs8866	17	65373979	G	C	NA	-0.120051	0.0223952	8.30E-08	exposure	TRUE	reported	PITPNC1
```
第二份文件是只有snp的一个txt，我将其命名为snp.txt,格式如下
> [!IMPORTANT]
> 这个没有行名的

```
rs148574467
rs6749854
rs62194936
rs3791454
rs7672991
rs8866
```

不知道怎么制作这两份文件话可以使用下面这个代码
现在你应该是已经获得了你需要分析的蛋白名称，你只需要把名字放到genelist中

> [!NOTE]
看不懂下面这个代码的话，可以看https://github.com/Lin-zikai/GWAS/blob/main/local_clump.md ，有详细介绍，其实就是做了一个clump而已

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
snplist <- exp[,"SNP"]
fwrite(snplist,"snp.txt",col.names = F)
```

##2.获取结局数据
现在根据前面我们确定需要的snp去结局文件中去获得这些snp

> [!CAUTION]
> 我是直接在装有这几百份文件的文件夹中直接写的sh脚本和运行，如果你存放位置不同的话需要你自己改路径哦

我将这个文件命名为extract_rs2.sh，和那些文件保存在一起,-P 8默认8个线程，自己调整
``` shell
# SNP 列表文件
SNP_LIST="snp.txt"

# 输出文件
OUTPUT_FILE="extracted_rows.tsv"

# 清空输出文件，确保开始时文件是空的
> "$OUTPUT_FILE"

# 读取 SNP 列表并构建正则表达式
SNP_PATTERN=$(awk '{printf "\\b"$0"\\b|"}' $SNP_LIST | sed 's/|$//')

# 使用 xargs 和 zgrep 在多线程中处理文件
# -P 参数指定并行进程的数量
# -I{} 用于替换每个输入项
find . -name "*.tsv.gz" | xargs -P 8 -I{} sh -c "zgrep -E '$SNP_PATTERN' '{}' | awk -v fname='{}' '{print \$0 \"\t\" fname}' >> '$OUTPUT_FILE'"

echo "Extraction complete. Results are in $OUTPUT_FILE"

```

## 当你一次性输入很多个基因的时候，那么find . -name "*.tsv.gz" | xargs可能会报错以下内容
``` R
/GPUFS/gyfyy_jxhe_1/User/zyt/Phewas/drug/extract_rs2.sh: line 18: /usr/bin/xargs: Argument list too long
Warning message:
In fread("extracted_rows.tsv") :
  File 'extracted_rows.tsv' has size 0. Returning a NULL data.table.
```

## 这时候可以每100个处理一次，extract_rs2.sh换成以下内容（-P 8默认16个线程，自己调整）

``` R
#!/bin/bash

SNP_LIST="snp.txt"

OUTPUT_FILE="extracted_rows.tsv"

> "$OUTPUT_FILE"

TOTAL_SNPS=$(wc -l < "$SNP_LIST")
BATCH_SIZE=100
NUM_BATCHES=$(( (TOTAL_SNPS + BATCH_SIZE - 1) / BATCH_SIZE ))

for (( i=0; i<NUM_BATCHES; i++ )); do
    START_LINE=$(( i * BATCH_SIZE + 1 ))
    END_LINE=$(( (i + 1) * BATCH_SIZE ))
    
    SNP_PATTERN=$(sed -n "${START_LINE},${END_LINE}p" $SNP_LIST | awk '{printf "\\b"$0"\\b|"}' | sed 's/|$//')
    
    find . -name "*.tsv.gz" | xargs -P 16 -I{} sh -c "zgrep -E '$SNP_PATTERN' '{}' | awk -v fname='{}' -v OFS='\t' '{print \$0, fname}' >> '$OUTPUT_FILE'"
    
    echo "Batch $((i+1)) of $NUM_BATCHES completed."
done

echo "Extraction complete. Results are in $OUTPUT_FILE"
```

使用R调用shell，可以把这几个R代码合并在一起就可以执行一个R脚本完成所有工作，只需要修改genelist
``` R
result <- system("/GPUFS/gyfyy_jxhe_1/User/lzk/drug/extract_rs2.sh", intern = TRUE)
write.csv(result,"result.csv")
```

##3.执行MR分析

``` R

exp_dat <- fread("expdat.csv")
out <- fread("extracted_rows.tsv")
sigene <- unique(exp_dat$id.exposure)
all <- data.table()
for (g in sigene) {
  exp2 <- exp_dat[exp_dat$id.exposure == g,]
  
  
  colnames(out) <- c("chrom"   ,      "pos"     ,      "ref"     ,      "alt"      ,     "rsids"     ,    "nearest_genes", "consequence"  ,
                     "pval"     ,     "beta"     ,     "sebeta" ,       "af"     ,       "case_af"   ,    "control_af"  ,  "tstat" ,"name")
  outlist <- unique(out$name) 
  out2 <- out[out$rsids == exp2$SNP]
  if (nrow(out2) == 0) {
    next
  }
  process_gene <- function(i) {
    outsnp <- out2[out2$name == i]
    if (nrow(outsnp) == 0) {
      return(NULL)
    }
    outsnp <- format_data(
      dat=outsnp,
      type = "outcome",
      header = TRUE,
      snps = exp2$SNP,
      snp_col = "rsids",
      beta_col = "beta",
      se_col = "sebeta",
      eaf_col = "af",
      effect_allele_col = "alt",
      other_allele_col = "ref",
      # samplesize_col = "all_meta_sample_N",
      pval_col = "pval",
      chr_col = "chrom",
      pos_col = "pos"
    )
    outsnp$id.outcome <- i
    if (nrow(outsnp) == 0) {
      return(NULL)
    }
    dat <- harmonise_data(exposure_dat=exp2, outcome_dat=outsnp,action = 1)
    if (nrow(dat) == 0) {
      return(NULL)
    }
    res <- tryCatch({
        mr(dat) 
    }, error = function(e) {
        cat("Error in MR process for gene", i, ": ", e$message, "\n")
        return(NULL)
    })
    
    if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) {
        return(NULL)
    }
    
    return(res)
}
  
  
  reslist <- mclapply(outlist, process_gene, mc.cores = 8)
  
  saveRDS(reslist,"mr_res.rds")
  
  result <- data.frame()
  for (i in 1:length(reslist)) {
    if (is.null(nrow(reslist[[i]])) == T) {
      next
    }
    r1 <- reslist[[i]]
    r1$gene  <-  i
    result <- rbind(result,r1)
  }
  
  all <- rbind(all,result)
}


a <- fread("phenotypes.tsv")
result <- all 
head(a)
head(result)
result$phenocode<-substring(result$id.outcome,3)
result$phenocode <- gsub("\\.tsv\\.gz", "", result$phenocode)
# 假设您的数据框叫做 df
integer_part <- as.integer(a$phenocode)
decimal_part <- sub("^[^.]*\\.?", "", a$phenocode)

# 格式化整数部分
formatted_integer_part <- sprintf("%03d", integer_part)

# 重新组合整数和小数部分
a$phenocode <- ifelse(nchar(decimal_part) > 0, 
                      paste0(formatted_integer_part, ".", decimal_part), 
                      formatted_integer_part)
a <- a[,c(11,13,15)]
all <- left_join(result,a,by="phenocode")

# 计算 -log10(p-value)
all$minus_log10_pval <- -log10(all$pval)

write.csv(all,"result.csv")
```

result.csv文件里就是所有的phewas结果啦



再送你一个曼哈顿图
``` R
# 绘制曼哈顿图
ggplot(all, aes(x = reorder(category, minus_log10_pval), y = minus_log10_pval, color = category)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = 'Category', y = '-log10(p-value)', title = 'Manhattan Plot by Category') +
  scale_color_hue(l = 50) # 使用不同颜色代表不同类别
```



