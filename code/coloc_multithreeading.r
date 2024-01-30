library(data.table)
library(dplyr)
library(coloc)
library(MungeSumstats)
library(parallel)



disease <- fread("cough/data/cough_37.txt")

id.e <- "disease"
disease$N <- 400000
expo <- disease
# expo<-dplyr::select(disease,SNP=rsid,CHR,BP=POS,A1=ALT,A2=REF,FRQ=all_meta_AF,BETA=inv_var_meta_beta,SE=inv_var_meta_sebeta,P=inv_var_meta_p,N)
expo <- expo%>%dplyr::mutate(Z=BETA/SE)
# expo <- liftover(expo,convert_ref_genome="hg19",ref_genome ="hg38",chrom_col = "CHR",
#                  start_col = "BP")
# expo <- expo[expo$IMPUTATION_gen_build == T,]
# expo <- expo[,1:11]
expo <- na.omit(expo, cols = "SNP")
expo <- as_tibble(expo)
expo <- na.omit(expo, cols = "P")
expo[expo$P == 0,]$P <- 1e-50


gene <- fread("drug/Whole_Blood.txt.gz")
genelist <- unique(gene$Gene)
id.o <- "gene"
outcome<-dplyr::select(gene,Gene,SNP,CHR=Chr,BP,A1,A2,BETA=b,SE=SE,P=p)
outcome$N <- 670
# outcome <- outcome%>%dplyr::mutate(Z=BETA/SE)
outcome <- na.omit(outcome, cols = "SNP")
outcome[outcome$P == 0,]$P <- 1e-50


coloc_gene <- function(i) {
  
  
  outcomegene <- outcome[outcome$Gene == i,]
  outcomegene <- outcomegene[!duplicated(outcomegene$SNP),]
  outcomegene <- as_tibble(outcomegene)
  
  expogene <- expo[expo$SNP %in% outcomegene$SNP,]
  outcomegene <- outcomegene[outcomegene$SNP %in% expogene$SNP ,]
  outcomegene <- outcomegene %>%
  left_join(expogene %>% select(SNP, FRQ), by = "SNP")
  # 如果 outcomegene 为空，则跳过此次循环
  if (nrow(outcomegene) == 0) {
    return(NULL)
  }
  
  
  
  datalist_names <- c("snp","position","beta","type","MAF","N","s","pvalues")
  eqtl_data_list <- vector("list", length(datalist_names))
  names(eqtl_data_list) <- datalist_names
  eqtl_data_list$snp <- expogene$SNP
  eqtl_data_list$position <- expogene$BP
  eqtl_data_list$beta <- expogene$BETA
  eqtl_data_list$pvalues <- expogene$P
  eqtl_data_list$MAF <- expogene$FRQ
  eqtl_data_list$s <- 0.07
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
  gwas_data_list$type <- "quant"
  gwas_data_list$N <- outcomegene$N

  
  abf_res <- coloc.abf(dataset1=eqtl_data_list,
                       dataset2=gwas_data_list)
  abf_res$gene <- i
  return(abf_res) 
}

reslist <- mclapply(genelist, coloc_gene, mc.cores = 4)

saveRDS(reslist,"drug/cough_coloc.rds")

result <- data.frame()
for (i in 1:length(reslist)) {
  result <- rbind(result,c(reslist[[i]]$summary,reslist[[i]]$gene))
}

colnames(result) <- c("nsnps",	"PP.H0.abf",	"PP.H1.abf",	"PP.H2.abf",	"PP.H3.abf",	"PP.H4.abf",	"gene")
saveRDS(result,"drug/cough_coloc_result.rds")

