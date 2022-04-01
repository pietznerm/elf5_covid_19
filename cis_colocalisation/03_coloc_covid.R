#!/usr/bin/env Rscript

## script to run coloc for SOMAmer and COVID GWAS
## Maik Pietzner 15/07/2021
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4//Studies//People/Maik/COVID19_release6/02_naive_cis_coloc/")

## --> import parameters <-- ##

## aptamer
soma  <- args[1]
## chromosome
chr.s <- as.numeric(args[2])
## start position of the region
pos.s <- as.numeric(args[3])
## end position of the region
pos.e <- as.numeric(args[4])

#-----------------------------------------#
##-- 	       import protein data       --##
#-----------------------------------------#

## read the relevant data
res.soma        <- paste0("zcat ~/rds/rds-rjh234-mrc-epid/Studies/People/Ellie/SomaLogic/GWAS_discovery_subsets/Meta-Analysis/Fenland_auto_chrX_filtered/output/",
                          soma,"_Fenland_MA_auto_chrX_filtered.txt.gz",
                          " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                          " '$18 == chr && $19 >= low && $19 <= upp {print $0}' -")
res.soma        <- data.table::fread(cmd = res.soma, sep = "\t", header = F, data.table = F)
## add names
names(res.soma) <- c("rsid", "MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "StdErr", "Pvalue", "Direction", "HetSq", "HetChiSq",
                     "HetDf", "HetPVal", "TotalSampleSize", "chr", "pos")
## keep only needed
res.soma        <- res.soma[, c("rsid", "MarkerName", "Freq1", "Allele1", "Allele2", "Effect", "StdErr", "Pvalue", "TotalSampleSize", "chr", "pos")]

## apply MAF filter
res.soma$MAF <- ifelse(res.soma$Freq1 > .5, 1 - res.soma$Freq1, res.soma$Freq1)
# res.soma     <- subset(res.soma, MAF >= .01)

#-----------------------------------------#
##-- 	        import trait data        --##
#-----------------------------------------#

## --> HGI data release : 01/06/2021 <-- ##
res.hgi <- lapply(c("A2", "B1", "B2", "C2"), function(x){
  ## get a data header first
  hd  <- paste0("zcat ../01_stats_download/COVID19_HGI_", x, "_ALL_leave_23andme_20210607.b37.txt.gz", " | head -1")
  hd  <- unlist(data.table::fread(cmd=hd, sep="\t", header=F, data.table = F))
  ## read the data
  tmp <- paste0("zcat ../01_stats_download/COVID19_HGI_", x, "_ALL_leave_23andme_20210607.b37.txt.gz",
                " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                " '$1 == chr && $2 >= low && $2 <= upp {print $0}' -")
  print(tmp)
  tmp         <- data.table::fread(cmd=tmp, sep="\t", header=F, data.table = F)
  if(nrow(tmp) > 0){
    ## rename first entry
    hd[1]       <- "chr"
    ## add header
    names(tmp)  <- hd
    ## reduce to what is needed
    tmp         <- tmp[, c("chr", "POS", "REF", "ALT", "SNP", "all_meta_N", "all_inv_var_meta_beta", "all_inv_var_meta_sebeta", "all_inv_var_meta_p", 
                           "all_inv_var_het_p", "all_inv_var_meta_effective", "all_meta_AF", "all_inv_var_meta_cases", "all_inv_var_meta_controls", "rsid")]
    # ## add names
    names(tmp)  <- c("chr", "pos", "Allele2", "Allele1", "SNP", "N_studies", "Effect", "StdErr", "Pvalue", "het_p", "TotalSampleSize", "Freq1", "n.case", "n.controls", "rsid")
    ## add phenotype
    tmp$outcome <- x
    print(summary(tmp))
    
    ## possibly drop SNPs with signficantly less observations
    
    return(tmp)
  }
})
## combine and reshape
res.hgi      <- do.call(rbind, res.hgi)

## drop rsid
res.hgi$rsid <- NULL

## store minimum p-value
p.min        <- min(res.hgi$Pvalue, na.rm=T)

#-----------------------------------------#
##-- 	            run coloc            --##
#-----------------------------------------#

require(coloc)

## proceed only if at least some evidence in the region
if(p.min < 1e-5){
  
  #-----------------------------------------#
  ##--           prepare data            --##
  #-----------------------------------------#
  
  ## rename alleles
  res.hgi[, c("Allele1", "Allele2")] <- t(apply(res.hgi[, c("Allele1", "Allele2")], 1, function(k){
    if(nchar(k[1]) > 1){
      return(c("i", "d"))
    }else if(nchar(k[2]) > 1){
      return(c("d", "i"))
    }else{
      return(tolower(k))
    }
  }))
  
  ## now reshape
  res.hgi   <- reshape(res.hgi, idvar=c("chr", "pos", "Allele1", "Allele2", "SNP"), timevar="outcome", direction="wide")
  
  ## combine with results from SomaLogic (careful: INDELs are coded differently)
  res.all   <- merge(res.soma, res.hgi, by=c("chr", "pos"), suffixes = c(".soma", ".covid"))
  
  #-----------------------------------------#
  ##--          import LD-matrix         --##
  #-----------------------------------------#
  
  ## read the dosage file
  require(data.table)
  tmp           <- fread(paste0("tmp_input//tmp.", soma, ".", chr.s, ".", pos.s, ".", pos.e, ".dosage"), sep=" ", header=T, data.table=F)
  ## transpose
  rownames(tmp) <- tmp$rsid
  ## store allele information to realign effect estimates afterwards
  tmp.info      <- tmp[, 1:6]
  tmp           <- t(tmp[,-c(1:6)])
  ## retransform to data set (keep X in mind for variable names)
  tmp           <- data.frame(ID_1=rownames(tmp), tmp)
  
  ## create another column to info to map names from the SNP data set
  tmp.info$id         <- sapply(tmp.info$rsid, function(x) ifelse(substr(x, 1, 2) == "rs", x, paste0("X", gsub(":", ".", x))))
  ## edit some IDs (X-chromosome)
  tmp.info$id         <- gsub("XX", "X", tmp.info$id)
  tmp.info$id         <- gsub("XAffx-", "Affx.", tmp.info$id)
  ## create MarkerName column as well
  tmp.info$MarkerName <- apply(tmp.info, 1, function(x){
    paste0("chr", as.numeric(x[1]), ":", as.numeric(x[4]), "_", paste(sort(x[5:6]), collapse = "_"))
  })
  
  ## rename
  pheno <- tmp
  
  #-----------------------------------------#
  ##--            run coloc              --##
  #-----------------------------------------#
  
  ## add to ease mapping of LD matrix
  res.all   <- merge(res.all, tmp.info[, c("MarkerName", "id")], by="MarkerName")
  
  ## order by position
  res.all   <- res.all[order(res.all$pos),]
  
  ## now run coloc for all examples
  res.coloc <- lapply(c("A2", "B1", "B2", "C2"), function(x){
    
    ## get the needed columns
    tmp <- na.omit(res.all[, c("rsid", "Effect", "StdErr", "TotalSampleSize", "MAF", paste(c("Effect", "StdErr", "n.case", "n.controls"), x, sep="."))])
    print(tmp[which.max(abs(tmp$Effect/tmp$StdErr)),])
    
    ## compute s by maximum
    s   <- max(tmp[, paste0("n.case.",x)])/(max(tmp[, paste0("n.case.",x)]) + max(tmp[, paste0("n.controls.",x)]))
    print(s)
    
    ## restrict to SNPs with at least 80% cases
    ct  <- max(tmp[, paste0("n.case.",x)])*.8
    tmp <- tmp[which(tmp[, paste0("n.case.",x)] >= ct), ]
    print(dim(tmp))
    
    ## run coloc
    res <- coloc.abf(dataset1=list(beta=tmp[, "Effect"], varbeta=tmp[, "StdErr"]^2, N=max(tmp[, "TotalSampleSize"]), type="quant", sdY=1),
                     ## obtained numbers from different cohorts
                     dataset2=list(beta=tmp[, paste0("Effect.",x)], varbeta=tmp[, paste0("StdErr.",x)]^2, 
                                   N=max(tmp[, paste0("n.case.",x)] + tmp[, paste0("n.controls.",x)]), 
                                   type="cc", 
                                   s=s),
                     MAF=tmp[, "MAF"])
    
    ## store information
    res              <- data.frame(t(res$summary))
    res$pheno        <- soma
    res$outcome      <- x
    res$chr          <- chr.s
    res$region_start <- pos.s
    res$region_end   <- pos.e
    
    ## add strongest hits for each
    ii  <- tmp$rsid[which.max(abs(tmp$Effect/tmp$StdErr))]
    jj  <- tmp$rsid[which.max(abs(tmp[, paste0("Effect.",x)]/tmp[, paste0("StdErr.",x)]))]
    
    ## add to the data
    foo <- subset(res.all, rsid == ii)[, c("MarkerName", "rsid", "Allele1.soma", "Allele2.soma", "Effect", "StdErr", "Pvalue", 
                                                         "Allele1.covid", "Allele2.covid", paste(c("Effect", "StdErr", "Pvalue"), x, sep=".")), drop=F]
    ## create output
    tmp        <- cbind(res, foo[1,])
    tmp$signal <- 1:nrow(tmp)
    
    ## add LD-check
    tmp$ld.check <- cor(pheno[,tmp.info$id[which(tmp.info$rsid == ii)]], pheno[,tmp.info$id[which(tmp.info$rsid == jj)]])^2
    ## check whether protein signal is preserved in the overlap
    tmp$ld.sens  <- cor(pheno[,tmp.info$id[which(tmp.info$rsid == ii)]], pheno[,tmp.info$id[which(tmp.info$rsid == res.soma$rsid[which.max(abs(res.soma$Effect/res.soma$StdErr))])]])^2
    
    ## edit some names
    names(tmp) <- gsub(x, "outcome", names(tmp))
    
    return(tmp)
    
  })
  res.coloc <- do.call(rbind, res.coloc)
  
  ## write to file
  write.table(res.coloc, paste("output/results.coloc", soma, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
  
  #-----------------------------------------#
  ##--           do some plotting        --##
  #-----------------------------------------#
  
  ## restrict to informative once
  tmp <- subset(res.coloc, PP.H4.abf > .25 | ld.check > .8)
  tmp <- unique(tmp[, 1:8])
  
  ## start plotting
  if(nrow(tmp) > 0){
    for(j in 1:nrow(tmp)){
      source("scripts/plot_coloc_results.R")
      png(paste0("graphics/", soma, ".", tmp[j,8], ".", chr.s, ".", pos.s, ".", pos.e, ".png"), width=16, height=8, units="cm", res=200)
      par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, cex.main=.6, font.main=2)
      ## more complex layout for gene assignment
      layout(matrix(c(1,1,1,2,3,4),3,2), heights = c(.43,.37,.2))
      plot.regional.coloc(res.all, tmp[j,], pheno, tmp[j,8], soma, tmp.info)
      dev.off()
    }
  }
}else{
  cat("found no evidence for", soma, "on chr", chr.s, "between", pos.s, "and", pos.e, "\n")
  ## write to file
  write.table(data.frame(pheno=soma, outcome=NA), paste("output/results.coloc", soma, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
}

