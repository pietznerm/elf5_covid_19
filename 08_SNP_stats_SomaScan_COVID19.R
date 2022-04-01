#!/usr/bin/env Rscript

## script to obtain effect estimates for protein and COVID-19 
## Maik Pietzner 21/03/2022
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
## snp to query
crs   <- args[5]
## split
crs   <- strsplit(crs, ",")[[1]]

#-----------------------------------------#
##-- 	       import protein data       --##
#-----------------------------------------#

## read the relevant data
res.soma        <- paste0("zcat ~/rfs/rfs-epid-rfs-mpB3sSsgAn4//Studies//People/Ellie/SomaLogic/GWAS_discovery_subsets/Meta-Analysis/Fenland_auto_chrX_filtered/output/",
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
##--              store                --##
#-----------------------------------------#

## write data needed to file
write.table(subset(res.all, rsid %in% crs), paste("look_ups/lookup", soma, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)
