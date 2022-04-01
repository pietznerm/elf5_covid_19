#!/usr/bin/env Rscript

## script to obtain effect estimates for Olink proteins and COVID-19 
## Maik Pietzner 22/03/2022
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)

setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4//Studies//People/Maik/COVID19_release6/02_naive_cis_coloc/")

## --> import parameters <-- ##

## aptamer
olink  <- args[1]
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
if(chr.s == 23){
  res.olink        <- paste0("zcat ~/rfs/rfs-epid-rfs-mpB3sSsgAn4//Studies/People/Ellie/Olink/Fenland_485/GWAS/output/formatted/combined_chrX/",
                             olink,"_chrX_forMA.gz",
                             " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                             " '$1 == chr && $3 >= low && $3 <= upp {print $0}' -")
  res.olink        <- data.table::fread(cmd = res.olink, sep = "\t", header = F, data.table = F)
  
}else{
  res.olink        <- paste0("zcat ~/rfs/rfs-epid-rfs-mpB3sSsgAn4//Studies/People/Ellie/Olink/Fenland_485/GWAS/output/formatted/combined/",
                             olink, "_withMarkerName.out.gz",
                             " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                             " '$1 == chr && $3 >= low && $3 <= upp {print $0}' -")
  res.olink        <- data.table::fread(cmd = res.olink, sep = " ", header = F, data.table = F)
  
}
## add names
names(res.olink)          <- c("chr", "rsid", "pos", "Allele2", "Allele1", "Freq1", "info", "Effect", "StdErr", "tval", "log10p", "MarkerName")
## apply MAF filter
res.olink$MAF             <- ifelse(res.olink$Freq1 > .5, 1 - res.olink$Freq1, res.olink$Freq1)
# res.olink     <- subset(res.olink, MAF >= .01)
res.olink$TotalSampleSize <- 485
## compute p-value
res.olink$Pvalue          <- 10^-res.olink$log10p

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

## combine with results from olinkLogic (careful: INDELs are coded differently)
res.all   <- merge(res.olink, res.hgi, by=c("chr", "pos"), suffixes = c(".olink", ".covid"))

#-----------------------------------------#
##--              store                --##
#-----------------------------------------#

## write data needed to file
write.table(subset(res.all, rsid %in% crs), paste("look_ups/lookup", olink, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names=F)

