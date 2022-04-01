###############################################
#### ELF5 and Covid-19 and other outcomes  ####
#### Maik Pietzner              19/08/2021 ####
###############################################

rm(list=ls())
setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4/Studies/People/Maik/COVID19_release6/04_elf5/data/")
options(stringsAsFactors = F)
load(".RData")

##############################################
####   obtain the summary stats needed    ####
##############################################

## load most recent collation of results
res.pqtls <- read.table("../../../cis_pQTL_pipeline/input/SomaLogic.Targets.cis.regions.txt", sep="\t", header=T)

## aptamer of interest
soma      <- "res_invn_X13457_33"

## subset the results to aptamer of interest
res.pqtls <- subset(res.pqtls, pheno == soma) ## one cis

## define 1MB region to extract
chr.s     <- res.pqtls$chr
pos.m     <- res.pqtls$pos
pos.s     <- res.pqtls$region_start
pos.e     <- res.pqtls$region_end

## obtain stats
source("../scripts/get_somamer_stats.R")

## ELF5
res.elf5  <- res.region.soma(soma, chr.s, pos.s, pos.e)

## CAT (other gene highlighted)
res.cat   <- res.region.soma("res_invn_X3488_64", chr.s, pos.s, pos.e)

##############################################
####        SNPs for LD information       ####
##############################################

## run as as a bash job
system(paste("./../scripts/02_get_SNPs_LD", chr.s, pos.s, pos.e))

## read the dosage file
require(data.table)
tmp           <- fread(paste0("tmp.", chr.s, ".", pos.s, ".", pos.e, ".dosage"), sep=" ", header=T, data.table=F)
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
## reduce to SNPs in the stats
tmp.info            <- subset(tmp.info, MarkerName %in% res.elf5$MarkerName)

## rename
pheno <- tmp
rm(tmp); gc()

##############################################
####               cis-eQTL               ####
##############################################

## function to do so
source("../scripts/coloc_eQTL_GTEx_v8.R")

## ELF5
res.gtex.elf5 <- coloc.eQTL(res.elf5, chr.s, pos.e, pos.s, res.pqtls$ensembl_gene_id, pheno, tmp.info, r.t = T) 
## divide
stats.elf5    <- res.gtex.elf5[[2]]
res.gtex.elf5 <- res.gtex.elf5[[1]]

## CAT
res.gtex.cat  <- coloc.eQTL(res.cat, chr.s, pos.e, pos.s, "ENSG00000121691", pheno, tmp.info, r.t = T) 
## divide
stats.cat     <- res.gtex.cat[[2]]
res.gtex.cat  <- res.gtex.cat[[1]]

##############################################
####                PheWAS                ####
##############################################

## only for ELF5
source("../scripts/coloc_phewas.R")
res.phewas    <- phewas.coloc(res.elf5, chr.s, pos.s, pos.e, pheno, tmp.info)

##############################################
####            HyprColoc HGI             ####
##############################################

#--------------------------------#
##--     import HGI stats     --##
#--------------------------------#

## --> HGI data release : 01/06/2021 <-- ##
res.hgi <- lapply(c("A2", "B1", "B2", "C2"), function(x){
  ## get a data header first
  hd  <- paste0("zcat ../../01_stats_download/COVID19_HGI_", x, "_ALL_leave_23andme_20210607.b37.txt.gz", " | head -1")
  hd  <- unlist(data.table::fread(cmd=hd, sep="\t", header=F, data.table = F))
  ## read the data
  tmp <- paste0("zcat ../../01_stats_download/COVID19_HGI_", x, "_ALL_leave_23andme_20210607.b37.txt.gz",
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

## rename alleles
res.hgi[, c("Allele1", "Allele2")] <- t(apply(res.hgi[, c("Allele1", "Allele2")], 1, function(k){
  if(nchar(k[1]) > 1){
    return(c("I", "D"))
  }else if(nchar(k[2]) > 1){
    return(c("D", "I"))
  }else{
    return(k)
  }
}))

## now reshape
res.hgi   <- reshape(res.hgi, idvar=c("chr", "pos", "Allele1", "Allele2", "SNP"), timevar="outcome", direction="wide")

## create MarkerName to ease merging
res.hgi$MarkerName <- apply(res.hgi[, c("chr", "pos", "Allele1", "Allele2")], 1, function(x){
  return(paste0("chr", as.numeric(x[1]), ":", as.numeric(x[2]), "_", paste(sort(toupper(x[3:4])), collapse = "_")))
})

#---------------------------------------#
##--  create one molecular data set  --##
#---------------------------------------#

## keep only information really needed
res.all        <- stats.elf5[, c("MarkerName", "chr", "pos", "Allele1.gtex", "Allele2.gtex", "MAF", "slope.Lung", "slope_se.Lung", "pval_nominal.Lung", "N.Lung",
                                 "Effect.aligned", "StdErr", "Pvalue", "TotalSampleSize")]

## edit names
names(res.all) <- c("MarkerName", "chr", "pos", "Allele1", "Allele2", "MAF", "Effect.ELF5.Lung", "StdErr.ELF5.Lung", "Pvalue.ELF5.Lung", "TotalSampleSize.ELF5.Lung",
                    "Effect.ELF5.protein", "StdErr.ELF5.protein", "Pvalue.ELF5.protein", "TotalSampleSize.ELF5.protein")

## same for CAT
res.tmp        <- stats.cat[, c("MarkerName", "chr", "pos", "Allele1.gtex", "Allele2.gtex", "MAF", 
                                "slope.Lung", "slope_se.Lung", "pval_nominal.Lung", "N.Lung",
                                "slope.Whole_Blood", "slope_se.Whole_Blood", "pval_nominal.Whole_Blood", "N.Whole_Blood",
                                "Effect.aligned", "StdErr", "Pvalue", "TotalSampleSize")]

## edit names
names(res.tmp) <- c("MarkerName", "chr", "pos", "Allele1.gtex", "Allele2.gtex", "MAF", 
                    "Effect.CAT.Lung", "StdErr.CAT.Lung", "Pvalue.CAT.Lung", "TotalSampleSize.CAT.Lung",
                    "Effect.CAT.Whole_Blood", "StdErr.CAT.Whole_Blood", "Pvalue.CAT.Whole_Blood", "TotalSampleSize.Whole_Blood",
                    "Effect.CAT.protein", "StdErr.CAT.protein", "Pvalue.CAT.protein", "TotalSampleSize.CAT.protein")

## combine
res.all        <- merge(res.all, res.tmp)
## test whether effect slignment is needed
identical(res.all$Allele1, res.all$Allele1.gtex)
## drop not needed alleles
res.all$Allele1.gtex <- res.all$Allele2.gtex <- NULL

#---------------------------------------#
##--      combine with HGI data      --##
#---------------------------------------#

## add Covid stats
res.all <- merge(res.all, res.hgi, by=c("MarkerName", "chr", "pos"), suffix=c(".protein", ".covid"))

## align everything to the protein stats
for(j in c("A2", "B1", "B2", "C2")){
  res.all[, paste0("Effect.", j)] <- ifelse(res.all$Allele1.protein == res.all$Allele1.covid, res.all[, paste0("Effect.", j)], -res.all[, paste0("Effect.", j)])
}
## delete what is no longer needed
res.all$Allele1.covid <- res.all$Allele2.covid <- NULL

## drop NAs
res.all <- na.omit(res.all)

#----------------------------------------------#
##--              run hyprcoloc             --##
#----------------------------------------------#

## create label
lab.plot <- data.frame(name=gsub("Effect\\.", "", grep("^Effect\\.", names(res.all), value=T)), 
                       label=c("ELF5 transcript - Lung",
                               "ELF5 protein - Blood",
                               "CAT transcript - Lung",
                               "CAT transcript - Blood",
                               "CAT protein - Blood",
                               "Severe COVID-19 vs. population",
                               "Hospitalised/not hospitalised COVID-19",
                               "Hospitalised COVID-19 vs. population",
                               "COVID-19 vs. population"))

## add SNP id
res.all  <- merge(res.all, tmp.info[, c("MarkerName", "rsid")])

## import function to do so
source("../scripts/run_hyprcoloc.R")
pp.hyper <- run.hyprcoloc(res.all, lab.plot$name, c(0,0,0,0,0,1,1,1,1))

#----------------------------------------------#
##--            do some plotting            --##
#----------------------------------------------#


## load function to do so
source("../scripts/plot_hyprcoloc.R")

## ELF5
png("../graphics/ELF5.COVID19.20220318.png", width=8, height=10, res=900, units = "cm")
par(mar=c(.1,1.5,.75,.5), cex.axis=.5, cex.lab=.5, bty="l", tck=-.01, mgp=c(.6,0,0), xaxs="i", lwd=.5)
## add layout
layout(matrix(1:6,6,1), heights = c(rep(.3,5),.2))
## smaller region
plot.hyprcoloc(subset(res.all, pos >= pos.s+1e5 & pos <= pos.e-1e5), 
               subset(pp.hyper, p1 == 1e-4 & p2 == .98 & thr == .5), 2, lab.plot, pheno, tmp.info)
dev.off()

## CAT
png("../graphics/CAT.expression.20220318.png", width=8, height=6, res=900, units = "cm")
par(mar=c(.1,1.5,.75,.5), cex.axis=.5, cex.lab=.5, bty="l", tck=-.01, mgp=c(.6,0,0), xaxs="i", lwd=.5)
## add layout
layout(matrix(1:3,3,1), heights = c(rep(.3,3),.2))
plot.hyprcoloc(subset(res.all, pos >= pos.s+1e5 & pos <= pos.e-1e5), 
               subset(pp.hyper, p1 == 1e-4 & p2 == .98 & thr == .5), 
               1, lab.plot, pheno, tmp.info,
               a.vars = "rs766826")
dev.off()

##############################################
####          add lung function           ####
##############################################

## import conditional stats and plot (conditioned on chr11:35308988_C_T and chr11:35033491_A_G)
res.fcv            <- data.table::fread("FEV1_FVC.11.34000340.36000340.cma.cojo", sep="\t", header=T, data.table=F)
## reduce to information needed
res.fcv            <- res.fcv[, c("Chr", "SNP", "refA", "bC", "bC_se", "pC")]
## rename
names(res.fcv)     <- c("chr", "MarkerName", "refA", "Effect.FCV", "StdErr.FCV", "Pvalue.FCV")

## add to all the other traits
res.all            <- merge(res.all, res.fcv)
## align effect estimates
res.all$Effect.FCV <- ifelse(res.all$Allele1.protein == res.all$refA, res.all$Effect.FCV, -res.all$Effect.FCV)
## delete what is no longer needed
res.all$refA       <- NULL

## omit NAs
res.all            <- na.omit(res.all)

#----------------------------------------------#
##--              run hyprcoloc             --##
#----------------------------------------------#

## extend the labels
lab.plot      <- rbind(lab.plot, data.frame(name="FCV", label="Lung function (FEV1/FVC)"))

## run hyprcoloc
pp.hyper.lung <- run.hyprcoloc(res.all, lab.plot$name, c(0,0,0,0,0,1,1,1,1,0))

#----------------------------------------------#
##--            do some plotting            --##
#----------------------------------------------#

## ELF5
png("../graphics/ELF5.COVID19.Lung.20220318.png", width=8, height=11, res=900, units = "cm")
par(mar=c(.1,1.5,.75,.5), cex.axis=.5, cex.lab=.5, bty="l", tck=-.01, mgp=c(.6,0,0), xaxs="i", lwd=.5)
## add layout
layout(matrix(1:7,7,1), heights = c(rep(.3,6),.2))
plot.hyprcoloc(subset(res.all, pos >= pos.s+1e5 & pos <= pos.e-1e5), subset(pp.hyper.lung, p1 == 1e-4 & p2 == .98 & thr == .5), 2, lab.plot, pheno, tmp.info)
dev.off()
