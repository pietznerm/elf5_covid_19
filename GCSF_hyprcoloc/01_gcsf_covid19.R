###############################################
#### G-CSF and Covid-19 and other outcomes ####
#### Maik Pietzner              15/10/2021 ####
###############################################

rm(list=ls())
setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4/Studies/People/Maik/COVID19_release6/06_gcsf/data/")
options(stringsAsFactors = F)
load(".RData")

##############################################
####      obtain summary stats needed     ####
##############################################

## load regional meta-data
olink.targets              <- read.delim("../../../COVID19_v2/02_naive_cis_coloc/input/Olink.mapping.genomic.position.txt", sep="\t", header=T)
olink.targets$pheno        <- paste0("invn_res_X", tolower(olink.targets$varname)) 
## define regions
olink.targets$region_start <- olink.targets$start_position - 5e5
olink.targets$region_end   <- olink.targets$end_position   + 5e5
## edit start
olink.targets$region_start <- ifelse(olink.targets$region_start < 0, 1, olink.targets$region_start)

## aptamer of interest
olink                      <- "invn_res_Xg_csf_neu"

## subset the results to aptamer of interest
olink.targets              <- subset(olink.targets, pheno == olink) 

## define 1MB region to extract
chr.s                      <- olink.targets$chr
## use consistent region around the protein-encoding gene
pos.s                      <- olink.targets$region_start
pos.e                      <- olink.targets$region_end

#-----------------------------------------#
##-- 	       import protein data       --##
#-----------------------------------------#

## import Olink stats
res.olink                 <- paste0("zcat ~/rfs/rfs-epid-rfs-mpB3sSsgAn4//Studies/People/Ellie/Olink/Fenland_485/GWAS/output/formatted/combined/",
                                    olink, "_withMarkerName.out.gz",
                                    " | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                                    " '$1 == chr && $3 >= low && $3 <= upp {print $0}' -")
res.olink                 <- data.table::fread(cmd = res.olink, sep = " ", header = F, data.table = F)
## add names
names(res.olink)          <- c("chr", "rsid", "pos", "Allele2", "Allele1", "Freq1", "info", "Effect", "StdErr", "tval", "log10p", "MarkerName")
## apply MAF filter
res.olink$MAF             <- ifelse(res.olink$Freq1 > .5, 1 - res.olink$Freq1, res.olink$Freq1)
res.olink$TotalSampleSize <- 485
## compute p-value
res.olink$Pvalue          <- 10^-res.olink$log10p

#-----------------------------------------#
##--       import genotype data        --##
#-----------------------------------------#

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
tmp.info            <- subset(tmp.info, MarkerName %in% res.olink$MarkerName)

## rename
pheno <- tmp
rm(tmp); gc()

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
res.hgi            <- reshape(res.hgi, idvar=c("chr", "pos", "Allele1", "Allele2", "SNP"), timevar="outcome", direction="wide")

## create MarkerName to ease merging
res.hgi$MarkerName <- apply(res.hgi[, c("chr", "pos", "Allele1", "Allele2")], 1, function(x){
  return(paste0("chr", as.numeric(x[1]), ":", as.numeric(x[2]), "_", paste(sort(toupper(x[3:4])), collapse = "_")))
})

#--------------------------------#
##--    combine everything    --##
#--------------------------------#

res.all <- merge(res.olink, res.hgi, by=c("MarkerName", "chr", "pos"), suffix=c(".protein", ".covid"))

## align everything to the protein stats
for(j in c("A2", "B1", "B2", "C2")){
  res.all[, paste0("Effect.", j)] <- ifelse(res.all$Allele1.protein == res.all$Allele1.covid, res.all[, paste0("Effect.", j)], -res.all[, paste0("Effect.", j)])
}
## delete what is no longer needed
res.all$Allele1.covid <- res.all$Allele2.covid <- NULL

## drop NAs
res.all <- na.omit(res.all)

## rename
names(res.all)[c(9,10,15)] <- paste0(names(res.all)[c(9,10,15)], ".proteins")

#----------------------------------------------#
##--              run hyprcoloc             --##
#----------------------------------------------#

require(hyprcoloc)

## create matrices for input
betas           <- as.matrix(res.all[, grep("^Effect\\.", names(res.all), value=T)])
ses             <- as.matrix(res.all[, grep("^StdErr\\.", names(res.all), value=T)])
## add names
colnames(betas) <- colnames(ses) <- gsub("Effect\\.", "", colnames(betas))

## store SNP-IDs
snp.id          <- res.all$rsid

## run over a grid of priors
prior.grid      <- expand.grid(p1= 1e-4, p2 = c(0.95, 0.98, 0.99, 0.999), thr = c(0.5,0.6,0.7,0.8,0.9), stringsAsFactors = F)

## run hypercoloc
pp.hyper        <- lapply(1:nrow(prior.grid), function(x){
  print(prior.grid[x,])
  ## run hypcoloc with set of priors
  tmp <- hyprcoloc(betas, ses, trait.names = colnames(betas), snp.id=snp.id,
                   # bb.selection = "align", ## allow for multiple variants
                   ## define one binary outcome
                   binary.outcomes = c(0,1,1,1,1),
                   # trait.cor = cor_mat[colnames(betas),colnames(betas)],
                   prior.1 = prior.grid$p1[x],
                   prior.2 = prior.grid$p2[x],
                   reg.thresh = prior.grid$thr[x],
                   align.thresh = prior.grid$thr[x],
                   snpscores = T)
  # print(tmp)
  ## get results
  tmp     <- subset(tmp$results, traits != "None")
  if(nrow(tmp) > 0){
    ## add settings
    tmp$p1  <- prior.grid$p1[x]
    tmp$p2  <- prior.grid$p2[x]
    tmp$thr <- prior.grid$thr[x]
    return(tmp)
  }
})
## combine into one data frame
pp.hyper <- do.call(rbind, pp.hyper)
## protein and C2 most consistent

#----------------------------------------------#
##--            do some plotting            --##
#----------------------------------------------#

## create label
lab.plot       <- data.frame(name=colnames(betas), 
                             label=c("G-CSF Protein - Blood",
                                     "Severe COVID-19 vs. population",
                                     "Hospitalised/not hospitalised COVID-19",
                                     "Hospitalised COVID-19 vs. population",
                                     "COVID-19 vs. population"))

## load function to do so
source("../scripts/plot_hyprcoloc.R")
png("../graphics/CSF.COVID19.png", width=8, height=8, res=300, units = "cm")
par(mar=c(.1,1.5,.75,.5), cex.axis=.5, cex.lab=.5, bty="l", tck=-.01, mgp=c(.6,0,0), xaxs="i", lwd=.5)
## add layout
layout(matrix(1:4,4,1), heights = c(rep(.3,3),.2))
plot.hyprcoloc(subset(res.all, pos >= pos.s+1e5 & pos <= pos.e-1e5), subset(pp.hyper, p1 == 1e-4 & p2 == .98 & thr == .5), 1, lab.plot, pheno, tmp.info)
dev.off()

##############################################
####            look into PheWAS          ####
##############################################

## import function to do so (default priors!)
source("../scripts/coloc_phewas.R")
res.csf <- phewas.coloc(res.olink, chr.s, pos.s, pos.e, pheno, tmp.info)
subset(res.csf, PP.H4.abf > .8 & ld.check.top > .8)
save.image()

## write to file for supplemental table
write.table(subset(res.csf, PP.H4.abf > .8 & ld.check.top > .8 & nsnps > 500), "Results.coloc.G.CSF.Olink.20220321.txt", sep="\t", row.names=F)

##############################################
####            look into GTEx            ####
##############################################

## import function to do so
source("../scripts/coloc_eQTL_GTEx_v8.R")

## CSF
res.gtex.csf <- coloc.eQTL(res.olink, chr.s, pos.e, pos.s, olink.targets$ensembl_gene_id, pheno, tmp.info, r.t = T) 
## divide
stats.csf    <- res.gtex.csf[[2]]
res.gtex.csf <- res.gtex.csf[[1]]
## seems to be no signal

##############################################
####    import core blood cell traits     ####
##############################################

## get subset of singnificant traits
sig.phe      <- subset(res.csf, PP.H4.abf > .8 & ld.check.top > .8 & nsnps > 500)

#------------------------------------------------#
##-- blood cell traits from Astle et al. 2016 --##
#------------------------------------------------#

## load package to query the database
require(ieugwasr)

## now get those
res.blood <- lapply(subset(sig.phe, nchar(id.ieu) == 16)$id.ieu, function(x){
  
  ## obtain summary stats and trait info
  reg        <- paste0(chr.s, ":", pos.s,"-", pos.e)
  res        <- associations(reg, x)
  ## just to make sure everything will work fine
  res        <- subset(res, !is.na(beta))
  ## make unique, since some data processing error
  res        <- unique(res)
  
  ## get meta information
  return(res)
  
})
## combine into one data frame
res.blood <- do.call(rbind, res.blood)

## reshape
library(tidyr)
res.blood <- pivot_wider(res.blood, id_cols = c("chr", "position", "rsid", "ea", "nea"),
                               values_from = c("beta", "se", "p", "n"),
                               names_from = "id", names_sep = ".")
## convert to data frame
res.blood <- data.frame(res.blood)

## drop ambiguous IDs
ii              <- table(res.blood$rsid)
ii              <- names(ii[ii>1]) 
if(length(ii) > 0){
  res.blood <- subset(res.blood, !(rsid %in% ii))
}

#------------------------------------------------#
##--           combine to do hyprcoloc        --##
#------------------------------------------------#

## dom some renaming
names(res.blood) <- gsub("beta", "Effect", names(res.blood))
names(res.blood) <- gsub("^se", "StdErr", names(res.blood))
names(res.blood) <- gsub("^p\\.", "Pvalue.", names(res.blood))
names(res.blood) <- gsub("^n\\.", "N.", names(res.blood))
## drop rsID
res.blood$rsid   <- NULL

## combine
res.all          <- merge(res.all, res.blood, by.x=c("chr", "pos"), by.y=c("chr", "position"))

#----------------------------------------------#
##--              run hyprcoloc             --##
#----------------------------------------------#

require(hyprcoloc)

## create matrices for input
betas           <- as.matrix(res.all[, grep("^Effect\\.", names(res.all), value=T)])
ses             <- as.matrix(res.all[, grep("^StdErr\\.", names(res.all), value=T)])
## add names
colnames(betas) <- colnames(ses) <- gsub("Effect\\.", "", colnames(betas))

## store SNP-IDs
snp.id          <- res.all$rsid

## run over a grid of priors
prior.grid      <- expand.grid(p1= 1e-4, p2 = c(0.95, 0.98, 0.99, 0.999), thr = c(0.5,0.6,0.7,0.8,0.9), stringsAsFactors = F)

## run hypercoloc
pp.hyper.blood  <- lapply(1:nrow(prior.grid), function(x){
  print(prior.grid[x,])
  ## run hypcoloc with set of priors
  tmp <- hyprcoloc(betas, ses, trait.names = colnames(betas), snp.id=snp.id,
                   # bb.selection = "align", ## allow for multiple variants
                   ## define one binary outcome
                   binary.outcomes = c(0,1,1,1,1, rep(0,9)),
                   # trait.cor = cor_mat[colnames(betas),colnames(betas)],
                   prior.1 = prior.grid$p1[x],
                   prior.2 = prior.grid$p2[x],
                   reg.thresh = prior.grid$thr[x],
                   align.thresh = prior.grid$thr[x],
                   snpscores = T)
  # print(tmp)
  ## get results
  tmp     <- subset(tmp$results, traits != "None")
  if(nrow(tmp) > 0){
    ## add settings
    tmp$p1  <- prior.grid$p1[x]
    tmp$p2  <- prior.grid$p2[x]
    tmp$thr <- prior.grid$thr[x]
    return(tmp)
  }
})
## combine into one data frame
pp.hyper.blood <- do.call(rbind, pp.hyper.blood)

#----------------------------------------------#
##--            do some plotting            --##
#----------------------------------------------#

## add blood trait label
tmp            <- subset(sig.phe, nchar(id.ieu) == 16)[, c("id.ieu", "trait")]
names(tmp)     <- c("name", "label")
## edit to match to label
tmp$name       <- gsub("-", ".", tmp$name)

## create label
lab.plot       <- rbind(lab.plot, tmp)

## --> stacked regional association plot <-- ##

## load function to do so
source("../scripts/plot_hyprcoloc.R")
png("../graphics/CSF.COVID19.blood.cell.traits.20220321.png", width=8, height=16, res=300, units = "cm")
par(mar=c(.1,1.5,.75,.5), cex.axis=.5, cex.lab=.5, bty="l", tck=-.01, mgp=c(.6,0,0), xaxs="i", lwd=.5)
## add layout
layout(matrix(1:13,13,1), heights = c(rep(.3,12),.3))
plot.hyprcoloc(subset(res.all, pos >= pos.s+1e5 & pos <= pos.e-1e5), subset(pp.hyper.blood, p1 == 1e-4 & p2 == .98 & thr == .5), 1, lab.plot, pheno, tmp.info)
dev.off()

#-------------------------------------#
##-- revised version for Manuscript --#
#-------------------------------------#

## create temporary file for plotting
tmp        <- subset(pp.hyper.blood, p1 == 1e-4 & p2 == .98 & thr == .5)
tmp$traits <- "proteins, ebi.a.GCST004610, ebi.a.GCST004614, ebi.a.GCST004629, A2, C2"

source("../scripts/plot_hyprcoloc.R")
png("../graphics/CSF.COVID19.blood.cell.traits.manuscript.20220321.png", width=8, height=12, res=300, units = "cm")
par(mar=c(.1,1.5,.75,.5), cex.axis=.5, cex.lab=.5, bty="l", tck=-.01, mgp=c(.6,0,0), xaxs="i", lwd=.5)
## add layout
layout(matrix(1:7,7,1), heights = c(rep(.3,6),.3))
plot.hyprcoloc(subset(res.all, pos >= pos.s+2e5 & pos <= pos.e-2e5), tmp, 1, lab.plot, pheno, tmp.info)
dev.off()
