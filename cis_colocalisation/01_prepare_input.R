################################################
#### cis-pQTL coloc COVID-19 GWAS release 6 ####
#### Maik Pietzner               15/07/2021 ####
################################################

rm(list=ls())
setwd("~/rfs/rfs-epid-rfs-mpB3sSsgAn4//Studies/People/Maik/COVID19_release6/02_naive_cis_coloc/input/")
options(stringsAsFactors = F)
load(".RData")

##############################################
####        load SomaLogic targets        ####
##############################################

## import cis-regions from SomaScan assay
regions.soma <- read.table("/home/mdp50/rds/rds-rjh234-mrc-epid/Studies/People/Maik/cis_pQTL_pipeline/input/SomaLogic.Targets.cis.regions.txt", sep="\t", header=T) 

##############################################
####        create input for coloc        ####
##############################################

## subet to promising regions
coloc.regions        <- regions.soma
coloc.regions$Pvalue <- as.numeric(coloc.regions$Pvalue) 
coloc.regions        <- subset(coloc.regions, Pvalue < 1e-5)

## create input file
write.table(coloc.regions[, c("pheno", "chr", "region_start", "region_end")], "coloc.regions.txt", sep="\t", row.names=F, col.names=F, quote=F)

##############################################
####            import results            ####
##############################################

## collect output
ii                 <- dir("../output/")

## import coloc results
res.coloc <- lapply(ii, function(x){
  ## read results
  tmp              <- read.table(paste0("../output/", x), sep="\t", header=T)
  ## add some identifier
  tmp$region_start <- as.numeric(strsplit(x, "\\.")[[1]][5]) 
  tmp$region_end   <- as.numeric(strsplit(x, "\\.")[[1]][6])
  tmp$chr          <- as.numeric(strsplit(x, "\\.")[[1]][4])
  ## indicate no coloc
  if(is.na(tmp$outcome[1])){
    tmp$outcome <- "no evidence"
  }
  return(tmp)
})
require(plyr)
res.coloc <- do.call(rbind.fill, res.coloc)

## combine to see results
res.coloc <- merge(coloc.regions, res.coloc, by=c("pheno", "chr", "region_start", "region_end"), all = T, suffixes = c(".soma", ".covid"))

## drop colocs with no preserved signal
res.coloc           <- subset(res.coloc, ld.sens >= .8)

## create some columns to ease selection of targets
res.coloc$P4.shared <- ifelse(res.coloc$PP.H4.abf > .7, 1, 0)
res.coloc$LD.same   <- ifelse(res.coloc$ld.check > .8, 1, 0)
res.coloc$candidate <- ifelse(res.coloc$P4.shared == 1 | res.coloc$LD.same == 1, 1, 0)

## extract the interesting bits
ii  <- unique(subset(res.coloc, candidate == 1)$id) ## 35
## extract the interesting bits
tmp <- subset(res.coloc, id %in% ii)
## write to file
write.table(tmp, "Candidate.proteins.coloc.COVID19.20210802.txt", sep="\t", row.names=F)

##############################################
####          load Olink targets          ####
##############################################

olink.targets              <- read.delim("../../../COVID19_v2/02_naive_cis_coloc/input/Olink.mapping.genomic.position.txt", sep="\t", header=T)
table(olink.targets$chromosome_name)
olink.targets$chr          <- ifelse(olink.targets$chromosome_name == "X", 23, as.numeric(olink.targets$chromosome_name))
## omit targets w/o genomic position
olink.targets              <- subset(olink.targets, !is.na(chr))
## add ID
olink.targets$pheno        <- paste0("invn_res_X", tolower(olink.targets$varname)) 
## define regions
olink.targets$region_start <- olink.targets$start_position - 5e5
olink.targets$region_end   <- olink.targets$end_position   + 5e5
## edit start
olink.targets$region_start <- ifelse(olink.targets$region_start < 0, 1, olink.targets$region_start)

## create input file
write.table(olink.targets[, c("pheno", "chr", "region_start", "region_end")], "olink.coloc.regions.txt", sep="\t", row.names=F, col.names=F, quote=F)

#---------------------------------#
##--        load results       --##
#---------------------------------#
 
ii        <- dir("../output_olink/")

## collate results
res.olink <- lapply(ii, function(x){
  ## load results
  tmp <- read.table(paste0("../output_olink/", x), sep="\t", header=T)
  ## add some identifier
  tmp$region_start <- as.numeric(strsplit(x, "\\.")[[1]][5]) 
  tmp$region_end   <- as.numeric(strsplit(x, "\\.")[[1]][6])
  tmp$chr          <- as.numeric(strsplit(x, "\\.")[[1]][4])
  ## indicate no coloc
  if(is.na(tmp$outcome[1])){
    tmp$outcome <- "no evidence"
  }
  return(tmp)
})  
require(plyr)
res.olink <- do.call(rbind.fill, res.olink)

## combine to see results
res.olink <- merge(olink.targets[, c("pheno", "chr", "region_start", "region_end")], res.olink, by=c("pheno", "chr", "region_start", "region_end"), 
                   all = T, suffixes = c(".olink", ".covid"))

## drop colocs with no preserved signal
res.olink           <- subset(res.olink, ld.sens >= .8)

## create some columns to ease selection of targets
res.olink$P4.shared <- ifelse(res.olink$PP.H4.abf > .7, 1, 0)
res.olink$LD.same   <- ifelse(res.olink$ld.check > .8, 1, 0)
res.olink$candidate <- ifelse(res.olink$P4.shared == 1 | res.olink$LD.same == 1, 1, 0)

## extract the interesting bits
ii  <- unique(subset(res.olink, candidate == 1)$pheno)
## extract the interesting bits
tmp <- subset(res.olink, pheno %in% ii)
## write to file
write.table(tmp, "Candidate.Olink.proteins.coloc.COVID19.20210803.txt", sep="\t", row.names=F)

#############################################
#### clarify SomaScan and Olink proteins ####
#############################################

#-------------------#
##--   SomaScan  --##
#-------------------#

## create list of SL proteins with successful COVID-19 results
res.soma.overall <- lapply(unique(res.coloc$id), function(x){
  
  ## get all results
  tmp <- subset(res.coloc, id == x)
  
  ## get all variants
  vf  <- unique(tmp$rsid.covid)
  
  ## return information needed
  return(data.frame(tmp[1, c("pheno", "Target", "SeqId", "TargetFullName", "chr", "region_start", "region_end")], 
                    cand.rsid = paste(vf, collapse = ","), coloc.num = nrow(tmp)))
  
})
res.soma.overall <- do.call(rbind, res.soma.overall)
## N = 1126 genomic regions/aptamers

## write to file for query
write.table(unique(res.soma.overall[, c("pheno", "chr", "region_start", "region_end", "cand.rsid")]), "Variant.query.SomaScan.COVID19.20220321.txt", 
            sep="\t", row.names=F, col.names = F, quote = F)

## --> import results <-- ##

## import output
ii           <- dir("../look_ups/")
## restrict to SomaScan proteins
ii           <- grep("res_invn_X", ii, value=T)
## import
res.soma.tmp <- lapply(ii, function(x){
  
  ## import results
  tmp              <- read.table(paste0("../look_ups/", x), sep="\t", header=T)
  
  ## restrict to the variant with most information
  tmp$ind          <- apply(tmp[, grep("Effect", names(tmp)), drop=F], 1, function(k) sum(is.na(k)))
  tmp$ind.p        <- apply(tmp[, grep("Pvalue\\.", names(tmp)), drop=F], 1, min, na.rm=T)

  ## restrict to the most significant with (almost) complete observations
  tmp              <- tmp[which(tmp$ind == min(tmp$ind)), , drop=F]
  tmp              <- tmp[which(tmp$ind.p == min(tmp$ind.p)), , drop=F]
  
  ## record meta data data
  x                <- strsplit(x, "\\.")[[1]]
  tmp$pheno        <- x[2]
  tmp$region_start <- as.numeric(x[4])
  tmp$region_end   <- as.numeric(x[5])
  
  ## return
  return(tmp)
})
## combine into one data frame
res.soma.tmp <- do.call(rbind, res.soma.tmp)

## align effect estimates
for(j in c("A2", "B1", "B2", "C2")){
  res.soma.tmp[, paste0("Effect.", j)] <- ifelse(res.soma.tmp$Allele1.soma == res.soma.tmp$Allele1.covid, res.soma.tmp[, paste0("Effect.", j)], -res.soma.tmp[, paste0("Effect.", j)])
}

## delete what is no longer needed
res.soma.tmp$Allele1.covid <- res.soma.tmp$Allele2.covid <- NULL
res.soma.tmp[, paste0("Freq1.", c("A2", "B1", "B2", "C2"))] <- NULL 

## add more detailed information
res.soma.overall <- merge(res.soma.overall, res.soma.tmp)

## add coloc results
tmp              <- reshape(res.coloc[, c("pheno", "Ensembl.Gene.ID", "Gene.Name.Name", "chr", "region_start", "region_end", "PP.H3.abf", "PP.H4.abf", "ld.check", "outcome")], 
                            idvar=c("pheno", "Ensembl.Gene.ID", "Gene.Name.Name", "chr", "region_start", "region_end"),
                            direction="wide", timevar = "outcome")

## add to variant output
res.soma.overall <- merge(res.soma.overall, tmp)

#-------------------#
##--    Olink    --##
#-------------------#

## create similar ID for olink
res.olink$id      <- paste(res.olink$pheno, res.olink$chr, res.olink$region_start, res.olink$region_end, sep=".") 

## create list of SL proteins with successful COVID-19 results
res.olink.overall <- lapply(unique(res.olink$id), function(x){
  
  ## get all results
  tmp <- subset(res.olink, id == x)
  
  ## get the most frequently selected variant
  vf  <- unique(tmp$rsid)
  
  ## return information needed
  return(data.frame(tmp[1, c("pheno", "chr", "region_start", "region_end")], 
                    cand.rsid = paste(vf, collapse = ","), coloc.num = nrow(tmp)))
  
})
res.olink.overall <- do.call(rbind, res.olink.overall)
## N = 372 genomic regions/olink assays

## write to file for query
write.table(unique(res.olink.overall[, c("pheno", "chr", "region_start", "region_end", "cand.rsid")]), "Variant.query.Olink.COVID19.20220321.txt", 
            sep="\t", row.names=F, col.names = F, quote = F)

## --> import results <-- ##

## import output
ii            <- dir("../look_ups/")
## restrict to SomaScan proteins
ii            <- grep("invn_res", ii, value=T)
## import
res.olink.tmp <- lapply(ii, function(x){
  
  ## import results
  tmp              <- read.table(paste0("../look_ups/", x), sep="\t", header=T)
  
  ## restrict to the variant with most information
  tmp$ind          <- apply(tmp[, grep("Effect", names(tmp)), drop=F], 1, function(k) sum(is.na(k)))
  tmp$ind.p        <- apply(tmp[, grep("Pvalue\\.", names(tmp)), drop=F], 1, min, na.rm=T)
  
  ## restrict to the most significant with (almost) complete observations
  tmp              <- tmp[which(tmp$ind == min(tmp$ind)), , drop=F]
  tmp              <- tmp[which(tmp$ind.p == min(tmp$ind.p)), , drop=F]
  
  ## record meta data data
  x                <- strsplit(x, "\\.")[[1]]
  tmp$pheno        <- x[2]
  tmp$region_start <- as.numeric(x[4])
  tmp$region_end   <- as.numeric(x[5])
  
  ## return
  return(tmp)
})
## combine into one data frame
res.olink.tmp               <- do.call(rbind, res.olink.tmp)

## edit alleles
res.olink.tmp$Allele1.olink <- tolower(gsub("TRUE", "T", res.olink.tmp$Allele1.olink))
res.olink.tmp$Allele2.olink <- tolower(gsub("TRUE", "T", res.olink.tmp$Allele2.olink))

## align effect estimates
for(j in c("A2", "B1", "B2", "C2")){
  res.olink.tmp[, paste0("Effect.", j)] <- ifelse(res.olink.tmp$Allele1.olink == res.olink.tmp$Allele1.covid, res.olink.tmp[, paste0("Effect.", j)], -res.olink.tmp[, paste0("Effect.", j)])
}

## delete what is no longer needed
res.olink.tmp$Allele1.covid <- res.olink.tmp$Allele2.covid <- NULL
res.olink.tmp[, paste0("Freq1.", c("A2", "B1", "B2", "C2"))] <- NULL 

## add more detailed information
res.olink.overall <- merge(res.olink.overall, res.olink.tmp)

## add coloc results
tmp               <- reshape(res.olink[, c("pheno", "chr", "region_start", "region_end", "PP.H3.abf", "PP.H4.abf", "ld.check", "outcome")], idvar=c("pheno", "chr", "region_start", "region_end"),
                             direction="wide", timevar = "outcome")

## add to variant output
res.olink.overall <- merge(res.olink.overall, tmp)

#############################################
#### combine SomaScan and Olink proteins ####
#############################################

## add information on SomaScan
res.soma.overall$platform_id      <- res.soma.overall$SeqId
res.soma.overall$platform         <- "SomaScan_v4"
names(res.soma.overall)           <- gsub("\\.soma", "", names(res.soma.overall))

## add information on Olink
res.olink.overall$ind             <- NULL
res.olink.overall                 <- merge(res.olink.overall, olink.targets)
res.olink.overall$platform_id     <- res.olink.overall$OlinkID
res.olink.overall$platform        <- "Olink_PEA"
res.olink.overall$Target          <- res.olink.overall$Assay
res.olink.overall$Ensembl.Gene.ID <- res.olink.overall$ensembl_gene_id
res.olink.overall$Gene.Name.Name  <- res.olink.overall$hgnc_symbol
names(res.olink.overall)          <- gsub("\\.olink", "", names(res.olink.overall))

#-------------------------------------#
##--     create combined table     --##
#-------------------------------------#

## intersection of names
ii                            <- intersect(names(res.soma.overall), names(res.olink.overall)) 
res.combined                  <- rbind(res.soma.overall[, ii], res.olink.overall[, ii])

## avoid that Excel maltreats gene names
res.combined$Gene.Name.Name   <- paste0("'", res.combined$Gene.Name.Name)
res.combined$Target           <- paste0("'", res.combined$Target)

## define order 
jj                            <- c("platform_id", "platform", "Target", "Gene.Name.Name", "Ensembl.Gene.ID", "chr", "region_start", "region_end", "rsid", "pos", "Freq1", "Allele1", "Allele2",
                                   "Effect", "StdErr", "Pvalue", "coloc.num", 
                                   "PP.H3.abf.A2", "PP.H4.abf.A2", "ld.check.A2", "n.case.A2", "n.controls.A2", "Effect.A2", "StdErr.A2", "Pvalue.A2", "het_p.A2",
                                   "PP.H3.abf.B1", "PP.H4.abf.B1", "ld.check.B1", "n.case.B1", "n.controls.B1", "Effect.B1", "StdErr.B1", "Pvalue.B1", "het_p.B1",
                                   "PP.H3.abf.B2", "PP.H4.abf.B2", "ld.check.B2", "n.case.B2", "n.controls.B2", "Effect.B2", "StdErr.B2", "Pvalue.B2", "het_p.B2",
                                   "PP.H3.abf.C2", "PP.H4.abf.C2", "ld.check.C2", "n.case.C2", "n.controls.C2", "Effect.C2", "StdErr.C2", "Pvalue.C2", "het_p.C2")

## write to file
write.table(res.combined[,jj], "Results.cis.coloc.SomaLogic.Olink.all.20220322.txt", sep="\t", row.names = F)
