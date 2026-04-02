#!/usr/bin/env Rscript
library(haplo.stats)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-p", "--prefix"), type = "character", default = NULL,
              help = "PLINK text file prefix (.map and .ped)"),
  make_option(c("-m", "--mhs"), type = "character", default = NULL,
              help = "Selected MHs file path"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output genotype file path")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list, 
                          description = "MHs phasing script - Infer MH genotypes from SNP genotype data",
                          usage = "Rscript %prog -p <PLINK_prefix> -m <MHs_file> -o <output_file>")
opt <- parse_args(opt_parser)

# Check required parameters
if (is.null(opt$prefix) || is.null(opt$mhs) || is.null(opt$output)) {
  cat("Error: Missing required parameters\n")
  print_help(opt_parser)
  quit(status = 1)
}

# Extract parameters and construct file paths
prefix <- opt$prefix
map_file <- paste0(prefix, ".map")
ped_file <- paste0(prefix, ".ped")
mhs_file <- opt$mhs
output_file <- opt$output

# Function to get genotyped data.frame, GTbyHS = GenotypeByHaplo.stats
GTbyHS <- function(MH_colname, label, geno){
  save.em <- haplo.em(geno=geno, locus.label= label, miss.val=c(0,NA,'0'))
  hap1codeAll <- sapply(1:nrow(geno), function(i) {idx = save.em$subj.id == i; save.em$hap1code[idx][which.max(save.em$post[idx])];})
  hap2codeAll <- sapply(1:nrow(geno), function(i) {idx = save.em$subj.id == i; save.em$hap2code[idx][which.max(save.em$post[idx])];})
  MH_label <- sapply(1:nrow(save.em$haplotype), function(x) paste(save.em$haplotype[x,], collapse = '-'))
  Genotype1 <- factor(hap1codeAll, levels=1:length(MH_label),labels=MH_label)
  Genotype2 <- factor(hap2codeAll, levels=1:length(MH_label),labels=MH_label)
  GenotypeAll <- data.frame(Genotype1, Genotype2)
  names(GenotypeAll) <- MH_colname
  rownames(GenotypeAll) <- rownames(geno)
  return(GenotypeAll)
}

# Load SNP genotype data
dfMap <- read.table(map_file, header=FALSE, col.names=c('chr', 'rsID', 'cM', 'bp'))

rsIDdup <- c(rbind(paste0(dfMap$rsID, ".a1"), paste0(dfMap$rsID, ".a2")))
pedHeader <- c(c('FID', 'IID', 'father', 'mother', 'sex', 'phe'), rsIDdup)

dfPed <- read.table(ped_file, header=FALSE, colClasses = "character")
colnames(dfPed) <- pedHeader

# Estimate MH genotypes using maximum likelihood
MHs <- readLines(mhs_file)

MHs_genotype <- dfPed[, 1:2]

for(rownum in 1:length(MHs)){
  if(rownum%%3 == 1){
    MHname <- sub("^>", "", MHs[rownum])
    MH_colname <- c(paste0(MHname,".A1"),paste0(MHname,".A2"))
  }
  if(rownum%%3 == 0){
    label <- unlist(strsplit(MHs[rownum],split="\t"))
    selectcol <- c(rbind(paste0(label, ".a1"), paste0(label, ".a2")))
    geno_sel <- subset(dfPed, select = selectcol)
    geno_tem <- GTbyHS(MH_colname=MH_colname, label=label, geno=geno_sel)
    MHs_genotype <- cbind(MHs_genotype, geno_tem)
  }
}
write.table(MHs_genotype, file=output_file, sep="\t", row.names=FALSE, quote=FALSE)
