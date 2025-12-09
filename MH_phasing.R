#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

library(haplo.stats)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("Usage: Rscript MH_phasing.R <input_tsv> <output_tsv>")
}

input_file <- args[1]
output_file <- args[2]

geno_data <- read.csv(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
locus_label_all <- colnames(geno_data)
locus_label <- gsub("\\.a1$", "", locus_label_all[seq(1, length(locus_label_all), by = 2)])

save.em <- haplo.em(geno = geno_data, locus.label = locus_label, miss.val = c(0, NA))
hap1codeAll <- sapply(1:nrow(geno_data), function(i) {idx = save.em$subj.id == i; save.em$hap1code[idx][which.max(save.em$post[idx])];})
hap2codeAll <- sapply(1:nrow(geno_data), function(i) {idx = save.em$subj.id == i; save.em$hap2code[idx][which.max(save.em$post[idx])];})
MH_label <- sapply(1:nrow(save.em$haplotype), function(x) paste(save.em$haplotype[x,], collapse = '-'))
Genotype1 <- factor(hap1codeAll, levels=1:length(MH_label),labels=MH_label)
Genotype2 <- factor(hap2codeAll, levels=1:length(MH_label),labels=MH_label)
GenotypeAll <- data.frame(Genotype1, Genotype2)
names(GenotypeAll) <- c("temp.A1", "temp.A2")

write.table(GenotypeAll, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Haplotype analysis completed. Output written to:", output_file, "\n")
