library(optparse)
library(adegenet)
library(ggplot2)
library(plot3D)
library(tidyr)
library(dplyr)

parser <- OptionParser(description = "PCA using adegenet")
parser <- add_option(parser, c("-i", "--input"), type = "character",
                     help = "Input file path")
parser <- add_option(parser, c("-o", "--output_prefix"), type = "character",
                     help = "Output file prefix")
args <- parse_args(parser)

#args <- parse_args(parser, args = c("-i", "try_MHs_SNPs_genotypes_for_adegenet_PCA.txt", "-o", "test_pca_out"))

input_file <- args$input
output_prefix <- args$output_prefix

# PCA using adegenet
dataRaw <- read.table(input_file, header = TRUE, sep = " ", check.names = FALSE)
popDf <- dataRaw[, c(1, 2)]
genotypeDf <- dataRaw[, -c(1, 2)]
rownames(genotypeDf) <- dataRaw[, 1]
adeObjectRaw <- new("genind", tab = genotypeDf)
pop(adeObjectRaw) <- popDf$population
pcaResultRaw <- dudi.pca(adeObjectRaw, nf = 3, scannf = FALSE)
PC_contribution <- pcaResultRaw$eig / sum(pcaResultRaw$eig)

# Draw 3D PCA plot
## Map colors to populations
palette <- c("#AC291F", "#0B5DA6", "#D8731E", "#1D753D", "#6460A1",
             "#0D9488", "#64748B", "#BE185D", "#0E7490", "#7C3AED")
pop_vec <- pop(adeObjectRaw)
pop_vec <- factor(pop_vec, levels = sort(levels(pop_vec)))
n_levels <- nlevels(pop_vec)
if (n_levels > 10) {
  stop("Too many populations: the palette supports at most 10 populations.")
}
pop_colors <- palette[as.integer(pop_vec)]
## Draw 3D plot and save as TIFF (600 DPI)
tiff_filename <- paste0(output_prefix, "_PCA_3dplot.tiff")
tiff(tiff_filename, width = 8, height = 7, units = "in", res = 600)
plot3D::scatter3D(
  x = pcaResultRaw$li[,1],
  y = pcaResultRaw$li[,2],
  z = pcaResultRaw$li[,3],
  pch = 21,      
  cex = 1.5,     
  col = NA,        
  bg = pop_colors,
  xlab = paste0("PC1 (", formatC(signif(PC_contribution[1] * 100, 3), format = "f", digits = 2), "%)"),  
  ylab = paste0("PC2 (", formatC(signif(PC_contribution[2] * 100, 3), format = "f", digits = 2), "%)"),
  zlab = paste0("PC3 (", formatC(signif(PC_contribution[3] * 100, 3), format = "f", digits = 2), "%)"),
  ticktype = "detailed",
  bty = "b2", 
  box = TRUE, 
  theta = -40,   
  phi = 20,      
  d = 10,
  colkey = FALSE,     
  cex.lab = 0.9  
)
legend("right", legend = levels(pop_vec), pch = 16, col = palette[seq_len(n_levels)], bty = "n", cex = 0.9, pt.cex = 1.5, xpd = TRUE)
dev.off()
