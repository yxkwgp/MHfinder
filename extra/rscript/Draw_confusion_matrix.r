Sys.setenv(R_FUTURE_ARGUMENTS = "ignore")

library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-p", "--predFile"), type = "character", help = "predFile")
parser <- add_option(parser, c("-t", "--trueFile"), type = "character", help = "trueFile")
parser <- add_option(parser, c("-o", "--outputFile"), type = "character", help = "outputFile")
args <- parse_args(parser)

library(cvms)
library(caret)
library(ggplot2)
library(ggnewscale)
library(ggimage)

Draw_confusion_matrix <- function(yPred, yTrue) {
    conf_mat <- confusion_matrix(targets = yTrue, predictions = yPred)
    p <- plot_confusion_matrix(conf_mat$`Confusion Matrix`[[1]],
                      font_counts = font(size =4.5,
                                         angle = 0,
                                         color = "Black"),
                      class_order = rev(levels(yTrue)),
                      add_zero_shading = FALSE,
                      add_arrows = FALSE,
                      rm_zero_text = FALSE,
                      add_normalized = FALSE,
                      add_col_percentages = FALSE,
                      add_row_percentages = FALSE,
                      add_sums = TRUE,
                      palette = "Greys",
                      sums_settings = sum_tile_settings(palette = "Purples"),
                      darkness = 0.65) + 
        theme(
          plot.title = element_text(size = 12, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(size = 8, margin = margin(t = 0, r = 2, b = 0, l = 0)),
          axis.title.x = element_text(size = 10, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 10, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.ticks = element_blank(),
          panel.spacing = unit(0, "cm"))
    return(p)
}

yTrue <- readLines(args$trueFile)
yTrue <- factor(yTrue)
yPred <- readLines(args$predFile)
yPred <- factor(yPred)

p <- Draw_confusion_matrix(yPred, yTrue)
plot_height = 1 + (length(levels(yTrue)) + 1) * 0.5
plot_width = 1 + (length(levels(yTrue)) + 1) * 0.5
ggsave(args$outputFile, p, width = plot_width, height = plot_height, dpi = 600)
