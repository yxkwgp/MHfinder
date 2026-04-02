library(optparse)
library(purrr)
library(readr)
library(readxl)
library(writexl)
library(dplyr)
library(ggplot2)
library(patchwork)

parser <- OptionParser(description = "AUC and ACC decline with number of MHs")
parser <- add_option(parser, c("-d", "--input_dir"), type = "character",
                     help = "Result directory of random kept MHs")
parser <- add_option(parser, c("-i", "--max_n_mh_input_prefix"), type = "character",
                     help = "Result file prefix of LOOCV of all MHs")
parser <- add_option(parser, c("-o", "--output_prefix"), type = "character",
                     help = "Output file prefix")
args <- parse_args(parser)
#args <- parse_args(parser, args = c("-d", "../08.Naive_Bayes_random_kept_MHs", "-i", "Naive_Bayes_48MHs", "-o", "Random_kept_48MHs"))

input_dir <- args$input_dir
input_dir <- sub('/$', '', input_dir)
input_dir <- paste0(input_dir, '/')
max_n_mh_input_prefix <- args$max_n_mh_input_prefix
output_prefix <- args$output_prefix

# Get populations from max N MHs result
max_n_mh_yTrue_file <- paste0(max_n_mh_input_prefix, '_yTrue.txt')
populations <- readLines(max_n_mh_yTrue_file) %>% as.vector() %>% unique() %>% sort()
col_names <- c("mean", populations)

# Get AUC and ACC of max N MHs
evaluation_matrix_top_df <- read_excel(paste0(max_n_mh_input_prefix, '_evaluation_metrics.xlsx')) %>% as.data.frame()
rownames(evaluation_matrix_top_df) <- evaluation_matrix_top_df[, 1]
evaluation_matrix_top_df <- evaluation_matrix_top_df[, -1]
AUC_top <- as.numeric(evaluation_matrix_top_df$AUC)
names(AUC_top) <- populations
ACC_top <- readLines(paste0(max_n_mh_input_prefix, '_accuracy.txt')) %>% as.numeric()

# Get AUC and ACC of random kept MHs
all_result_dirs <- list.dirs(input_dir, full.names = FALSE, recursive = FALSE)
n_upper_bound <- max(as.integer(vapply(strsplit(all_result_dirs, "_"), function(x) x[2], character(1))), na.rm = TRUE)
AUC_random_df <- map_dfr(1:n_upper_bound, function(i) {
    file_path <- file.path(input_dir, paste0("kept_", i), "AUC.txt")
    read_tsv(file_path, col_names = col_names, show_col_types = FALSE) %>%
        mutate(n_MHs = i)
})
ACC_random_df <- map_dfr(1:n_upper_bound, function(i) {
    file_path <- file.path(input_dir, paste0("kept_", i), "ACC.txt")
    acc_values <- as.numeric(readLines(file_path))
    data.frame(ACC = acc_values, n_MHs = i)
})

# Map color to populations
pop_colors <- c("#AC291F", "#0B5DA6", "#D8731E", "#1D753D", "#6460A1",
             "#0D9488", "#64748B", "#BE185D", "#0E7490", "#7C3AED")
if(length(names(populations)) > 10) {
  stop("Too many populations: the palette supports at most 10 populations.")
}
color_map <- setNames(pop_colors[1:length(populations)], populations)

# Function of Draw boxplot
plot_breaks <- c(n_upper_bound + 1, seq(floor((n_upper_bound + 1) / 10) * 10, 10, by = -10), 1)
plot_auc_decline <- function(population, color) { 
    p <- ggplot(AUC_random_df, aes(x = factor(n_MHs), y = !!sym(population))) +
        geom_boxplot(fill = color, outlier.size = 0.5) +
        geom_point(aes(x = as.character(n_upper_bound + 1), y = AUC_top[population]), color = color, size = 1.5) +
        annotate("text", x = as.character(round((n_upper_bound + 1) / 3 * 2)), y = 0.75, label = population, size = 10, fontface = "bold") +
        labs(
            x = "Number of Markers",
            y = "AUC"
        ) +
        coord_transform(ylim = c(0.5, 1)) +
        theme_classic() +
        theme(
            axis.title = element_text(size = 25),
            axis.text = element_text(size = 20)
        ) +
        scale_x_discrete(limits = c(as.character(n_upper_bound + 1), rev(levels(factor(AUC_random_df$n_MHs)))),
                         breaks = plot_breaks)
    
    return(p)
}

# Draw AUC decline plot for each population and store in a list
auc_decline_plot_list <- lapply(populations, function(pop) plot_auc_decline(pop, color_map[pop]))
names(auc_decline_plot_list) <- populations

p_ACC <- ggplot(ACC_random_df, aes(x = factor(n_MHs), y = ACC)) +
    geom_boxplot(fill = "grey60", outlier.size = 0.5) +
    geom_point(aes(x = as.character(n_upper_bound + 1), y = ACC_top), color = "grey60", size = 1.5) +
    labs(
        x = "Number of Markers",
        y = "Accuracy"
    ) +
    coord_transform(ylim = c(0, 1)) +
    theme_classic() +
    theme(
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20)
    ) +
    scale_x_discrete(limits = c(as.character(n_upper_bound + 1), rev(levels(factor(ACC_random_df$n_MHs)))),
                     breaks = plot_breaks)

# Combine plots and save
auc_decline_plot_list <- c(auc_decline_plot_list, list(ACC = p_ACC))
combined_plot <- wrap_plots(auc_decline_plot_list, ncol = 2) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(
    plot.tag = element_text(face = "bold", size = 20),
    plot.tag.position = c(0.05, 1.05),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)
  )
plot_width = 14
if(length(auc_decline_plot_list) %% 2 == 0) {
  plot_height = length(auc_decline_plot_list) / 2 * 6
} else {
  plot_height = (length(auc_decline_plot_list) + 1) / 2 * 6
}
ggsave(paste0(output_prefix, "_AUC_ACC_decline_with_n_MHs.tiff"), combined_plot, width = plot_width, height = plot_height, dpi = 600)

# Calculate the median of each population and its decline by the number of MHs
median_by_n_MHs <- AUC_random_df %>%
    group_by(n_MHs) %>%
    summarise(
        across(all_of(populations), median, .names = "Median AUC ({.col})"), 
        `Mininum Mean AUC` = min(mean),
        .groups = "drop"
    ) %>%
    bind_rows(tibble(n_MHs = n_upper_bound + 1, !!!setNames(AUC_top[populations], paste0("Median AUC (", populations, ")")), mean_minimum = mean(AUC_top))) %>%
    arrange(n_MHs) %>%
    mutate(
        across(starts_with("Median AUC ("), ~ lead(.x) - .x, .names = "{.col} Reduction"),
        `Minimum AUC Reduction` = lead(`Mininum Mean AUC`) - `Mininum Mean AUC`
    ) %>%
    rename_with(~ sub("^(.+?) \\((.+)\\) Reduction$", "\\1 Reduction (\\2)", .x), ends_with(" Reduction")) %>%
    arrange(desc(n_MHs)) %>%
    rename(`N MH` = n_MHs) %>%
    select(`N MH`, `Mininum Mean AUC`, `Minimum AUC Reduction`, !!!map(populations, ~ c(paste0("Median AUC (", .x, ")"), paste0("Median AUC Reduction (", .x, ")"))) %>% flatten_chr())

write_xlsx(median_by_n_MHs, path = paste0(output_prefix, "_AUC_median_by_n_MHs.xlsx"))
