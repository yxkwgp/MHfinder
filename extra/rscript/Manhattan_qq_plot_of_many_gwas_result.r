library(optparse)
library(data.table)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)

parser <- OptionParser(description = "Manhattan and QQ plot for many GWAS results")
parser <- add_option(parser, c("-p", "--p_value_threshold"), type = "numeric",
                     help = "P-value threshold (float)")
parser <- add_option(parser, c("-g", "--gwas_result_files"), type = "character",
                     help = "One or more GWAS result file paths, comma-separated (e.g. file1.txt,file2.txt)")
parser <- add_option(parser, c("-o", "--output_prefix"), type = "character",
                     help = "Output file prefix")
args <- parse_args(parser)

#args <- parse_args(parser, args = c("-p", "5e-8", "-g", "SCGWAS_AFR_merge.assoc.logistic, SCGWAS_AMR_merge.assoc.logistic, SCGWAS_EAS_merge.assoc.logistic, SCGWAS_EUR_merge.assoc.logistic, SCGWAS_SAS_merge.assoc.logistic", "-o", "try_manhattan_qq"))

p_value_threshold <- args$p_value_threshold
gwas_result_files <- trimws(strsplit(args$gwas_result_files, ",")[[1]])
output_prefix <- args$output_prefix

# Read each GWAS result file and preprocessing
gwas_result_list <- lapply(gwas_result_files, data.table::fread)
names(gwas_result_list) <- vapply(gwas_result_files, function(filename) strsplit(basename(filename), "_")[[1]][2], character(1L))
gwas_result_list <- lapply(gwas_result_list, function(df) dplyr::mutate(df, log10P_minus = -log10(P)))
min_P = min(sapply(gwas_result_list, function(df) min(df$P, na.rm =TRUE)))
CHR_len_df = bind_rows(gwas_result_list) %>%
    group_by(CHR) %>%
    summarise(CHR_len = max(BP)-min(BP))
CHR_POS_df <- CHR_len_df  %>% 
  mutate(total = cumsum(as.numeric(CHR_len)) - as.numeric(CHR_len)) %>%
  select(-CHR_len)
n_CHR = length(unique(CHR_len_df$CHR))
X_axis_CHR_df <- bind_rows(gwas_result_list) %>%
    right_join(CHR_POS_df, by="CHR") %>%
    mutate(BPcum = BP + total) %>%
    group_by(CHR) %>%
    summarize(center=(max(BPcum) + min(BPcum)) / 2 ) %>%
    filter(CHR %in% c(1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 18, 22))

# Map color to populations
pop_colors <- c("#AC291F", "#0B5DA6", "#D8731E", "#1D753D", "#6460A1",
             "#0D9488", "#64748B", "#BE185D", "#0E7490", "#7C3AED")
if(length(names(gwas_result_list)) > 10) {
  stop("Too many populations: the palette supports at most 10 populations.")
}
color_map <- setNames(pop_colors[1:length(names(gwas_result_list))], names(gwas_result_list))

# Function: Draw Manhattan plot
draw_manhattan_plot <- function(gwas_data, p_value_threshold, chr_pos_df, x_axis_df, n_chr, chr_color) {
    ## Process data
    SNP_pos_df <- gwas_data %>%
        right_join(chr_pos_df, by="CHR") %>%
        mutate(BPcum = BP + total) %>%
        filter(P < 0.01) %>%
        mutate(significant = ifelse(P < p_value_threshold, "yes", "no"))
    
    p <- ggplot(SNP_pos_df, aes(x=BPcum, y=log10P_minus)) +
        geom_point(aes(color=as.factor(CHR)), size=1) +
        scale_color_manual(values = rep(chr_color, n_chr)) +
        geom_point(data=subset(SNP_pos_df, significant=="yes"), color="red", size=1.5) +
        geom_text_repel(data=subset(SNP_pos_df, significant=="yes"), aes(label=SNP), size=4.5, ylim = c(-log10(p_value_threshold), NA), segment.color = "red") +
        scale_x_continuous(label = x_axis_df$CHR, breaks = x_axis_df$center, expand = c(0.03, 0.03)) +
        scale_y_continuous(breaks = c(2, 5, 10, 20, 50, 100, 200), expand = c(0.03, 0)) +
        coord_transform(y = 'log10', ylim = c(2, 200)) +
        geom_hline(yintercept = -log10(p_value_threshold), color = 'red', linewidth = 0.5) + 
        labs(x = "Chromosome", y = expression(-log[10](italic(p)))) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text = element_text(size = 18),
            axis.text.y = element_text(size = 18, angle = 90, hjust = 0.5),
            axis.title = element_text(size = 20),
            plot.margin = margin(t = 20, r = 5, b = 5, l = 5)
        )
    
    return(p)
}

# Draw Manhattan plot for each GWAS result
manhattan_plot_list <- lapply(names(gwas_result_list), function(pop_name) {
  chr_color <- unname(c(color_map[pop_name], "grey60"))
  draw_manhattan_plot(gwas_result_list[[pop_name]], p_value_threshold, CHR_POS_df, X_axis_CHR_df, n_CHR, chr_color)
})
names(manhattan_plot_list) <- names(gwas_result_list)

# Draw QQ plot
gwas_result_merged <- bind_rows(gwas_result_list, .id = "population") %>%
  filter(!is.na(P), P > 0) %>%
  group_by(population) %>%
  arrange(P) %>%
  rename(observed = log10P_minus) %>%
  mutate(expected = -log10(ppoints(n()))) %>%
  ungroup()

plot_qq <- ggplot(gwas_result_merged, aes(x = expected, y = observed, color = population)) +
    geom_point(size = 1.5) +
    geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 0.8) +
    scale_color_manual(values = color_map) +
    coord_transform(ylim = c(0, 200)) +
    labs(
        x = expression(Expected ~ -log[10](italic(p))),
        y = expression(Observed ~ -log[10](italic(p))),
        color = "Population"
    ) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle = 90, hjust = 0.5),
        axis.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.position = c(0.1, 0.9),
        legend.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5)
    )

# Patchwork and save
manhattan_plot_list <- c(manhattan_plot_list, list(QQ = plot_qq))
combined_plot <- wrap_plots(manhattan_plot_list, ncol = 2) +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(
    plot.tag = element_text(face = "bold", size = 20),
    plot.tag.position = c(0.05, 1.05),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)
  )
plot_width = 12
if(length(manhattan_plot_list) %% 2 == 0) {
  plot_height = length(manhattan_plot_list) / 2 * 6
} else {
  plot_height = (length(manhattan_plot_list) + 1) / 2 * 6
}
ggsave(paste0(output_prefix, "_manhattan_QQ_plots.tiff"), combined_plot, width = plot_width, height = plot_height, dpi = 600)
