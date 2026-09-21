
#ADMIXTURE VISULIZATION
#Load library
library(tidyverse)
library(patchwork)
library(stringr)
library(dplyr)
library(ggplot2)
library(grid)

#Working directory
setwd(PATH)


#Read mapping file 
mapping_lines <- readLines("mapping_file.txt")
mapping_lines <- trimws(mapping_lines)

mapping <- data.frame(
    full_name = sapply(strsplit(mapping_lines, ":"), `[`, 1),
    new_name  = sapply(strsplit(mapping_lines, ":"), `[`, 2),
    sample    = sapply(strsplit(mapping_lines, ":"), `[`, 3),
    stringsAsFactors = FALSE
) %>%
mutate(species = case_when(
    grepl("kisambo",       full_name, ignore.case = TRUE) ~ "E. kisambo",
    grepl("hildebrandtii", full_name, ignore.case = TRUE) ~ "E. hildebrandtii",
    grepl("bubalinus",     full_name, ignore.case = TRUE) ~ "E. bubalinus",
    grepl("sclavoi",       full_name, ignore.case = TRUE) ~ "E. sclavoi",
    TRUE ~ "Unknowns"
))

cat("─── Mapping file loaded ───\n")
print(head(mapping))
cat("\nSamples per species:\n")
print(table(mapping$species))



k1 to 10 ALL SAMPLES
cv_data_all <- read_table("cv_errors_clean.txt",
                           col_names = c("K", "rep", "CV_error"))

cv_summary_all <- cv_data_all %>%
    group_by(K) %>%
    summarize(mean = mean(CV_error), sd = sd(CV_error)) %>%
    arrange(K)

print(cv_summary_all, n = Inf)

p_cv_all <- ggplot(cv_summary_all, aes(x = K, y = mean)) +
    geom_line(color = "black", linewidth = 0.1) +
    #geom_point(size = 0.1, color = "black") +
    scale_x_continuous(breaks = cv_summary_all$K) +
    labs(
        x = "Number of clusters (K)",
        y = "Mean cross-validation error (CV)"
    ) +
    theme_classic() +
    theme(
        text          = element_text(family = "Helvetica"),
        axis.title    = element_text(size = 2.5, face = "bold"),
        axis.text     = element_text(size = 3, color = "black"),
        plot.title    = element_text(size = 5, face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
       axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1)  
    )

print(p_cv_all)

ggsave("CV_error_plot_allsamples.png", p_cv_all, width = 1.5, height = 1, dpi = 600)
ggsave("CV_error_plot_allsamples.pdf", p_cv_all, width = 1.5, height = 1, dpi = 600)


# K=4- all samples, rep3
q4_pruned <- read.table("gardenonly_pruned_clean_maf01_K4_rep3.Q")
fam_pruned <- read.table("gardenonly_pruned_clean_maf01.fam")
q4_pruned$sample <- fam_pruned$V2
q4_pruned <- q4_pruned %>%
    left_join(mapping %>% select(sample, species, new_name), by = "sample") %>%
    mutate(sample = ifelse(!is.na(new_name), new_name, sample))
q4_pruned_long <- q4_pruned %>%
    pivot_longer(cols = starts_with("V"), names_to = "cluster", values_to = "proportion")

cluster4_colors <- c(
    "V1" = "#E69F00",
    "V2" = "#0072B2",
    "V3" = "#CC79A7",
    "V4" = "#59A14F")


p_k4_pruned <- ggplot(q4_pruned_long, aes(x = sample, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", width = 1, linewidth = 0.1, color = "black") +
    facet_grid(~species, scales = "free_x", space = "free", switch = "x") +
    scale_fill_manual(values = cluster4_colors) +
    theme_void() +
    labs(y = "Ancestry proportion", x = NULL, fill = "Cluster") +
    theme(
        text             = element_text(family = "Helvetica"),
        axis.text.x      = element_text(size = 4, angle = 90, hjust = 1, vjust = 1, color = "black",
                                          margin = margin(t = -15)),
        axis.text.y      = element_text(size = 12, color = "black"),
        plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.y     = element_text(size = 14, angle = 90, margin = margin(r = 10)),
        axis.ticks.length.x = unit(0, "pt"),
        strip.text       = element_text(size = 12, face = "italic"),
        strip.placement  = "outside",
        panel.spacing.x  = unit(0, "pt"),
        legend.position  = "none"
    )
print(p_k4_pruned)
ggsave("admixture_K4_allsamples.png", plot = p_k4_pruned, width = 16, height = 7, dpi = 300, bg = "white")

# K=5 - all samples, rep3
q5_pruned <- read.table("gardenonly_pruned_clean_maf01_K5_rep3.Q")
fam_pruned <- read.table("gardenonly_pruned_clean_maf01.fam")
q5_pruned$sample <- fam_pruned$V2
q5_pruned <- q5_pruned %>%
    left_join(mapping %>% select(sample, species, new_name), by = "sample") %>%
    mutate(sample = ifelse(!is.na(new_name), new_name, sample))
q5_pruned_long <- q5_pruned %>%
    pivot_longer(cols = starts_with("V"), names_to = "cluster", values_to = "proportion")

cluster5_colors <- c(
    "V1" = "#E69F00",
    "V2" = "#59A14F",
    "V3" = "#0072B2",
    "V4" = "red",
    "V5" = "#CC79A7"
)

p_k5_pruned <- ggplot(q5_pruned_long, aes(x = sample, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", width = 1, linewidth = 0.1, color = "black") +
    facet_grid(~species, scales = "free_x", space = "free", switch = "x") +
    scale_fill_manual(values = cluster5_colors) +
    theme_void() +
    labs(y = "Ancestry proportion", x = NULL, fill = "Cluster") +
    theme(
        text             = element_text(family = "Helvetica"),
        axis.text.x      = element_text(size = 4, angle = 90, hjust = 1, vjust = 1, color = "black",
                                          margin = margin(t = -15)),
        axis.text.y      = element_text(size = 12, color = "black"),
        plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.y     = element_text(size = 14, angle = 90, margin = margin(r = 10)),
        axis.ticks.length.x = unit(0, "pt"),
        strip.text       = element_text(size = 12, face = "italic"),
        strip.placement  = "outside",
        panel.spacing.x  = unit(0, "pt"),
        legend.position  = "none"
    )
print(p_k5_pruned)
ggsave("admixture_K5_allsamples.png", plot = p_k5_pruned, width = 16, height = 7, dpi = 300, bg = "white")


#K=6 REP7
q6_pruned <- read.table("gardenonly_pruned_clean_maf01_K6_rep7.Q")
fam_pruned <- read.table("gardenonly_pruned_clean_maf01.fam")
q6_pruned$sample <- fam_pruned$V2
q6_pruned <- q6_pruned %>%
    left_join(mapping %>% select(sample, species, new_name), by = "sample") %>%
    mutate(sample = ifelse(!is.na(new_name), new_name, sample))
q6_pruned_long <- q6_pruned %>%
    pivot_longer(cols = starts_with("V"), names_to = "cluster", values_to = "proportion")

cluster6_colors <- c(
    "V1" = "#CC79A7",
    "V2" = "#59A14F",
    "V3" = "#E69F00",
    "V4" = "#0072B2",
    "V5" = "red",
"V6" = "cyan")

p_k6_pruned <- ggplot(q6_pruned_long, aes(x = sample, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", width = 1, linewidth = 0.1, color = "black") +
    facet_grid(~species, scales = "free_x", space = "free", switch = "x") +
    scale_fill_manual(values = cluster6_colors) +
    theme_void() +
    labs(y = "Ancestry proportion", x = NULL, fill = "Cluster") +
    theme(
        text             = element_text(family = "Helvetica"),
        axis.text.x      = element_text(size = 4, angle = 90, hjust = 1, vjust = 1, color = "black",
                                          margin = margin(t = -15)),
        axis.text.y      = element_text(size = 12, color = "black"),
        plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title.y     = element_text(size = 14, angle = 90, margin = margin(r = 10)),
        axis.ticks.length.x = unit(0, "pt"),
        strip.text       = element_text(size = 12, face = "italic"),
        strip.placement  = "outside",
        panel.spacing.x  = unit(0, "pt"),
        legend.position  = "none"
    )
print(p_k6_pruned)
ggsave("admixture_K6_allsamples.png", plot = p_k6_pruned, width = 16, height = 7, dpi = 300, bg = "white")

