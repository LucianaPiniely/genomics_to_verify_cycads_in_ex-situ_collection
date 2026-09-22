#PCoA
#We had two subset for PCoA so I did for both subset

#subset 1 (All samples)
#Read PLINK files
plink_data <- read.plink("gardenonly_pruned.bed",
                          "gardenonly_pruned.bim",
                          "gardenonly_pruned.fam")

# extract genotype matrix and assign sample names
gt <- as(plink_data$genotypes, "numeric")
rownames(gt) <- plink_data$fam$member
cat("Genotype matrix:", nrow(gt), "samples x", ncol(gt), "SNPs\n")

# Read mapping file (my samples has different names so I had to rename)
tipnames <- readLines("mapping_file.txt")
tipnames <- trimws(tipnames)
mapping <- data.frame(
    full_name  = sapply(strsplit(tipnames, ":"), `[`, 1),
    new_name   = sapply(strsplit(tipnames, ":"), `[`, 2),
    sample     = sapply(strsplit(tipnames, ":"), `[`, 3),
    stringsAsFactors = FALSE
) %>%
mutate(species = case_when(
    grepl("kisambo",       full_name, ignore.case = TRUE) ~ "E. kisambo",
    grepl("hildebrandtii", full_name, ignore.case = TRUE) ~ "E. hildebrandtii",
    grepl("bubalinus",     full_name, ignore.case = TRUE) ~ "E. bubalinus",
    grepl("sclavoi",       full_name, ignore.case = TRUE) ~ "E. sclavoi",
    grepl("kanga",         full_name, ignore.case = TRUE) ~ "E. kanga",
    TRUE ~ "Unknowns"
))

#Impute missing and run PCoA
for (i in 1:ncol(gt)) {
    gt[is.na(gt[,i]), i] <- mean(gt[,i], na.rm = TRUE)
}

# compute distance matrix — samples are rows
dist_mat <- dist(gt, method = "euclidean")
pcoa_result <- pcoa(dist_mat)

coords <- as.data.frame(pcoa_result$vectors[, 1:2])
colnames(coords) <- c("PCoA1", "PCoA2")
coords$sample <- rownames(coords)
var_explained <- pcoa_result$values$Relative_eig * 100

# join with mapping and rename samples
coords <- coords %>%
    left_join(mapping %>% select(sample, species, new_name), by = "sample") %>%
    mutate(sample = ifelse(!is.na(new_name), new_name, sample)) %>%
    select(-new_name)
coords$species[is.na(coords$species)] <- "Unknowns"

species_levels <- c("E. bubalinus", "E. hildebrandtii", "E. kisambo", "E. sclavoi", "Unknowns")

coords$species <- factor(coords$species, levels = species_levels)
coords_sub2$species <- factor(coords_sub2$species, levels = species_levels)

print(head(coords$sample))
print(table(coords$species))

# PERMANOVA
# ─────────────────────────────────────────────────────────────────────────────
group <- factor(coords$species)
adonis_result <- vegan::adonis2(dist_mat ~ group, permutations = 999)
R2 <- round(as.numeric(adonis_result$R2[1]), 3)
pvalue <- signif(adonis_result$`Pr(>F)`[1], 2)
print(adonis_result)
cat("R2:", R2, "\n")
cat("P-value:", pvalue, "\n")

# Plot
samples_to_label <- c("NBG_U6", "HEID_H1", "MBG_S1", "DG_S6", "BOGA_S1", "HNT_H14", "DBG_H1")

p_A <- ggplot(coords, aes(x = PCoA1, y = PCoA2, color = species)) +
    geom_point(size = 1, alpha = 0.8) +
    geom_text_repel(
        data = coords %>% filter(sample %in% samples_to_label),
        aes(label = sample),
        size = 2,
        max.overlaps = Inf,
        show.legend = FALSE
    ) +
    annotate("text",
        x = Inf, y = Inf,
        label = paste0("Adonis R² = ", R2, "\nP-value = ", pvalue),
        hjust = 1.1, vjust = 1.5,
        size = 2,
        color = "black") +
    scale_color_manual(values = c(
        "E. kisambo"       = "#E69F00",
        "E. hildebrandtii" = "#59A14F",
        "E. bubalinus"     = "#0072B2",
        "E. sclavoi"       = "#CC79A7",
        "Unknowns"         = "red"
    )) +
    labs(
        x = paste0("PCoA1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PCoA2 (", round(var_explained[2], 1), "%)"),
         title = " PCoA of subset A samples",
        color = "Species"
    ) +
    coord_cartesian(xlim = c(-6, 6), ylim = c(-4, 15)) +
    theme_bw() +
    theme( text = element_text(family = "Times New Roman"), 
        legend.position = "right",
        legend.title    = element_text(size = 10),
        legend.text     = element_text(size = 10, color = "black", face = "italic"),
        axis.title      = element_text(size = 10, face = "bold"),
       plot.title = element_text(size = 10, color = "black", face = "bold"),
        axis.text       = element_text(size = 8, color = "black"),
        plot.subtitle   = element_text(size = 8, hjust = 0.5, color = "black"),
         panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("PCoA_allsamples.png",
       width = 10,
       height = 10,
       dpi = 600)


#Subset 2 (Those that didn't show clear separation on subset 1)
#Read PLINK files
plink_subset2 <- read.plink("gardenonly_pruned_subset.bed",
                             "gardenonly_pruned_subset.bim",
                             "gardenonly_pruned_subset.fam")

gt_sub2 <- as(plink_subset2$genotypes, "numeric")
rownames(gt_sub2) <- plink_subset2$fam$member
cat("Subset2 genotype matrix:", nrow(gt_sub2), "samples x", ncol(gt_sub2), "SNPs\n")

# Read mapping file
tipnames <- readLines("mapping_file.txt")
tipnames <- trimws(tipnames)
mapping <- data.frame(
    full_name  = sapply(strsplit(tipnames, ":"), `[`, 1),
    new_name   = sapply(strsplit(tipnames, ":"), `[`, 2),
    sample     = sapply(strsplit(tipnames, ":"), `[`, 3),
    stringsAsFactors = FALSE
) %>%
mutate(species = case_when(
    grepl("kisambo",       full_name, ignore.case = TRUE) ~ "E. kisambo",
    grepl("hildebrandtii", full_name, ignore.case = TRUE) ~ "E. hildebrandtii",
    grepl("bubalinus",     full_name, ignore.case = TRUE) ~ "E. bubalinus",
    grepl("sclavoi",       full_name, ignore.case = TRUE) ~ "E. sclavoi",
    grepl("kanga",         full_name, ignore.case = TRUE) ~ "E. kanga",
    TRUE ~ "Unknowns"
))

#Impute missing and run PCoA
for (i in 1:ncol(gt_sub2)) {
    gt_sub2[is.na(gt_sub2[,i]), i] <- mean(gt_sub2[,i], na.rm = TRUE)
}

dist_mat_sub2 <- dist(gt_sub2, method = "euclidean")
pcoa_sub2 <- pcoa(dist_mat_sub2)

coords_sub2 <- as.data.frame(pcoa_sub2$vectors[, 1:2])
colnames(coords_sub2) <- c("PCoA1", "PCoA2")
coords_sub2$sample <- rownames(coords_sub2)
var_explained_sub2 <- pcoa_sub2$values$Relative_eig * 100

# join with mapping
coords_sub2 <- coords_sub2 %>%
    left_join(mapping %>% select(sample, species, new_name), by = "sample") %>%
    mutate(sample = ifelse(!is.na(new_name), new_name, sample)) %>%
    select(-new_name)
coords_sub2$species[is.na(coords_sub2$species)] <- "Unknowns"
species_levels <- c("E. bubalinus", "E. hildebrandtii", "E. kisambo", "E. sclavoi", "Unknowns")

coords$species <- factor(coords$species, levels = species_levels)
coords_sub2$species <- factor(coords_sub2$species, levels = species_levels)

print(table(coords_sub2$species))

#PERMANOVA
group_sub2 <- factor(coords_sub2$species)
adonis_sub2 <- vegan::adonis2(dist_mat_sub2 ~ group_sub2, permutations = 999)
R2_sub2 <- round(as.numeric(adonis_sub2$R2[1]), 3)
pvalue_sub2 <- signif(adonis_sub2$`Pr(>F)`[1], 2)
cat("R2:", R2_sub2, "\n")
cat("P-value:", pvalue_sub2, "\n")

#Plot
samples_to_label_sub2 <- c("NBG_U6", "HEID_H1", "DG_S6", "LPC_H2", "NBG_H4")

p_B <- ggplot(coords_sub2, aes(x = PCoA1, y = PCoA2, color = species)) +
    geom_point(size = 1, alpha = 0.8) +
    geom_text_repel(
        data = coords_sub2 %>% filter(sample %in% samples_to_label_sub2),
        aes(label = sample),
        size = 2,
        max.overlaps = Inf,
        show.legend = FALSE
    ) +
    annotate("text",
        x = Inf, y = Inf,
        label = paste0("Adonis R² = ", R2_sub2, "\nP-value = ", pvalue_sub2),
        hjust = 1.1, vjust = 1.5,
        size = 2,
        color = "black") +
    scale_color_manual(values = c(
        "E. kisambo"       = "#E69F00",
        "E. hildebrandtii" = "#59A14F",
        "E. bubalinus"     = "#0072B2",
        "E. sclavoi"       = "#CC79A7",
        "Unknowns"         = "red"
    )) +
    labs(
        x = paste0("PCoA1 (", round(var_explained_sub2[1], 1), "%)"),
        y = paste0("PCoA2 (", round(var_explained_sub2[2], 1), "%)"),
        color = "Species",
        title = "PCoA of subset B samples "
    ) +
    theme_bw() +
    coord_cartesian(xlim = c(-6, 10), ylim = c(-6, 6)) +
    theme(
        text             = element_text(family = "Helvetica"),
        legend.position  = "right",
        legend.title     = element_text(size = 10),
        legend.text      = element_text(size = 10, color = "black", face = "italic"),
        plot.title = element_text(size = 10, color = "black", face = "bold"),
        axis.title       = element_text(size = 10, face = "bold"),
       rplot.title = element_text(size = 10, color = "black", face = "bold", hjust = 0.5),
plot.tag.position = "top",
        axis.text        = element_text(size = 8, color = "black"),
        plot.subtitle    = element_text(size = 10, hjust = 0.5, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )

ggsave("PCoA_subset2.png",
       width = 13,
       height = 10,
       dpi = 600)

       #END