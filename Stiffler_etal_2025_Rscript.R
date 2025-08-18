############### Stiffler et al. 2025 Statistcal Analysis and Figure Generation R Script ###############
############### Prophage induction alters microbiome composition and functional capacity in a Sargassum-derived biofilm ###############

rm(list = ls())
library(readxl); library(pheatmap); library(tidyr); library(ggplot2); library(ggExtra); library(dplyr); library(ggpubr);
library(vegan); library(pairwiseAdonis); library(genoPlotR); library(stringr); library(tidyverse); library(IRanges)

############### Figure 1: Mitomycin C induces pelagic Sargassum biofilm-associated prophages ###############
##### Figure 1b: Virus-to-host Ratio Heatmap #####
VHR <- read_excel('/path/to/Supplemental_DataSheet4.xlsx')

VHR$prophage_host_ratio <- as.numeric(VHR$prophage_host_ratio)

VHR[is.na(VHR)] <- 0

VHR_less_1 <- VHR %>% filter(VHR$prophage_host_ratio < 1) 
VHR_great_1 <- VHR %>% filter(VHR$prophage_host_ratio > 1) 

VHR_stats <- VHR %>% group_by(prophage,sample_type) %>% summarize(
  n = n(),
  xbar = mean(prophage_host_ratio),
  median = median(prophage_host_ratio),
  sd = sd(prophage_host_ratio),
  se = sd/sqrt(n))

VHR_wide_stats <- VHR_stats %>% select(sample_type, prophage, median) %>% pivot_wider(names_from = prophage, values_from = median)

VHR_wide_stats_matrix <- as.matrix(VHR_wide_stats[,2:12])

rownames(VHR_wide_stats_matrix) = sapply(VHR_wide_stats$sample_type,function(x)
  strsplit(as.character(x),split = "\\\\")[[1]][1])

pheatmap(mat = VHR_wide_stats_matrix, scale = "column",
         fontsize = 8, angle_col = 45, cellheight = 10, cellwidth = 15,
         cluster_rows=T, cluster_col=T, na_col = 'white',
         color=colorRampPalette(c("#cefae4","#e0dd6c","#ba5619","#2a061b"))(50))

VHR_wide <- VHR %>% select(sample_type, replicate, prophage, prophage_host_ratio) %>% pivot_wider(names_from = prophage, values_from = prophage_host_ratio)

VHR_matrix <- as.matrix(VHR_wide[,3:13])

VHR.dist<-vegdist(VHR_matrix, method='bray')
VHR.dist

set.seed(36)
VHR.dif <-adonis2(VHR.dist~as.factor(VHR_wide$sample_type), data=VHR_wide, permutations=9999)
VHR.dif

pairwise.adonis2(VHR.dist~sample_type, data=VHR_wide, permutations=9999, methods='bray')

VHR.simper <- simper(VHR_matrix, VHR_wide$sample_type, permutations = 9999)
VHR.simper.summary <- summary(VHR.simper, ordered = TRUE, digits = max(3,getOption("digits") - 3))

##### Figure 1c: Prophage genome plots #####
phage_annotations <- read_excel('/path/to/Supplemental_DataSheet5.xlsx')
head(phage_annotations)

phage_annotations  <- phage_annotations  %>% mutate(col = ifelse(group == 'Hypothetical','#DDDDDD',col))
phage_annotations  <- phage_annotations  %>% mutate(fill = ifelse(group == 'Hypothetical','#DDDDDD',fill))

phage_annotations  <- phage_annotations  %>% mutate(col = ifelse(group == 'Replication','#ba5619',col))
phage_annotations  <- phage_annotations  %>% mutate(fill = ifelse(group == 'Replication','#ba5619',fill))

phage_annotations  <- phage_annotations  %>% mutate(col = ifelse(group == 'Temperate life cycle','#664C5B',col))
phage_annotations  <- phage_annotations  %>% mutate(fill = ifelse(group == 'Temperate life cycle','#664C5B',fill))

phage_annotations  <- phage_annotations  %>% mutate(col = ifelse(group == 'Host evasion','#90AFA0',col))
phage_annotations  <- phage_annotations  %>% mutate(fill = ifelse(group == 'Host evasion','#90AFA0',fill))

phage_annotations  <- phage_annotations  %>% mutate(col = ifelse(group == 'Other','#cefae4',col))
phage_annotations  <- phage_annotations  %>% mutate(fill = ifelse(group == 'Other','#cefae4',fill))

phage_annotations  <- phage_annotations  %>% mutate(col = ifelse(group == 'Packaging and Lysis','#e0dd6c',col))
phage_annotations  <- phage_annotations  %>% mutate(fill = ifelse(group == 'Packaging and Lysis','#e0dd6c',fill))

phage_annotations  <- phage_annotations  %>% mutate(col = ifelse(group == 'Structural','#948816',col))
phage_annotations  <- phage_annotations  %>% mutate(fill = ifelse(group == 'Structural','#948816',fill))


####### Fig. 1c Prophage genome plots #######
#seperate into the different genomes 
k127_162855 <- phage_annotations %>% filter(name == "MC-12PostInd-BF1_megahit_k127_555745_flag_0_multi_68.4010_len_7549_extended_partial")
k127_567814 <- phage_annotations %>% filter(name == "CT-12PostInd-BF4_megahit_k127_431539_flag_0_multi_14.0000_len_69873_extended_partial")
k127_233705 <- phage_annotations %>% filter(name == "MC-12PostInd-BF1_megahit_k127_558051_flag_0_multi_34.4574_len_15183")
k127_173951 <- phage_annotations %>% filter(name == "CT-12PostInd-BF1_megahit_k127_173951_flag_0_multi_10693.0828_len_47076|provirus_16859_46428")

# pull out the coumns necessary for genoPlotR 
# name, start, end, strand, col, fill, lty, lwd, pch, cex, gene_type
# and convert them to a data frame
k127_162855_df2 <- data.frame(k127_162855[,2:12])
k127_567814_df2 <- data.frame(k127_567814[,2:12])
k127_233705_df2 <- data.frame(k127_233705[,2:12])
k127_173951_df2 <- data.frame(k127_173951[,2:12])

# turn the data frame above into the dna_seg object and then a string
# that genoPlotR uses 
dna_seg1.1 <- dna_seg(k127_162855_df2)
dna_segs.1 <- list(dna_seg1.1)

dna_seg1.2 <- dna_seg(k127_567814_df2)
dna_segs.2 <- list(dna_seg1.2)

dna_seg1.3 <- dna_seg(k127_233705_df2)
dna_segs.3 <- list(dna_seg1.3)

dna_seg1.4 <- dna_seg(k127_173951_df2)
dna_segs.4 <- list(dna_seg1.4)

plot_gene_map(dna_segs.1,
              comparisons = NULL,
              tree = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              legend = NULL,
              annotations = NULL, 
              annotation_height = 1, 
              annotation_cex = 0.8, 
              seg_plots=NULL,    # user-defined plots
              seg_plot_height=3, # height of plots (in lines)
              seg_plot_height_unit="lines", # unit of preceding
              seg_plot_yaxis=3, # if non-null or non false, ticks
              seg_plot_yaxis_cex=scale_cex,
              xlims = NULL,
              offsets = NULL,
              minimum_gap_size = 0.05,
              fixed_gap_length = FALSE,
              limit_to_longest_dna_seg = TRUE,
              main = NULL, 
              main_pos = "centre", 
              dna_seg_labels = NULL, 
              dna_seg_label_cex=1,
              dna_seg_label_col="black",
              gene_type = NULL,
              arrow_head_len = 200,
              dna_seg_line = TRUE,
              scale = TRUE, 
              dna_seg_scale = TRUE,
              n_scale_ticks=7,
              scale_cex=0.6,
              override_color_schemes = FALSE,
              plot_new=TRUE,
              debug = 0)

plot_gene_map(dna_segs.2,
              comparisons = NULL,
              tree = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              legend = NULL,
              annotations = NULL, 
              annotation_height = 1, 
              annotation_cex = 0.8, 
              seg_plots=NULL,    # user-defined plots
              seg_plot_height=3, # height of plots (in lines)
              seg_plot_height_unit="lines", # unit of preceding
              seg_plot_yaxis=3, # if non-null or non false, ticks
              seg_plot_yaxis_cex=scale_cex,
              xlims = NULL,
              offsets = NULL,
              minimum_gap_size = 0.05,
              fixed_gap_length = FALSE,
              limit_to_longest_dna_seg = TRUE,
              main = NULL, 
              main_pos = "centre", 
              dna_seg_labels = NULL, 
              dna_seg_label_cex=1,
              dna_seg_label_col="black",
              gene_type = NULL,
              arrow_head_len = 200,
              dna_seg_line = TRUE,
              scale = TRUE, 
              dna_seg_scale = TRUE,
              n_scale_ticks=7,
              scale_cex=0.6,
              override_color_schemes = FALSE,
              plot_new=TRUE,
              debug = 0)

plot_gene_map(dna_segs.3,
              comparisons = NULL,
              tree = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              legend = NULL,
              annotations = NULL, 
              annotation_height = 1, 
              annotation_cex = 0.8, 
              seg_plots=NULL,    # user-defined plots
              seg_plot_height=3, # height of plots (in lines)
              seg_plot_height_unit="lines", # unit of preceding
              seg_plot_yaxis=3, # if non-null or non false, ticks
              seg_plot_yaxis_cex=scale_cex,
              xlims = NULL,
              offsets = NULL,
              minimum_gap_size = 0.05,
              fixed_gap_length = FALSE,
              limit_to_longest_dna_seg = TRUE,
              main = NULL, 
              main_pos = "centre", 
              dna_seg_labels = NULL, 
              dna_seg_label_cex=1,
              dna_seg_label_col="black",
              gene_type = NULL,
              arrow_head_len = 200,
              dna_seg_line = TRUE,
              scale = TRUE, 
              dna_seg_scale = TRUE,
              n_scale_ticks=7,
              scale_cex=0.6,
              override_color_schemes = FALSE,
              plot_new=TRUE,
              debug = 0)

plot_gene_map(dna_segs.4,
              comparisons = NULL,
              tree = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              legend = NULL,
              annotations = NULL, 
              annotation_height = 1, 
              annotation_cex = 0.8, 
              seg_plots=NULL,    # user-defined plots
              seg_plot_height=3, # height of plots (in lines)
              seg_plot_height_unit="lines", # unit of preceding
              seg_plot_yaxis=3, # if non-null or non false, ticks
              seg_plot_yaxis_cex=scale_cex,
              xlims = NULL,
              offsets = NULL,
              minimum_gap_size = 0.05,
              fixed_gap_length = FALSE,
              limit_to_longest_dna_seg = TRUE,
              main = NULL, 
              main_pos = "centre", 
              dna_seg_labels = NULL, 
              dna_seg_label_cex=1,
              dna_seg_label_col="black",
              gene_type = NULL,
              arrow_head_len = 200,
              dna_seg_line = TRUE,
              scale = TRUE, 
              dna_seg_scale = TRUE,
              n_scale_ticks=7,
              scale_cex=0.6,
              override_color_schemes = FALSE,
              plot_new=TRUE,
              debug = 0)

############### Figure 2: Viral fractional abundance changes upon induction, and select temperate phages experience a substantial increase ###############
##### Figure 2b: Viral fractional abundance NMDS plot #####
Biofilm_virus_abundances <- read_excel("/path/to/Supplementak_DataSheet6.xlsx")

Biofilm_virus_abundances_wide <- Biofilm_virus_abundances %>% select(sample_type, sample, genome_id, FRACTIONAL_ABUNDANCE) %>% 
  pivot_wider(names_from = genome_id, values_from = FRACTIONAL_ABUNDANCE)

Biofilm_virus_abundances_wide2 <- Biofilm_virus_abundances %>% select(sample, genome_id, FRACTIONAL_ABUNDANCE) %>% 
  pivot_wider(names_from = sample, values_from = FRACTIONAL_ABUNDANCE)

matrix.biofilm.viruses <- as.matrix(Biofilm_virus_abundances_wide[,3:336]) 

dist.biofilm.viruses <-vegdist(matrix.biofilm.viruses, method='bray')
dist.biofilm.viruses

dif.biofilm.viruses <-adonis2(dist.biofilm.viruses~as.factor(Biofilm_virus_abundances_wide$sample_type), data=Biofilm_virus_abundances_wide, permutations=9999)
dif.biofilm.viruses

simper.biofilm.viruses <- simper(matrix.biofilm.viruses, Biofilm_virus_abundances_wide$sample_type, permutations = 9999)
simper.summary.biofilm.viruses <- summary(simper.biofilm.viruses, ordered = TRUE, digits = max(3,getOption("digits") - 3))
simper.summary.biofilm.viruses

simper_summary_biofilm_viruses_MC <- (simper.summary.biofilm.viruses[["MC-12PostInd_MC-PreInd"]])
simper_summary_biofilm_viruses_MC_CT <- (simper.summary.biofilm.viruses[["CT-12PostInd_MC-12PostInd"]])

simper_biofilm_viruses_MC <- simper_summary_biofilm_viruses_MC %>% filter(simper_summary_biofilm_viruses_MC$p < 0.05) %>% rownames_to_column(var = "Genome_ID") %>% select(Genome_ID)

simper_biofilm_viruses_MC_CT <- simper_summary_biofilm_viruses_MC_CT %>% filter(simper_summary_biofilm_viruses_MC_CT$p < 0.05) %>% rownames_to_column(var = "Genome_ID") %>% select(Genome_ID)

simper_biofilm_viruses <- rbind(simper_biofilm_viruses_MC, simper_biofilm_viruses_MC_CT)

nmds.2000 <- metaMDS(matrix.biofilm.viruses, permutations = 9999, distance = "bray")
nmds.2000
plot(nmds.2000)
goodness(nmds.2000)
stressplot(nmds.2000)

data.2000.scores.sample = as.data.frame(scores(nmds.2000)$sites)
data.2000.scores.sample
data.2000.scores.sample$sample_type = Biofilm_virus_abundances_wide$sample_type
data.2000.scores.sample

data.2000.scores.phage = as.data.frame(scores(nmds.2000)$species)
data.2000.scores.phage

plot_2000 <- ggplot() + 
  geom_point(data=data.2000.scores.sample,aes(x=NMDS1,y=NMDS2, colour =sample_type), size=3) + 
  stat_ellipse(data=data.2000.scores.sample,aes(x=NMDS1,y=NMDS2, colour =sample_type), level =0.95) +
  coord_equal() + theme_bw() + scale_colour_manual(values=c('CT-PreInd'="#b9e1cd",'MC-PreInd'="#e0dd6c",'CT-12PostInd'="#ba5619",'MC-12PostInd'="#2a061b"))
plot_2000
##### Figure 2c: Change in viral fractional abundance boxplots #####
viral_change_abundance <- Biofilm_virus_abundances_wide2 %>% summarize(
  genome_id = Biofilm_virus_abundances_wide2$genome_id,
  MC5_delta = ((Biofilm_virus_abundances_wide2$`MC-12PostInd-BF5` - Biofilm_virus_abundances_wide2$`MC-PreInd-BF5`)/Biofilm_virus_abundances_wide2$`MC-PreInd-BF5`)*100,
  MC4_delta = ((Biofilm_virus_abundances_wide2$`MC-12PostInd-BF4` - Biofilm_virus_abundances_wide2$`MC-PreInd-BF4`)/Biofilm_virus_abundances_wide2$`MC-PreInd-BF4`)*100,
  MC3_delta = ((Biofilm_virus_abundances_wide2$`MC-12PostInd-BF3` - Biofilm_virus_abundances_wide2$`MC-PreInd-BF3`)/Biofilm_virus_abundances_wide2$`MC-PreInd-BF3`)*100,
  MC2_delta = ((Biofilm_virus_abundances_wide2$`MC-12PostInd-BF2` - Biofilm_virus_abundances_wide2$`MC-PreInd-BF2`)/Biofilm_virus_abundances_wide2$`MC-PreInd-BF2`)*100,
  MC1_delta = ((Biofilm_virus_abundances_wide2$`MC-12PostInd-BF1` - Biofilm_virus_abundances_wide2$`MC-PreInd-BF1`)/Biofilm_virus_abundances_wide2$`MC-PreInd-BF1`)*100,
  CT5_delta = ((Biofilm_virus_abundances_wide2$`CT-12PostInd-BF5` - Biofilm_virus_abundances_wide2$`CT-PreInd-BF5`)/Biofilm_virus_abundances_wide2$`CT-PreInd-BF5`)*100,
  CT4_delta = ((Biofilm_virus_abundances_wide2$`CT-12PostInd-BF4` - Biofilm_virus_abundances_wide2$`CT-PreInd-BF4`)/Biofilm_virus_abundances_wide2$`CT-PreInd-BF4`)*100,
  CT3_delta = ((Biofilm_virus_abundances_wide2$`CT-12PostInd-BF3` - Biofilm_virus_abundances_wide2$`CT-PreInd-BF3`)/Biofilm_virus_abundances_wide2$`CT-PreInd-BF3`)*100,
  CT2_delta = ((Biofilm_virus_abundances_wide2$`CT-12PostInd-BF2` - Biofilm_virus_abundances_wide2$`CT-PreInd-BF2`)/Biofilm_virus_abundances_wide2$`CT-PreInd-BF2`)*100,
  CT1_delta = ((Biofilm_virus_abundances_wide2$`CT-12PostInd-BF1` - Biofilm_virus_abundances_wide2$`CT-PreInd-BF1`)/Biofilm_virus_abundances_wide2$`CT-PreInd-BF1`)*100
)

viral_change_abundance_simper <- simper_biofilm_viruses %>%
  left_join(viral_change_abundance %>% select(genome_id, MC5_delta, MC4_delta, MC3_delta, MC2_delta, 
                                              MC1_delta, CT5_delta, CT4_delta, CT3_delta, CT2_delta, CT1_delta), by = c("Genome_ID" = "genome_id"))

viral_change_abundance_simper_longer <- viral_change_abundance_simper %>% pivot_longer(!Genome_ID,  names_to = "replicate", values_to = "change_in_FA")

viral_change_abundance_simper_clean <- na.omit(viral_change_abundance_simper_longer)
viral_change_abundance_simper_clean2 <- viral_change_abundance_simper_clean %>% filter_all(all_vars(!is.infinite(.)))

MC_change_abundance_values <- viral_change_abundance_simper_clean2 %>% filter(replicate == "MC5_delta" | replicate == "MC4_delta"| replicate == "MC3_delta" |
                                                                                replicate == "MC2_delta" | replicate == "MC1_delta") %>% group_by(Genome_ID) %>% summarize(
                                                                                  sample_type = "Mitomycin_c",
                                                                                  n = n(), 
                                                                                  median = median(change_in_FA),
                                                                                  sd = sd(change_in_FA),
                                                                                  se = sd/sqrt(n),
                                                                                  median_plus_se = median + se,
                                                                                  median_minus_se = median - se)

CT_change_abundance_values <- viral_change_abundance_simper_clean2 %>% filter(replicate == "CT5_delta" | replicate == "CT4_delta"| replicate == "CT3_delta" |
                                                                                replicate == "CT2_delta" | replicate == "CT1_delta") %>% group_by(Genome_ID) %>% summarize(
                                                                                  sample_type = "No_treatment",
                                                                                  n = n(), 
                                                                                  median = median(change_in_FA),
                                                                                  sd = sd(change_in_FA),
                                                                                  se = sd/sqrt(n),
                                                                                  median_plus_se = median + se,
                                                                                  median_minus_se = median - se)

Combined_viral_median_change_abundance <- rbind(MC_change_abundance_values, CT_change_abundance_values)

# make sure difference change in mitomycin c samples is positive 
Change_abundnce_simper_phages_MC_positive_ids <- Combined_viral_median_change_abundance[Combined_viral_median_change_abundance[, 2] == "Mitomycin_c" & 
                                                                                          Combined_viral_median_change_abundance[, 4] > 0, ] %>% select(Genome_ID)

Change_abundance_phages_MC_pos <- Change_abundnce_simper_phages_MC_positive_ids %>%
  left_join(Combined_viral_median_change_abundance %>% select(Genome_ID,sample_type,n,median,sd,se,median_plus_se,median_minus_se), by = "Genome_ID") %>%
  na.omit(Change_abundance_phages_MC_pos)

# make sure standard error does not overlap 
# none of them overlap
se_range_overlap <- Change_abundance_phages_MC_pos %>%
  group_by(Genome_ID) %>%
  do({
    x <- .
    n <- nrow(x)
    overlaps <- matrix(FALSE, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          overlaps[i, j] <- x$median_minus_se[i] <= x$median_plus_se[j] &&
            x$median_minus_se[j] <= x$median_plus_se[i]
        }
      }
    }
    overlap_flag <- apply(overlaps, 1, any)
    x$overlap <- overlap_flag
    x
  }) %>%
  ungroup()
se_no_overlap <- se_range_overlap %>% group_by(Genome_ID) %>% filter(overlap == "FALSE") %>% select(Genome_ID)
se_no_overlap_derep <- se_no_overlap[!duplicated(se_no_overlap), ]

Phages_MC_pos_nose_overlap <- se_no_overlap_derep %>%
  left_join(Change_abundance_phages_MC_pos %>% select(Genome_ID,sample_type,n,median,sd,se,median_plus_se,median_minus_se), by = "Genome_ID")

Phages_MC_pos_nose_overlap_wide <- Phages_MC_pos_nose_overlap %>% select(Genome_ID, median, sample_type) %>% 
  pivot_wider(names_from = sample_type, values_from = median) %>% na.omit(Phages_MC_pos_nose_overlap)
head(Phages_MC_pos_nose_overlap_wide)

MC_greater_IDs <- as.data.frame(Phages_MC_pos_nose_overlap_wide$Genome_ID[Phages_MC_pos_nose_overlap_wide$Mitomycin_c > Phages_MC_pos_nose_overlap_wide$No_treatment])
MC_greater_IDs$Genome_ID <- MC_greater_IDs$`Phages_MC_pos_nose_overlap_wide$Genome_ID[Phages_MC_pos_nose_overlap_wide$Mitomycin_c > Phages_MC_pos_nose_overlap_wide$No_treatment]`

Phages_MC_pos_nose_overlap_great <- MC_greater_IDs %>%
  left_join(Phages_MC_pos_nose_overlap %>% select(Genome_ID,sample_type,n,median,sd,se,median_plus_se,median_minus_se), by = "Genome_ID")
Phages_MC_pos_nose_overlap_great$`Phages_MC_pos_nose_overlap_wide$Genome_ID[Phages_MC_pos_nose_overlap_wide$Mitomycin_c > Phages_MC_pos_nose_overlap_wide$No_treatment]` <- NULL

# get viruses that had infinite values introduced in the change in fractional abundance calculation
viruses_with_inf_values <- viral_change_abundance %>% filter(if_any(everything(), ~ is.infinite(.) | is.nan(.))) %>% select(genome_id)

Simper_phages_inf <- as.data.frame(intersect(unlist(viruses_with_inf_values), unlist(simper_biofilm_viruses)))

Simper_phages_inf$Genome_ID <- Simper_phages_inf$`intersect(unlist(viruses_with_inf_values), unlist(simper_biofilm_viruses))`

Simper_phages_inf$`intersect(unlist(viruses_with_inf_values), unlist(simper_biofilm_viruses))` <- NULL

Simper_phages_inf_values <- Simper_phages_inf %>%
  left_join(viral_change_abundance %>% select(genome_id, MC5_delta, MC4_delta, MC3_delta, MC2_delta, 
                                              MC1_delta, CT5_delta, CT4_delta, CT3_delta, CT2_delta, CT1_delta), by = c("Genome_ID" = "genome_id"))

Longer_simper_phages_inf <- Simper_phages_inf_values %>% pivot_longer(!Genome_ID,  names_to = "replicate", values_to = "change_in_FA")

clean_simper_phages_inf <- na.omit(Longer_simper_phages_inf) 
clean2_simper_phages_inf <- clean_simper_phages_inf %>% filter_all(all_vars(!is.infinite(.)))

MC_change_values_inf <- clean2_simper_phages_inf %>% filter(replicate == "MC5_delta" | replicate == "MC4_delta"| replicate == "MC3_delta" |
                                                              replicate == "MC2_delta" | replicate == "MC1_delta") %>% group_by(Genome_ID) %>% summarize(
                                                                sample_type = "Mitomycin_c",
                                                                n = n(), 
                                                                median = median(change_in_FA),
                                                                sd = sd(change_in_FA),
                                                                se = sd/sqrt(n),
                                                                median_plus_se = median + se,
                                                                median_minus_se = median - se)

CT_change_values_inf <- clean2_simper_phages_inf %>% filter(replicate == "CT5_delta" | replicate == "CT4_delta"| replicate == "CT3_delta" |
                                                              replicate == "CT2_delta" | replicate == "CT1_delta") %>% group_by(Genome_ID) %>% summarize(
                                                                sample_type = "No_treatment",
                                                                n = n(), 
                                                                median = median(change_in_FA),
                                                                sd = sd(change_in_FA),
                                                                se = sd/sqrt(n),
                                                                median_plus_se = median + se,
                                                                median_minus_se = median - se)

MC_change_values_inf_filt <- MC_change_values_inf %>% filter(!grepl("MC-12PostInd-BF1_megahit_k127_312308_flag_1_multi_13.0000_len_19863", Genome_ID)) %>% 
  filter(!grepl("CT-12PostInd-BF1_megahit_k127_195199_flag_0_multi_50.9236_len_13219", Genome_ID)) 
CT_change_values_inf_filt <- CT_change_values_inf %>% filter(!grepl("MC-12PostInd-BF1_megahit_k127_312308_flag_1_multi_13.0000_len_19863", Genome_ID)) %>% 
  filter(!grepl("CT-12PostInd-BF1_megahit_k127_195199_flag_0_multi_50.9236_len_13219", Genome_ID)) 

Combined_median_change_simper_phages_inf <- as.data.frame(rbind(MC_change_values_inf_filt, CT_change_values_inf_filt))

All_combined_simper_phages <- rbind(Combined_median_change_simper_phages_inf,Phages_MC_pos_nose_overlap_great)
All_combined_simper_phages$Genome_ID <- gsub("_flag.*", "", All_combined_simper_phages$Genome_ID)
All_combined_simper_phages$Genome_ID <- gsub("_megahit", "", All_combined_simper_phages$Genome_ID)

mycolors2 <- (c("#ba5619","#b9e1cd"))
ggplot(All_combined_simper_phages, aes(fill =sample_type, x=median, y= reorder(Genome_ID, -median))) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(x=median, y=Genome_ID, 
                    xmin=median-se, xmax=median+se), position=position_dodge(.9), width = 0.5) +
  scale_fill_manual(values= mycolors2) + theme_bw() 

##### Figure 2d: Viral genome plots #####
#seperate into the different genomes 
k127_137292 <- phage_annotations %>% filter(name == "MC-12PostInd-BF3_megahit_k127_137292_flag_0_multi_26.9193_len_3350_extended_partial")
k127_317508 <- phage_annotations %>% filter(name == "MC-12PostInd-BF4_megahit_k127_317508_flag_0_multi_92.9303_len_13319_extended_partial")
k127_564127 <- phage_annotations %>% filter(name == "MC-12PostInd-BF5_megahit_k127_564127_flag_0_multi_38.2793_len_2773_extended_partial")
k127_348120 <- phage_annotations %>% filter(name == "MC-12PostInd-BF1_megahit_k127_348120_flag_0_multi_14.9672_len_2233_extended_partial")

# pull out the coumns necessary for genoPlotR 
# name, start, end, strand, col, fill, lty, lwd, pch, cex, gene_type
# and convert them to a data frame
k127_137292_df <- data.frame(k127_137292[,2:12])
k127_317508_df <- data.frame(k127_317508[,2:12])
k127_564127_df <- data.frame(k127_564127[,2:12])
k127_348120_df <- data.frame(k127_348120[,2:12])

dna_seg1.5 <- dna_seg(k127_137292_df)
dna_segs.5 <- list(dna_seg1.5)

dna_seg1.6 <- dna_seg(k127_317508_df)
dna_segs.6 <- list(dna_seg1.6)

dna_seg1.7 <- dna_seg(k127_564127_df)
dna_segs.7 <- list(dna_seg1.7)

dna_seg1.8 <- dna_seg(k127_348120_df)
dna_segs.8 <- list(dna_seg1.8)

plot_gene_map(dna_segs.5,
              comparisons = NULL,
              tree = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              legend = NULL,
              annotations = NULL, 
              annotation_height = 1, 
              annotation_cex = 0.8, 
              seg_plots=NULL,    # user-defined plots
              seg_plot_height=3, # height of plots (in lines)
              seg_plot_height_unit="lines", # unit of preceding
              seg_plot_yaxis=3, # if non-null or non false, ticks
              seg_plot_yaxis_cex=scale_cex,
              xlims = NULL,
              offsets = NULL,
              minimum_gap_size = 0.05,
              fixed_gap_length = FALSE,
              limit_to_longest_dna_seg = TRUE,
              main = NULL, 
              main_pos = "centre", 
              dna_seg_labels = NULL, 
              dna_seg_label_cex=1,
              dna_seg_label_col="black",
              gene_type = NULL,
              arrow_head_len = 200,
              dna_seg_line = TRUE,
              scale = TRUE, 
              dna_seg_scale = TRUE,
              n_scale_ticks=7,
              scale_cex=0.6,
              override_color_schemes = FALSE,
              plot_new=TRUE,
              debug = 0)

plot_gene_map(dna_segs.6,
              comparisons = NULL,
              tree = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              legend = NULL,
              annotations = NULL, 
              annotation_height = 1, 
              annotation_cex = 0.8, 
              seg_plots=NULL,    # user-defined plots
              seg_plot_height=3, # height of plots (in lines)
              seg_plot_height_unit="lines", # unit of preceding
              seg_plot_yaxis=3, # if non-null or non false, ticks
              seg_plot_yaxis_cex=scale_cex,
              xlims = NULL,
              offsets = NULL,
              minimum_gap_size = 0.05,
              fixed_gap_length = FALSE,
              limit_to_longest_dna_seg = TRUE,
              main = NULL, 
              main_pos = "centre", 
              dna_seg_labels = NULL, 
              dna_seg_label_cex=1,
              dna_seg_label_col="black",
              gene_type = NULL,
              arrow_head_len = 200,
              dna_seg_line = TRUE,
              scale = TRUE, 
              dna_seg_scale = TRUE,
              n_scale_ticks=7,
              scale_cex=0.6,
              override_color_schemes = FALSE,
              plot_new=TRUE,
              debug = 0)

plot_gene_map(dna_segs.7,
              comparisons = NULL,
              tree = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              legend = NULL,
              annotations = NULL, 
              annotation_height = 1, 
              annotation_cex = 0.8, 
              seg_plots=NULL,    # user-defined plots
              seg_plot_height=3, # height of plots (in lines)
              seg_plot_height_unit="lines", # unit of preceding
              seg_plot_yaxis=3, # if non-null or non false, ticks
              seg_plot_yaxis_cex=scale_cex,
              xlims = NULL,
              offsets = NULL,
              minimum_gap_size = 0.05,
              fixed_gap_length = FALSE,
              limit_to_longest_dna_seg = TRUE,
              main = NULL, 
              main_pos = "centre", 
              dna_seg_labels = NULL, 
              dna_seg_label_cex=1,
              dna_seg_label_col="black",
              gene_type = NULL,
              arrow_head_len = 200,
              dna_seg_line = TRUE,
              scale = TRUE, 
              dna_seg_scale = TRUE,
              n_scale_ticks=7,
              scale_cex=0.6,
              override_color_schemes = FALSE,
              plot_new=TRUE,
              debug = 0)

plot_gene_map(dna_segs.8,
              comparisons = NULL,
              tree = NULL,
              tree_branch_labels_cex = NULL,
              tree_scale = FALSE,
              legend = NULL,
              annotations = NULL, 
              annotation_height = 1, 
              annotation_cex = 0.8, 
              seg_plots=NULL,    # user-defined plots
              seg_plot_height=3, # height of plots (in lines)
              seg_plot_height_unit="lines", # unit of preceding
              seg_plot_yaxis=3, # if non-null or non false, ticks
              seg_plot_yaxis_cex=scale_cex,
              xlims = NULL,
              offsets = NULL,
              minimum_gap_size = 0.05,
              fixed_gap_length = FALSE,
              limit_to_longest_dna_seg = TRUE,
              main = NULL, 
              main_pos = "centre", 
              dna_seg_labels = NULL, 
              dna_seg_label_cex=1,
              dna_seg_label_col="black",
              gene_type = NULL,
              arrow_head_len = 200,
              dna_seg_line = TRUE,
              scale = TRUE, 
              dna_seg_scale = TRUE,
              n_scale_ticks=7,
              scale_cex=0.6,
              override_color_schemes = FALSE,
              plot_new=TRUE,
              debug = 0)

############### Figure 3: Bacterial community composition shifts upon prophage induction  ###############
#######  bMAG relative abundance based figures (Supplemental Figure 1, Figure 3, and Supplemental Figure 2) ####### 
bMAG_relative_abundance <- read_excel('/path/to/Supplemental_DataSheet2.xlsx')
bMAG_relative_abundance_nobin11 <- bMAG_relative_abundance[bMAG_relative_abundance$Genomic_bins != "CT-12PostInd_bin.11", ]

bMAG_relative_abundance_long <- bMAG_relative_abundance %>% 
  select(Genomic_bins, `MC-12PostInd-BF1`,`MC-12PostInd-BF2`,`MC-12PostInd-BF3`,`MC-12PostInd-BF4`, `MC-12PostInd-BF5`,
         `CT-12PostInd-BF1`,`CT-12PostInd-BF2`, `CT-12PostInd-BF3`, `CT-12PostInd-BF4`, `CT-12PostInd-BF5`,
         `MC-PreInd-BF1`, `MC-PreInd-BF2`, `MC-PreInd-BF3`, `MC-PreInd-BF4`, `MC-PreInd-BF5`,
         `CT-PreInd-BF1`, `CT-PreInd-BF2`, `CT-PreInd-BF3`, `CT-PreInd-BF4`, `CT-PreInd-BF5`) %>% 
  pivot_longer(!Genomic_bins,  names_to = "replicate", values_to = "genomes_per_million_reads")

bMAG_relative_abundance_nobin11_long <- bMAG_relative_abundance_nobin11 %>% 
  select(Genomic_bins, `MC-12PostInd-BF1`,`MC-12PostInd-BF2`,`MC-12PostInd-BF3`,`MC-12PostInd-BF4`, `MC-12PostInd-BF5`,
         `CT-12PostInd-BF1`,`CT-12PostInd-BF2`, `CT-12PostInd-BF3`, `CT-12PostInd-BF4`, `CT-12PostInd-BF5`,
         `MC-PreInd-BF1`, `MC-PreInd-BF2`, `MC-PreInd-BF3`, `MC-PreInd-BF4`, `MC-PreInd-BF5`,
         `CT-PreInd-BF1`, `CT-PreInd-BF2`, `CT-PreInd-BF3`, `CT-PreInd-BF4`, `CT-PreInd-BF5`) %>% 
  pivot_longer(!Genomic_bins,  names_to = "replicate", values_to = "genomes_per_million_reads")

bMAG_relative_abundance_nobin11_long$sample_type <- sub("^([^-]+-[^-]+)-.*$", "\\1", bMAG_relative_abundance_nobin11_long$replicate)

bMAG_relative_abundance_flip <- bMAG_relative_abundance_long %>% pivot_wider(names_from = Genomic_bins, values_from = genomes_per_million_reads)
head(bMAG_relative_abundance_flip)
bMAG_relative_abundance_flip_nobin11 <- bMAG_relative_abundance_nobin11_long %>% pivot_wider(names_from = Genomic_bins, values_from = genomes_per_million_reads)
bMAG_relative_abundance_flip_nobin11$sample_type <- sub("^([^-]+-[^-]+)-.*$", "\\1", bMAG_relative_abundance_flip_nobin11$replicate)
head(bMAG_relative_abundance_flip_nobin11)

#######  Supplemental Figure 1 ####### 
#NMDS matrix for all bMAGs including CT-PostInd-bin.11 (Salmonella enterica)
NMDS_matrix <- as.matrix(bMAG_relative_abundance_flip[,2:30])
head(NMDS_matrix)
nmds.raw <- metaMDS(NMDS_matrix, permutations = 9999, distance = "bray", k = 2 ,maxit = 999, trymax = 500,wascores = TRUE)
nmds.raw
plot(nmds.raw)
goodness(nmds.raw)
stressplot(nmds.raw)
#NMDS matrix for exlcuding CT-PostInd-bin.11
NMDS_matrix_nobin11 <- as.matrix(bMAG_relative_abundance_flip_nobin11[,3:30])
head(NMDS_matrix_nobin11)
nmds.nobin11 <- metaMDS(NMDS_matrix_nobin11, permutations = 9999, distance = "bray", k = 2 ,maxit = 999, trymax = 500,wascores = TRUE)
nmds.nobin11
plot(nmds.nobin11)
goodness(nmds.nobin11)
stressplot(nmds.nobin11)

##### Figure 3a: bMAG fractional abundance NMDS plot #####
NMDS_matrix_flip <- as.matrix(bMAG_relative_abundance_nobin11[,9:28])
head(NMDS_matrix_flip)
nmds <- metaMDS(NMDS_matrix_flip, permutations = 9999, distance = "bray", k = 2 ,maxit = 999, trymax = 500,wascores = TRUE)
nmds
plot(nmds)
goodness(nmds)
stressplot(nmds)

# extract data scores for samples and bMAGs
data.scores.MAG = as.data.frame(scores(nmds)$sites)
bMAG_relative_abundance_nobin11$genus <- sub("(\\w+).*", "\\1", bMAG_relative_abundance_nobin11$GTDB_tk_taxonomy)
data.scores.MAG$Genomic_bins = bMAG_relative_abundance_nobin11$Genomic_bins
data.scores.MAG$Genus = bMAG_relative_abundance_nobin11$genus


sample_types_list <- c('MC-PostInd', 'MC-PostInd', 'MC-PostInd', 'MC-PostInd', 'MC-PostInd',
                       'CT-PostInd', 'CT-PostInd', 'CT-PostInd', 'CT-PostInd', 'CT-PostInd',
                       'MC-PreInd', 'MC-PreInd','MC-PreInd', 'MC-PreInd','MC-PreInd',
                       'CT-PreInd', 'CT-PreInd','CT-PreInd', 'CT-PreInd','CT-PreInd')

sample_types <- data.frame(sample_type = sample_types_list)
data.scores.sample = as.data.frame(scores(nmds)$species)
data.scores.sample$sample_type = sample_types$sample_type

# make a prettier NMDS plot :)
plot <- ggplot() + 
  geom_point(data=data.scores.MAG,aes(x=NMDS1,y=NMDS2,shape=Genus),size =2) +
  scale_shape_manual(values = c(0, 2, 4, 5, 8)) +
  geom_point(data=data.scores.sample,aes(x=NMDS1,y=NMDS2, colour =sample_type), size=3) + 
  stat_ellipse(data=data.scores.sample,aes(x=NMDS1,y=NMDS2, colour =sample_type), level =0.95) +
  coord_equal() + theme_bw() + scale_colour_manual(values=c('CT-PreInd'="#b9e1cd",'MC-PreInd'="#e0dd6c",'CT-PostInd'="#ba5619",'MC-PostInd'="#2a061b"))
plot

# PERMANOVA for bMAG relative abundance data 
bMAG_PERM_matrix <- as.matrix(bMAG_relative_abundance_flip_nobin11[,3:30])

bMAG.PERM.dist<-vegdist(bMAG_PERM_matrix, method='bray')
bMAG.PERM.dist

bMAG.PERM.dif <-adonis2(bMAG.PERM.dist~as.factor(bMAG_relative_abundance_flip_nobin11$sample_type), data=bMAG_relative_abundance_flip_nobin11, permutations=9999)
bMAG.PERM.dif

# Pairwise ADONIS for bMAG relative abundance data 
pairwise.adonis2(bMAG.PERM.dist~sample_type, data=bMAG_relative_abundance_flip_nobin11, permutations=9999, methods='bray')

# SIMPER analysis to identify bMAGs driving differences in samples
bMAG.simper <- simper(bMAG_PERM_matrix, bMAG_relative_abundance_flip_nobin11$sample_type, permutations = 9999, methods='bray')
bMAG.simper.summary <- summary(bMAG.simper, ordered = TRUE, digits = max(6,getOption("digits") - 6))

simper_summary_MC_MC <- (bMAG.simper.summary[["MC-12PostInd_MC-PreInd"]])
simper_summary_MC_MC <- simper_summary_MC_MC[simper_summary_MC_MC[, 7] < 0.05, ]

simper_summary_MC_CT <- (bMAG.simper.summary[["MC-12PostInd_CT-12PostInd"]])
simper_summary_MC_CT <- simper_summary_MC_CT[simper_summary_MC_CT[, 7] < 0.05, ]

####### Figure 3b: Change in growth of SIMPER identified bMAGs ####### 
# Calculate change in growth for bMAGS
bMAG_change_growth <- bMAG_relative_abundance_nobin11 %>% summarize(
  Genomic_bins = bMAG_relative_abundance_nobin11$Genomic_bins,
  MC5_delta = ((bMAG_relative_abundance_nobin11$`MC-12PostInd-BF5` - bMAG_relative_abundance_nobin11$`MC-PreInd-BF5`)/bMAG_relative_abundance_nobin11$`MC-PreInd-BF5`)*100,
  MC4_delta = ((bMAG_relative_abundance_nobin11$`MC-12PostInd-BF4` - bMAG_relative_abundance_nobin11$`MC-PreInd-BF4`)/bMAG_relative_abundance_nobin11$`MC-PreInd-BF4`)*100,
  MC3_delta = ((bMAG_relative_abundance_nobin11$`MC-12PostInd-BF3` - bMAG_relative_abundance_nobin11$`MC-PreInd-BF3`)/bMAG_relative_abundance_nobin11$`MC-PreInd-BF3`)*100,
  MC2_delta = ((bMAG_relative_abundance_nobin11$`MC-12PostInd-BF2` - bMAG_relative_abundance_nobin11$`MC-PreInd-BF2`)/bMAG_relative_abundance_nobin11$`MC-PreInd-BF2`)*100,
  MC1_delta = ((bMAG_relative_abundance_nobin11$`MC-12PostInd-BF1` - bMAG_relative_abundance_nobin11$`MC-PreInd-BF1`)/bMAG_relative_abundance_nobin11$`MC-PreInd-BF1`)*100,
  CT5_delta = ((bMAG_relative_abundance_nobin11$`CT-12PostInd-BF5` - bMAG_relative_abundance_nobin11$`CT-PreInd-BF5`)/bMAG_relative_abundance_nobin11$`CT-PreInd-BF5`)*100,
  CT4_delta = ((bMAG_relative_abundance_nobin11$`CT-12PostInd-BF4` - bMAG_relative_abundance_nobin11$`CT-PreInd-BF4`)/bMAG_relative_abundance_nobin11$`CT-PreInd-BF4`)*100,
  CT3_delta = ((bMAG_relative_abundance_nobin11$`CT-12PostInd-BF3` - bMAG_relative_abundance_nobin11$`CT-PreInd-BF3`)/bMAG_relative_abundance_nobin11$`CT-PreInd-BF3`)*100,
  CT2_delta = ((bMAG_relative_abundance_nobin11$`CT-12PostInd-BF2` - bMAG_relative_abundance_nobin11$`CT-PreInd-BF2`)/bMAG_relative_abundance_nobin11$`CT-PreInd-BF2`)*100,
  CT1_delta = ((bMAG_relative_abundance_nobin11$`CT-12PostInd-BF1` - bMAG_relative_abundance_nobin11$`CT-PreInd-BF1`)/bMAG_relative_abundance_nobin11$`CT-PreInd-BF1`)*100
)

bMAG_change_growth_longer <- bMAG_change_growth %>% pivot_longer(!Genomic_bins,  names_to = "replicate", values_to = "change_relative_abundance")

bMAG_change_growth_longer <- bMAG_change_growth %>% pivot_longer(!Genomic_bins,  names_to = "replicate", values_to = "delta_relative_abundance")

MC_change_growth <- bMAG_change_growth_longer %>% filter(replicate == "MC5_delta" | replicate == "MC4_delta"| replicate == "MC3_delta" |
  replicate == "MC2_delta" | replicate == "MC1_delta") %>% group_by(Genomic_bins) %>% summarize(
      sample_type = "MC",
      n = n(), 
      median = median(delta_relative_abundance),
      sd = sd(delta_relative_abundance),
      se = sd/sqrt(n))

CT_change_growth <- bMAG_change_growth_longer %>% filter(replicate == "CT5_delta" | replicate == "CT4_delta"| replicate == "CT3_delta" |
                    replicate == "CT2_delta" | replicate == "CT1_delta") %>% group_by(Genomic_bins) %>% summarize(
                        sample_type = "CT",
                        n = n(), 
                        median = median(delta_relative_abundance),
                        sd = sd(delta_relative_abundance),
                        se = sd/sqrt(n))

combined_median_change_growth <- rbind(MC_change_growth, CT_change_growth)

combined_median_change_SIMPR <- combined_median_change_growth %>% filter(Genomic_bins == 'MC-12PostInd_bin.4' | Genomic_bins == 'MC-12PostInd_bin.1' |
                                Genomic_bins == 'CT-PreInd_bin.6' | Genomic_bins == 'CT-12PostInd_bin.6' | Genomic_bins == 'MC-PreInd_bin.5' |
                                Genomic_bins == 'MC-PreInd_bin.2' | Genomic_bins == 'MC-PreInd_bin.1' | Genomic_bins == 'CT-PreInd_bin.5' |
                                Genomic_bins == 'CT-PreInd_bin.4' | Genomic_bins == 'CT-PreInd_bin.3' | Genomic_bins == 'MC-12PostInd_bin.5' |
                                Genomic_bins == 'MC-12PostInd_bin.2' | Genomic_bins == 'MC-PreInd_bin.4' | Genomic_bins == 'CT-12PostInd_bin.8' |
                                Genomic_bins == 'CT-12PostInd_bin.7' | Genomic_bins == 'CT-12PostInd_bin.3' | Genomic_bins == 'CT-12PostInd_bin.1' | 
                                Genomic_bins == 'CT-12PostInd_bin.5')

combined_median_change_SIMPR$sample_type <- factor(combined_median_change_SIMPR$sample_type, levels = c("CT", "MC"))
combined_median_change_SIMPR$sample_type <- factor(combined_median_change_SIMPR$sample_type, levels = rev(levels(combined_median_change_SIMPR$sample_type)))

mycolors2 <- (c("#b9e1cd","#ba5619"))
ggplot(combined_median_change_SIMPR, aes(fill =sample_type, x=median, y= 
  factor(Genomic_bins, level=c('MC-12PostInd_bin.4','MC-12PostInd_bin.1','CT-PreInd_bin.6','CT-12PostInd_bin.6',
  'MC-PreInd_bin.5', 'MC-PreInd_bin.2', 'MC-PreInd_bin.1','CT-PreInd_bin.5', 'CT-PreInd_bin.4',
  'CT-PreInd_bin.3','MC-12PostInd_bin.5','MC-12PostInd_bin.2','MC-PreInd_bin.4','CT-12PostInd_bin.8',
  'CT-12PostInd_bin.7','CT-12PostInd_bin.3','CT-12PostInd_bin.1', 'CT-12PostInd_bin.5')))) + 
  geom_bar(position="dodge", stat="identity") + geom_errorbar(aes(x=median, y=Genomic_bins, xmin=median-se, xmax=median+se), 
  position=position_dodge(.9), width = 0.5) + scale_fill_manual(values= mycolors2) + theme_bw() 

###### Sup. Fig. 2: bMAG relative abundance boxplots ##### 
mycolors4 <- (c("#b9e1cd","#ba5619","#e0dd6c", "#2a061b"))
SIMPER_bMAGs <- filter(bMAG_relative_abundance_nobin11_long, Genomic_bins %in% c(
  'CT-12PostInd_bin.8','CT-12PostInd_bin.6','CT-12PostInd_bin.3','MC-PreInd_bin.2','CT-PreInd_bin.4','MC-12PostInd_bin.5',
  'CT-12PostInd_bin.5','MC-12PostInd_bin.4','CT-12PostInd_bin.7','CT-12PostInd_bin.1','CT-PreInd_bin.5','MC-PreInd_bin.5',
   'MC-12PostInd_bin.2','MC-12PostInd_bin.1','CT-PreInd_bin.6','MC-PreInd_bin.1','CT-PreInd_bin.3','MC-PreInd_bin.4'))

SIMPER_boxplot <- ggplot(SIMPER_bMAGs, aes(x=sample_type, y=genomes_per_million_reads, 
  colour =factor(sample_type, level=c('CT-PreInd','CT-12PostInd','MC-PreInd','MC-12PostInd')))) + 
  scale_x_discrete(limits=c('CT-PreInd','CT-12PostInd','MC-PreInd','MC-12PostInd')) + geom_boxplot() + 
  theme_bw() + scale_colour_manual(values=mycolors4) + scale_y_continuous(trans ='log10') + removeGrid(x = TRUE) +
  stat_compare_means(aes(group = sample_type),comparisons = 
  list(c("MC-PreInd", "MC-12PostInd"), c("MC-12PostInd", "CT-12PostInd"), c("CT-PreInd", "MC-PreInd")), 
  method = "t.test", label = "p.signif") + facet_wrap(~Genomic_bins) 
SIMPER_boxplot

############### Figure 4: Prophage and host encoded AphA are divergent but have high structural similarity ###############
###### Figure 4a: AphA protein phylogeny tree #####
# The protein alignment and protein phylogeny tree was carried out using the online user interface version of MAFFT v7 
# Thee tree newick file was uploaded to iTOL for visualization, annotations, and downloaded for minor aesethic changes 
# in Adobe Illustrator

###### Figure 4b: Structural alignment of prophage and host encoded AphA #####
# Protein files for AphA from both prophage.8 and it's host CT-PostInd-bin.5 were uploaded to 
# AlphaFold Server all Alpha fold files for both structures were downloaded the protein structure files (.cif)
# Were uploaded to locally installed ChimeraX were the matchmaker command was used to align the two proteins
# AlphaFold Server: https://alphafoldserver.com/welcome
# ChimeraX Documentation: https://www.rbvi.ucsf.edu/chimerax/

###### Supplemental Figure 3: Structures for prophage.8 and CT-PostInd-bin.5 Apha proteins and AlphaFold #####
##### Supplemental Figure 3a: AphA structures (host = left and prophage = right) #####
##### Supplemental Figure 3b: Predicted aligned error for both AphA structures #####
# structures files (.cif) and predicted aligned error score files (.json) were uploaded and visualized with 
# ChimeraX. Structure colors represent pLDDT scores 

############### Figure 5: Bacterial metabolic capacity is altered post-induction ###############
pathways_full <- read_excel('/path/to/Supplemental_DataSheet7.xlsx')
head(pathways_full)

pathways_pathwise_great_70 <- filter(pathways_full, (pathwise_module_completeness>=0.70))
head(pathways_pathwise_great_70)
##### Figure 5a: bMAG % ANI tree with metabolism submodules present in bMAGs #####

pathways_submodules <- pathways_pathwise_great_70 %>% select(Genomic_bins, module_subcategory) %>% 
  group_by(Genomic_bins, module_subcategory) %>% summarize(
  module_subcategory_occurence =n())

submodules_wide <- pathways_submodules %>% select(Genomic_bins, module_subcategory, module_subcategory_occurence) %>% 
  pivot_wider(names_from = module_subcategory, values_from = module_subcategory_occurence)

submodules_wide_nona <- submodules_wide %>% replace(is.na(.), 0)
submodules_wide_presence <- submodules_wide_nona %>%
  mutate(across(where(is.numeric), ~ifelse(. >= 1, 1, .)))

# # table uploaded to iTOL to annotate bMAG ANI tree with submodule presence/absence data 
# write.table(submodules_wide_presence, file='/path/to/outpuy/Biofilm_submodule_bMAG_presence.tsv', quote=FALSE, sep='\t')

##### Figure 5b: Metabolism modules driving differences in post-induction samples #####
pathways_modules <- pathways_pathwise_great_70 %>% select(Genomic_bins, module_name)

pathway_modules_RA <- merge(pathways_modules, bMAG_relative_abundance_nobin11_long, by = "Genomic_bins", allow.cartesian = TRUE)

pathway_modules_RA_sums <- pathway_modules_RA %>% group_by(module_name, replicate) %>% summarize(
  module_sum = sum(genomes_per_million_reads)
)

modules_RA_wide <- pathway_modules_RA_sums %>% select(replicate, module_sum, module_name) %>% 
  pivot_wider(names_from = module_name, values_from = module_sum)

modules_RA_wide$sample_type <- str_extract(modules_RA_wide$replicate, "^[^-]+-[^-]+")

module_matrix <- as.matrix(modules_RA_wide[,2:133])

module.dist<-vegdist(module_matrix, method='bray') # breaks the dissimilarity matrix for some reason 

# PERMANOVA
submodule.dif <-adonis2(module.dist~as.factor(modules_RA_wide$sample_type), data=modules_RA_wide, permutations=9999)
submodule.dif

# Pairwise comparison
pairwise.adonis2(module.dist~sample_type, data=modules_RA_wide, permutations=9999, methods='bray')

#SIMPER
module.simper <- simper(module_matrix, modules_RA_wide$sample_type, permutations = 9999, methods='bray')
module.simper.summary <- summary(module.simper, ordered = TRUE, digits = max(3,getOption("digits") - 3))

simper_summary_MC_CT <- (module.simper.summary[["CT-12PostInd_MC-12PostInd"]])

simper_modules_MC_CT <- simper_summary_MC_CT %>% filter(simper_summary_MC_CT$p < 0.05) %>% rownames_to_column(var = "module_name") %>% select(module_name)

simper_modules_MC_CT_RA <- simper_modules_MC_CT %>%
  left_join(pathway_modules_RA_sums %>% select(module_name, replicate, module_sum), by = "module_name")

modules_MC_CT_RA_wide_simper <- simper_modules_MC_CT_RA %>% select(replicate, module_sum, module_name) %>% 
  pivot_wider(names_from = module_name, values_from = module_sum)

modules_simper_matrix <- as.matrix(modules_MC_CT_RA_wide_simper[,2:18])

rownames(modules_simper_matrix) = sapply(modules_MC_CT_RA_wide_simper$replicate,function(x)
  strsplit(as.character(x),split = "\\\\")[[1]][1])

pheatmap(mat = t(modules_simper_matrix), scale = "row",
         fontsize = 8, angle_col = 45, cellheight = 7, cellwidth = 7,
         cluster_rows=T, cluster_col=T, na_col = 'white',
         color=colorRampPalette(c("#cefae4","#e0dd6c","#ba5619","#2a061b"))(50))


############### Figure 6: Sargassum-derived biofilm community members are present in situ  ###############
##### Figure 6b: Biofilm induction temperate viruses are more abundant in situ  #####
Biofilm_viruses_2021_metagenomes <- read_excel ("/path/to/Supplemental_DataSheet8.xlsx")

Biofilm_viruses_abudance_medians <- Biofilm_viruses_abudance_info %>% group_by(Prediction, genome_id) %>% summarize(
  genome_median = median(FRACTIONAL_ABUNDANCE))

ggplot(Biofilm_viruses_abudance_medians, aes(x=genome_median, y=Prediction, fill = Prediction)) + 
  geom_violin() + geom_boxplot(width=0.1, fill = 'white') + scale_x_continuous(trans = 'log10') +
  scale_fill_manual(values= c("#ba5619", "#90AFA0")) + theme_bw()

wilcox.test(FRACTIONAL_ABUNDANCE ~ Prediction, data = Biofilm_viruses_abudance_info)

##### Figure 6c: Biofilm induction bMAGs are present in situ #####
bMAG_2021_metagenomes <- read_excel("/path/to/Supplemental_DataSheet9.xlsx")

bMAG_2021_metagenomes_longer <- bMAG_2021_metagenomes %>% select(Genomic_bins, `Sargasso-1`, `Sargasso-2`,
  `Sargasso-3`, `Sargasso-4`, `Sargasso-6`,`Seawater-Sterivex1`, `Seawater-Sterivex2`, `Seawater-Sterivex3`, 
  `Seawater-Sterivex4`, `Seawater-Sterivex5`) %>% 
  pivot_longer(!Genomic_bins,  names_to = "replicate", values_to = "genomes_per_million_reads")

bMAG_2021_metagenomes_longer$sample_type <- gsub("-.*", "", bMAG_2021_metagenomes_longer$replicate)

bMAG_2021_metagenome_stats <- bMAG_2021_metagenomes_longer %>% group_by(Genomic_bins, sample_type) %>% summarize(
  n = n(),
  xbar = mean(genomes_per_million_reads),
  median = median(genomes_per_million_reads),
  sd = sd(genomes_per_million_reads),
  se = sd/sqrt(n))

bMAG_2021_metagenome_stats_tax <- bMAG_2021_metagenome_stats %>%
  left_join(bMAG_2021_metagenomes %>% select(Genomic_bins, gtdb_tk_taxonomy), by = c("Genomic_bins" = "Genomic_bins"))

Sargassum <- bMAG_2021_metagenome_stats_tax %>% filter(sample_type == "Sargasso") 
Seawater <- bMAG_2021_metagenome_stats_tax %>% filter(sample_type == "Seawater") 

Sargassum_bMAG_tax_filt <- Sargassum %>% filter(xbar > 0.5)

Sargassum_05_bMAG_IDs <- Sargassum_bMAG_tax_filt %>% select(Genomic_bins)

Seawater_bMAG_tax_filt <- merge(Seawater, Sargassum_05_bMAG_IDs, by = "Genomic_bins") %>% 
  select(Genomic_bins, sample_type, n, xbar, median, sd, se, gtdb_tk_taxonomy)

mycolors7 = c('#948816',"#e0dd6c",'#ba5619','#755e6b', "#90AFA0",'#cefae4')

ggplot(Sargassum_bMAG_tax_filt, aes(x=xbar, y=reorder(Genomic_bins, -xbar), fill = gtdb_tk_taxonomy)) + 
  geom_bar(stat = "identity") + geom_errorbar(aes(x=xbar, y=Genomic_bins, xmin=xbar-se, xmax=xbar+se), 
  position=position_dodge(.9), width = 0.5) + theme_bw() + scale_fill_manual(values= mycolors7) + xlim(-0.5, 30)

ggplot(Seawater_bMAG_tax_filt, aes(x=xbar, y=factor(Genomic_bins, levels = c('MC-PreInd_bin.2','CT-PreInd_bin.4',
  'CT-12PostInd_bin.3', 'CT-PreInd_bin.8','CT-12PostInd_bin.10', 'CT-PreInd_bin.6', 'MC-12PostInd_bin.1',
  'MC-PreInd_bin.3','CT-12PostInd_bin.1','CT-PreInd_bin.2','MC-12PostInd_bin.5','MC-PreInd_bin.8',
  'MC-12PostInd_bin.7','CT-PreInd_bin.3', 'MC-PreInd_bin.1', 'CT-12PostInd_bin.5','CT-PreInd_bin.7')), fill = gtdb_tk_taxonomy)) + 
  geom_bar(stat = "identity") + geom_errorbar(aes(x=xbar, y=Genomic_bins, xmin=xbar-se, xmax=xbar+se), 
  position=position_dodge(.9), width = 0.5) + theme_bw() + scale_fill_manual(values= mycolors7) + xlim(-0.5, 30)


VHR.simper.summary