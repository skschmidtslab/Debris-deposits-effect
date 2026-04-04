#Processing of 16S data

#set the working directory
setwd("~/Documents/Steve Lab/Svalbard/16S_diversity_2021/")

#Packages--------
#load packages and functions
library(tidyverse)
library(ggpubr)
library(compositions)
library(zCompositions)
library(phyloseq)
library(ggalt)
library(gridExtra)
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("ggalt")
library(readxl)
library(Biostrings)
library(scales)
library(ggpubr)
library(grid)
library(reshape2)
library(tidyr)
library(readxl)
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library(gridExtra)
library(vegan)
library(ade4)
library(ggpubr)
library(microbiome)
library(ggplot2)
library(dplyr)
library(datasets)
library(ggrepel)
library(dada2)
library(vegan)
# install.packages("pairwise.adonis")
library(pairwiseAdonis)
#install.packages("dada2")
# devtools::install_github("benjjneb/dada2")
library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
theme_set(theme_bw())
`%notin%` <- function(x,y) !(x %in% y)

library(ade4)
library(ape)
library(Biostrings)
library(car)
library(cluster)
library(compositions)
library(dada2)
library(datasets)
library(farver)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(igraph)
library(IRanges)
library(lme4)
library(rstatix)
library(ShortRead)
library(sp)
library(stats)
library(SummarizedExperiment)
library(utils)
library(vegan)

#default is 999 permutations
library(pairwiseAdonis)
library(ShortRead)
library(rstatix)
library(DESeq2)
library(ggfortify)
# install.packages("ggfortify")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")


#####Creating phyloseq object ------
#Load the metadata:
Master_data <- read_excel("~/Documents/Steve Lab/Svalbard/18S_diversity_2021/SunSpears July 2021 master Richness_mod_stripes_labels   .xlsx")
Master_data
#Let's change the colunm name of sample by sampleID to have same column names in both dataframes
names(Master_data)[1] <- "sampleID"
#Let's remove colummns wich are not useful for this analysis
colnames(Master_data)
metadata<- subset(Master_data, select = c("sampleID", "distance...2", "stripe/site"))
names(metadata)[2] <- "distance"
names(metadata)[3] <- "stripe_site"
metadata <- subset(metadata, stripe_site %in% c("Stripe", "Off-stripe"))
metadata <- metadata %>%
  mutate(stripe_site = case_when(
    stripe_site == "Stripe" ~ "DDs",
    stripe_site == "Off-stripe" ~ "non-DDs",
    TRUE ~ stripe_site # Keep other values unchanged
  ))

metadata #One way to identify this is correct with the updated Stripe labels: 
#3 Off-stripes, 3 Stripes, and 4 Off-stripes in sequence

#We are going to work with 4 groups:
# 18S ESVs from 0 to 250 meters on the stripe
# 18S ESVs from 0 to 250 meters OFF the stripe
# 18S ESVs from 300 to 850 meters on the stripe
# 18S ESVs from 300 to 850 meters OFF the stripe
metadata <- subset(metadata, distance <= 850)

metadata$distance <- ifelse(metadata$distance <= 215, "0-215", metadata$distance)
metadata$distance <- ifelse(metadata$distance >= 315 & metadata$distance <= 850, "315-850", metadata$distance)
# Create a new column with the updated values based on conditions
metadata$sample_type <- ifelse(metadata$distance == "0-215" & metadata$stripe_site == "DDs", "0-215 DDs",
                               ifelse(metadata$distance == "0-215" & metadata$stripe_site == "non-DDs", "0-215 non-DDs",
                                      ifelse(metadata$distance == "315-850" & metadata$stripe_site == "DDs", "315-850 DDs",
                                             ifelse(metadata$distance == "315-850" & metadata$stripe_site == "non-DDs", "315-850 non-DDs", NA))))
metadata

otu_table_16 <- read.table("otu_table_16S.txt", sep="\t", header=TRUE, row.names=1)
# print updated column names
print(colnames(otu_table_16)) #168 samples, which comes from the seqtab_final.txt (DADA run)
otu_table_16
# Find the common sample IDs between the two tables
new_otu_table_16<- subset(otu_table_16, select = metadata$sampleID)

# tax_table <- read.table("seqtab_wTax_mctoolsr_preadj.txt", header=TRUE, row.names=1)
tax_table <- read.table("seqtax.txt", sep="\t",header = FALSE, row.names = 1)
#I have changed the results of Burkholderia-Caballeronia-Paraburkholderia into Burkholderia
colnames(tax_table) <-c("Kingdom","Phylum","Class","Order","Family","Genus")

# remove chloroplast, euks, and mitochondria
tax_table <- subset(tax_table,Kingdom != "Eukaryota")
tax_table <- subset(tax_table,Order != "Chloroplast")
tax_table <- subset(tax_table,Family != "Mitochondria")
tax_table <- subset(tax_table,Order != "Staphylococcales")
tax_table <- subset(tax_table,Order != "Enterobacterales")
tax_table <- subset(tax_table,Genus != "Streptococcus")
contam <- read.csv('contam_esvs_16s.csv',header = FALSE)
tax_table <- tax_table[rownames(tax_table) %notin% contam$V1,]

#I have here that the tax table have less numbers of ESV, so I have to remove ESVs from otu table

# new_tax_table_16 <- subset(tax_table, rownames(tax_table) %in% rownames(new_otu_table_16))
new_otu_table_16_2 <- subset(new_otu_table_16, rownames(new_otu_table_16) %in% rownames(tax_table))


# Create a new phyloseq object using the filtered sampledata
rownames(tax_table)
rownames(new_otu_table_16_2)
metadata <- as.data.frame(metadata)
row.names(metadata) <-metadata$sampleID
metadata$sampleID <-NULL

# Sort both data frames in the same order
taxa_sorted <- tax_table[rownames(new_otu_table_16_2), ]
rownames(taxa_sorted)
rownames(new_otu_table_16_2)
class(taxa_sorted)
class(new_otu_table_16_2)
class(metadata)
taxa_matrix = as.matrix(taxa_sorted)
OTU = otu_table(new_otu_table_16_2,taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(taxa_matrix)
DATA = sample_data(metadata)
physeq <- phyloseq(OTU,TAX,DATA)
#readyyyyyyyyyyy!



physeq_16_filtered <- prune_taxa(taxa_sums(physeq) > 0, physeq)
physeq <- physeq_16_filtered
physeq@sam_data
#Richness-----
a_my_comparisons <- list( c("0-215 DDs", "0-215 non-DDs"), c("315-850 DDs", "315-850 non-DDs"),
                          c("0-215 DDs", "315-850 DDs"), c("0-215 non-DDs", "315-850 non-DDs"))
colors <- c(
  "0-215 DDs" = "#ff686b",  
  "0-215 non-DDs" = "#4cc9f0",  
  "315-850 DDs" =  "#b30000",
  "315-850 non-DDs" = "#005f99")

symnum.args <- list(
  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("***", "**", "*", ".", "ns")
)

Richness_16S <- plot_richness(physeq, x = "sample_type", measures = c("Observed"), color = "sample_type") +
  geom_boxplot() +
  geom_jitter(size = 3.5, alpha = 0.6, width = 0.1)+
  scale_color_manual(
    name = "Sample groups",   # 👈 aquí le das el título que quieras
    values = colors
  )+   theme(legend.position = "right",
             axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), 
             plot.title = element_text(face = "bold"),
             plot.margin = margin(t = 5, r = 5, b = 5, l = 5)) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = a_my_comparisons, 
                     label = "p.signif", 
                     symnum.args = symnum.args,
                     vjust = 0.2,
                     size = 6.5) +
  labs(x = NULL, y = "Alpha diversity measure")+
  scale_x_discrete(expand = expansion(add = 0.5))+
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15))+
  theme(
    strip.text = element_text(size = 13, face = "bold"),   # 👈 tamaño del "Observed"
    strip.background = element_rect(fill = "white", color = "black")
  )+
  coord_cartesian(ylim = c(0, 870))
# +ggtitle("A) Alpha diversity in Prokaryotes")
Richness_16S

ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/Richness_16S.png", 
       Richness_16S, width = 4.5 , height = 5)

#####Beta diversity -----
pseq.rel <- microbiome::transform(physeq, "compositional")
otu <- abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

#run PERMANOVA of total microbial composition comparisons.
permanova <- adonis2(t(otu) ~ sample_type,
                     data = meta, permutations=99, method = "bray")
permanova #0.01 there are significantly differences between samples based on distance+stripe_site

# Perform pairwise PERMANOVA on Bray-Curtis dissimilarity matrix
DistBC = phyloseq::distance(pseq.rel, method = "bray")
metadata <- sample_data(pseq.rel)

set.seed(123) 
pairwise_adonis <- pairwise.adonis(
  DistBC, 
  metadata$sample_type, 
  perm = 1000)
df_pairwise_adonis <- as.data.frame(pairwise_adonis)
df_pairwise_adonis

# Create the table grob
table <- tableGrob(df_pairwise_adonis)
# Create the title grob with a blank grob on the right side
title <- arrangeGrob(
  textGrob("16S rRNA - Pairwise PERMANOVA", gp = gpar(fontsize = 16, fontface = "bold")),
  rectGrob(gp = gpar(fill = NA, col = NA)),  # Blank grob
  widths = c(0.9, 1),
  ncol = 2)
# Combine the title and table grobs
pairwise_permanova_16 <- grid.arrange(title, table, ncol = 1, heights = c(0.1, 1), top = 0.1, left = 0.1)

# Save the plot as a PNG file
ggsave("Pairwise_PERMANOVA_16.png", pairwise_permanova_16, width = 9.9, height = 2.2, dpi = 300)

####Dispersion------
#Let's also corroborate the distance as a factor of differeces with betadisper (Permdisper):
metadata_betad <- data.frame(sample_data(pseq.rel))
mod_betad <- betadisper(DistBC, metadata_betad$sample_type) #Betadisper for the 4 groups
anova_res <- anova(mod_betad)
print(anova_res)
perm_test <- permutest(mod_betad, pairwise = TRUE, permutations = 999)
print(perm_test) #permutation test

#Let's do it per group: 
#Pvalue <0.05 means dispersion is significantly different in the samples of that group.

sub_samples <- which(metadata_betad$distance == "0-215")
dist_sub <- as.dist(as.matrix(DistBC)[sub_samples, sub_samples])
group_sub <- metadata_betad$sample_type[sub_samples]
mod_sub <- betadisper(dist_sub, group_sub)
anova(mod_sub) #F=3.18, Pvalue=0.085 #Homogeneous dispersion

sub_samples <- which(metadata_betad$distance == "315-850")
dist_sub <- as.dist(as.matrix(DistBC)[sub_samples, sub_samples])
group_sub <- metadata_betad$sample_type[sub_samples]
mod_sub <- betadisper(dist_sub, group_sub)
anova(mod_sub) #F= 1.8242, Pvalue=0.1905 #Homogeneous dispersion

sub_samples <- which(metadata_betad$stripe_site == "DDs")
dist_sub <- as.dist(as.matrix(DistBC)[sub_samples, sub_samples])
group_sub <- metadata_betad$sample_type[sub_samples]
mod_sub <- betadisper(dist_sub, group_sub)
anova(mod_sub) #F=148.3, Pvalue= 0.000000000005235 #Different variance in samples inside of this group.

sub_samples <- which(metadata_betad$stripe_site == "non-DDs")
dist_sub <- as.dist(as.matrix(DistBC)[sub_samples, sub_samples])
group_sub <- metadata_betad$sample_type[sub_samples]
mod_sub <- betadisper(dist_sub, group_sub)
anova(mod_sub) #F= 8.7404 and Pvalue=0.006541 #Different variance in samples inside of this group.

#Interpretation: 
#In 0-215m samples: This validates that any differences found in PERMANOVA are due to location (centroid shifts) and not just differences in group variance.
#In 315-850m samples: This validates marginally that any differences found in PERMANOVA are due to location (centroid shifts) and not just differences in group variance.
#In DDs transect: Variation is not equal among there distances subgroups.###Here to be careful on the PERMANOVA interpretation.
#In non-DDs transect: Variation is not equal among there distances subgroups.###Here to be careful on the PERMANOVA interpretation.


