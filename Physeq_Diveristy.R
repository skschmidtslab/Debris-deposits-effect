#Processing for ASV in Euk
#Set working directory
# setwd("~/Artics_16S_18S")
#in my local computer
setwd("~/Documents/Steve Lab/Svalbard/18S_diversity_2021/")

#######Load packages ------##########
library(ggalt)
library(ggplot2)
library(gridExtra)
library(readxl)
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
library(ggpubr)
library(microbiome)
library(ggplot2)
library(dplyr)
# install.packages("pairwise.adonis")
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
library(ggplot2)
library(ggrepel)
library(igraph)
library(IRanges)
library(lme4)
library(stats)
library(SummarizedExperiment)
library(utils)
library(vegan)
library(pairwiseAdonis)
library(ShortRead)
library(rstatix)
library(DESeq2)
library(phyloseq)
library(sp)
library(GenomicRanges)
library(farver)

#######Import data ------###############
#Adam created a folder where he put the otutable of the 18S ESV DADA run. 
#The folder is svalbard_18S and this had a sub folder called 03_tabletax and 18s_barcodes.txt
#Dr.Steve sent me a file with the metadata (depth, distance, stripe or not stripe, etc) 

Master_data <- read_excel("SunSpears July 2021 master Richness_mod_stripes_labels   .xlsx")
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

metadata
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

########PHYLOSEQ OBJECT---------------#######
otu_table <- read.table("~/Documents/Steve Lab/Svalbard/svalbard_18S/03_tabletax/seqtab_final.txt", sep="\t", header=TRUE, row.names=1)
#Let's change the samples names removing part of the name "Undetermined_S0_L001_R1_001.fastq_"
colnames(otu_table) <- gsub("Undetermined_S0_L001_R1_001.fastq_", "", colnames(otu_table))
# print updated column names
print(colnames(otu_table)) #168 samples, which comes from the seqtab_final.txt (DADA run)
# Find the commun sample IDs between the two tables
new_otu_table<- subset(otu_table, select = metadata$sampleID)


tax_table <- read.table("seqtab_wTax_mctoolsr_preadj.txt", header=TRUE, row.names=1)
colnames(tax_table) <-c("Kingdom","Phylum","Class","Order","Family","Genus")
# the file comes from this file, but I just eliminate manually the Samples ID and just 
# leave the taxonomy names and the OTU numbers
# tax_table <- read.table("~/Documents/Svalbard/svalbard_18S/03_tabletax/seqtab_wTax_mctoolsr.txt", header=TRUE, row.names=1)
# colnames(tax_table) <- gsub("Undetermined_S0_L001_R1_001.fastq_", "", colnames(tax_table))
#Now les clean the OTU or the TAX table from contamination (thanks to Pacifica)
#4,267 ESVs
# remove chloroplast, euks, and mitochondria
library(dplyr)
# remove chloroplast, euks, and mitochondria
dim(tax_table)
sum(tax_table[, "Family"] == "Mitochondria", na.rm = TRUE)
tax_table1 <- tax_table[!(tax_table[, "Family"] %in% "Mitochondria"), ]
dim(tax_table1)
tax_table2 <- tax_table1[!(tax_table1[, "Kingdom"] %in% "Bacteria"), ]
dim(tax_table2)
tax_table3 <- tax_table2[!(tax_table2[, "Order"] %in% "Chloroplast"), ]
dim(tax_table3)
Phylum_to_remove <- c("Arthropoda", "Nematoda", "Vertebrata",
                      "Annelida")
tax_table4 <- tax_table3[!(tax_table3[, "Phylum"] %in% Phylum_to_remove), ]
dim(tax_table4)
Order_to_remove <- c("Liliopsida", "Caryophyllales", "Asterales",
                     "Brassicales", "Solanales", "Malvales", "Fabales",
                     "Sapindales", "Apiales", "Pinales", "Rosales")
tax_table5 <- tax_table4[!(tax_table4[, "Order"] %in% Order_to_remove), ]
dim(tax_table5)
#Other ASVs assigned in blast as: Solanum, and other flowering plants: ESV_14, ESV_338.
ASV_to_remove <- c("ESV_14", "ESV_338", "ESV_1011", "ESV_1029") #ASV_1011 is an arthropod, the others are flowering plants according to blast
tax_table6 <- tax_table5[!(rownames(tax_table5) %in% ASV_to_remove), ]
dim(tax_table6) #Done
#I also have a file for some other ASVS (the majority should be the same as the ones assigned to flowering plants..).
selected_esvs <- read.csv("~/Documents/Steve Lab/Svalbard/18S_diversity_2021/Selected_esvs_18s_v2.csv")
tax_table7 <- tax_table6[rownames(tax_table6) %notin% selected_esvs$ESV_ID,]  
dim(tax_table7) #4030 ASVs
#Merging otu table to not have the removed ASVs
otu_table_18S <- subset(new_otu_table, rownames(new_otu_table) %in% rownames(tax_table6)) 
dim(otu_table_18S)

#Creating a new phyloseq object using the filtered sampledata
rownames(tax_table6)
rownames(otu_table_18S)
metadata <- as.data.frame(metadata)
row.names(metadata) <-metadata$sampleID
metadata$sampleID <-NULL

#Sort both data frames in the same order
taxa_sorted <- tax_table6[rownames(otu_table_18S), ]
rownames(taxa_sorted)
rownames(otu_table_18S)
class(taxa_sorted)
class(otu_table_18S)
class(metadata)
taxa_matrix = as.matrix(taxa_sorted)
OTU = otu_table(otu_table_18S,taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(taxa_matrix)
DATA = sample_data(metadata)
physeq_18 <- phyloseq(OTU,TAX,DATA)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4030 taxa and 55 samples ]
# sample_data() Sample Data:       [ 55 samples by 3 sample variables ]
# tax_table()   Taxonomy Table:    [ 4030 taxa by 6 taxonomic ranks ]

#####Filtering physeq taxa sums ------
physeq_18@sam_data
physeq_18_filtered <- prune_taxa(taxa_sums(physeq_18) > 0, physeq_18)
see_asv<-as.data.frame(physeq_18_filtered@otu_table) #to check the asvs
physeq_18_F <- prune_samples(sample_sums(physeq_18_filtered) > 0, physeq_18_filtered)
physeq_18S <- physeq_18_F
see_asv2<-as.data.frame(physeq_18S@otu_table) #to check the asvs
sample_data(physeq_18S)


######physeq atributes------
#Let's know the features of our physloeq object
sample_variables(physeq_18S) #output: [1] "distance"    "stripe_site"
ntaxa(physeq_18S) #output: [1] 2044
nsamples(physeq_18S) #output: [1] 54
rank_names(physeq_18S) #output: [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus" 
#Everything looks good!but we don't have a phylogenetic tree.

######Alpha diveristy ------
a_my_comparisons <- list( c("0-215 DDs", "0-215 non-DDs"), c("315-850 DDs", "315-850 non-DDs"),
                          c("0-215 DDs", "315-850 DDs"), c("0-215 non-DDs", "315-850 non-DDs"))
symnum.args = list(cutpoints = c(0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
colors <- c(
  "0-215 DDs" = "#ff686b",      # rosado claro
  "0-215 non-DDs" = "#4cc9f0",   # celeste claro
  "315-850 DDs" =  "#b30000",
  "315-850 non-DDs" = "#005f99"  # azul intenso
)

symnum.args <- list(
  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("***", "**", "*", ".", "ns")
)

physeq_18S@sam_data$sample_type

###Plot - Richness:
positions_y <- c(247, 280, 307, 339)
Richness_18S <- plot_richness(physeq_18S, x = "sample_type", measures = c("Observed"), color = "sample_type") +
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
                     label.y = positions_y,
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
  coord_cartesian(ylim = c(0, 360))
Richness_18S
#Save plot
ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/Richness_18S_v2.png", 
       Richness_18S, width = 4.5 , height = 5)
