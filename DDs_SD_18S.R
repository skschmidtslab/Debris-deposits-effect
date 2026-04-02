#Heatmap for a comparison of debris deposits (DDs) and Supraglacial debris (SD) Top ASVs:
#This comparison is made to test if ASVs from 
#debris deposits could came from the supraglacial debris on top of the glacier.
#First I need to build a phyloseq object integrating both ASV seqtabs:
#DDs sequencing data was run in 2021.
#Supraglacial debris sequencing data was run in 2023.

#Setting the working directory
setwd("~/Documents/Steve Lab/Svalbard/DDs_SD/18S/")
#Load packages:
library(dada2)
library(readxl)
library(dplyr)
library(phyloseq)
theme_set(theme_bw())
`%notin%` <- function(x,y) !(x %in% y)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tidyverse)
library(tibble)
library(pheatmap)
library(ggplot2)
library(grid)

#List seqtabs with their paths
seqtab1 <- readRDS("~/Documents/Steve Lab/Svalbard/svalbard_18S/03_tabletax/seqtab_final.rds") #####path of the SD seqtab
seqtab2 <- readRDS("~/Documents/Steve Lab/Svalbard 2023/18S/seqtab_final (2).rds") #####path of the DDs seqtab
seqtab <- mergeSequenceTables(seqtab1, seqtab2) #Join seqtabs
seqtab_combined.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE) #Remove Chimeras
#Save it without chimeras:
saveRDS(seqtab_combined.nochim,"~/Documents/Steve Lab/Svalbard/DDs_SD/18S/seqtab_DDs_SD_combined.nochim.rds")

#To make a single taxonomy assigment for both:
###Taxonomy assignation: -----Run ONLY ONE TIME ----Then comment it and call it from file
# tax <- assignTaxonomy(seqtab_combined.nochim,
#                        "~/Documents/Steve Lab/Svalbard 2023/18S/pr2_version_5.0.0_SSU_dada2.fasta.gz",
#                       tryRC = TRUE,
#                       multithread = TRUE)
# table.fp <- "~/Documents/Steve Lab/Svalbard/DDs_SD/18S" ###Folder where I will locate my tax rds file
# head(tax)
# saveRDS(tax, paste0(table.fp, "/tax_combined_final.rds"))

#Now let's call the elements of the phyloseq object
#######Call Metadata------
Master_data_combined <- read_excel("metadata_svalbard2023_2021_combined.xlsx")
colnames(Master_data_combined)
names(Master_data_combined)[1] <- "sampleID"
#Subset (Only selecting the SD and DDs sample types) #Black_1 are samples from supraglacial debris, 
#Black_3 is the extension collected the same yaer and DDs is the debris deposits collected on the forefield.
Master_data_DDs_SD <- subset(Master_data_combined, Sample_set %in% c("Black_1", "Black_3", "DDs"))

#######Call Otutable------
seqtab.nochim <- readRDS("seqtab_DDs_SD_combined.nochim.rds")
seqtab.t <- as.data.frame(t(seqtab.nochim)) 
colnames(seqtab.t) ###This seqtab has all the samples from those runs (column names), we need to filter them.
#first, lets remove the prefixe "undetermined..." from all sample names:
cn <- colnames(seqtab.t) #Saving original names
cn[1:168] <- sub(".*fastq_", "", cn[1:168]) #Modifiy them until 168 sample.
colnames(seqtab.t) <- cn #Replace the object
colnames(seqtab.t) #See names now...it is cleaner.
#Here, let's select only the samples from DDs and SD:
Master_data_DDs_SD$sampleID #We want these ones.

seqtab.t<- subset(seqtab.t, select = Master_data_DDs_SD$sampleID) #Selecting the samples of this 1st sample set
dim(seqtab.t) #[1] 7159   25
#Remove ESVs which are 0 in all samples
rows_to_remove <- rowSums(seqtab.t[,]) == 0 #This will remove ESVs which are zero in at least one of your samples
#It is very important to remove those ESV as they are not part of any proportion of your samples.
seqtab.t <- seqtab.t[!rows_to_remove, ] #Ok, removed!
#(Here is time to look the ESVs on Blast and remove the taxas you consider are not having sense with your research
#Like potatos or E.coli in the glacier forefields).
dim(seqtab.t) #[1] 719  25

#######Pull out ESV repset-----
#This is giving numbers to each specific ESV
rep_set_ASVs <- as.data.frame(rownames(seqtab.t))
rep_set_ASVs <- mutate(rep_set_ASVs, ASV_ID = 1:n())
rep_set_ASVs <- mutate(rep_set_ASVs, ASV_ID = 1:nrow(rep_set_ASVs))
rep_set_ASVs$ASV_ID <- sub("^", "ASV_", rep_set_ASVs$ASV_ID)
rep_set_ASVs$ASV <- rep_set_ASVs$`rownames(seqtab.t)` 
rep_set_ASVs$`rownames(seqtab.t)` <- NULL
#Let's replace the ESV sequences by the ESV numbers in the table:
rownames(seqtab.t) <- rep_set_ASVs$ASV_ID 
seqtab.t$ASV_ID <- rownames(seqtab.t)
rownames(seqtab.t) <- NULL
seqtab.t <- seqtab.t[, c("ASV_ID", setdiff(colnames(seqtab.t), "ASV_ID"))]
rownames(seqtab.t) <- NULL
seqtab.t#Looks good. #I have 719 ESVs which are present in at least one of the 26 samples, blank or controls.
#Let's save it as our OTU table:
write.table(seqtab.t, "otutable_combined_2021_2023_DDs_SD", sep="\t", quote=TRUE, row.names=FALSE, col.names = TRUE)
write.csv(rep_set_ASVs, file = "rep_set_ASVs.csv") #Save que sequences with their ASV_ID, just in case to look them into Blast.

#######Tax Table------
tax <- readRDS("tax_combined_final.rds") #Reading the tax data
taxonomy <- as.data.frame(tax) #Converting into a dataframe
head(taxonomy) #As you can see this still have the ASVs that you need to remove.
#Lets do the same with the taxonomy table, ASV numbers
taxonomy$ASV <- as.factor(rownames(taxonomy)) #This is just putting the ASVs sequences into a new column
taxonomy <- merge(rep_set_ASVs, taxonomy, by = "ASV") #This is merging the two tables content

rownames(taxonomy) <- taxonomy$ASV_ID
rownames(seqtab.t) <- NULL
rownames(taxonomy) <- NULL
taxonomy <- subset(taxonomy, select = -c(ASV))
colnames(taxonomy)
colnames(taxonomy) <-c("ASV_ID","Domain","Supergroup","Kingdom", "Phylum","Class","Order", "Family", "Genus", "Species")
write.table(taxonomy, "taxonomy_combined_final.txt", sep="\t", quote=TRUE, row.names=FALSE, col.names = TRUE)

#######Transforming tables for phyloseq object----
otu_table_SDD_SD <- read.table("otutable_combined_2021_2023_DDs_SD", sep="\t", header=TRUE, row.names=1)
print(colnames(otu_table_SDD_SD)) #6 samples (3 exp + 1 blank and 2 controls)
otu_table_SDD_SD
tax_table <- read.table("taxonomy_combined_final.txt", sep="\t",header = TRUE, row.names = 1)
colnames(tax_table) <-c("Domain","Supergroup","Kingdom", "Phylum","Class","Order", "Family", "Genus", "Species")

#Review the classs names...it looks like the class is now the family, and the family a phyluim etc...
####################################

dim(tax_table)
# sum(tax_table[, "Family"] == "Mitochondria", na.rm = TRUE)
sum(grepl(":mito$", tax_table[, "Domain"]), na.rm = TRUE) #1 ...#But I noticed that sequences
#assigned to mitochondria in one of the seqtabs (DDs) are not assigned as that anymore...
#they actually don't have assignation...it looks like I will need to verify them in blast. 
#That happens when databases get updated but that is ok! I can also make a match of the sequences bymyself...after getting all the dataframes organized.
# Filtramos manteniendo solo las filas que NO terminan en ":mito"
tax_table1 <- tax_table[!grepl(":mito$", tax_table[, "Domain"]), ]
dim(tax_table1)
tax_table2 <- tax_table1[!(tax_table1[, "Domain"] %in% "Bacteria"), ]
dim(tax_table2)
tax_table3 <- tax_table2[!grepl(":chloroplast$", tax_table2[, "Domain"]), ]
dim(tax_table3)
Class_to_remove <- c("Arthropoda", "Nematoda", "Vertebrata",
                      "Annelida")
tax_table4 <- tax_table3[!(tax_table3[, "Class"] %in% Class_to_remove), ]
dim(tax_table4)
Order_to_remove <- c("Liliopsida", "Caryophyllales", "Asterales",
                     "Brassicales", "Solanales", "Malvales", "Fabales",
                     "Sapindales", "Apiales", "Pinales", "Rosales")
tax_table5 <- tax_table4[!(tax_table4[, "Order"] %in% Order_to_remove), ]
dim(tax_table5) #So, now I have to find...what to filter...Depending on the sequences that match the ones I filtered in the 18S processing of DDs.
########Generating a file for filtering out ------
#ASVs from DDs vs non-DDs:
Lines <- readLines("~/Documents/Steve Lab/Svalbard/svalbard_18S/03_tabletax/repset.fasta")
h_idx <- grep("^>", Lines)
fasta_df <- data.frame(
  header   = Lines[h_idx],
  sequence = Lines[h_idx + 1],
  stringsAsFactors = FALSE)
csv_df <- read.csv("~/Documents/Steve Lab/Svalbard/DDs_SD/18S/rep_set_ASVs.csv")
coincidences_IDs <- inner_join(fasta_df, csv_df, by = c("sequence" = "ASV"))
colnames(coincidences_IDs) <- c("ASV_ID_DDs_lonng", "sequence", "X", "ASV_ID_SD_DDs")
coincidences_IDs <- coincidences_IDs[, c("ASV_ID_SD_DDs", "sequence", "ASV_ID_DDs_lonng")]
message("Total matches: ", nrow(coincidences_IDs))
# I have 473 sequences which match with both seqtabs...
# From there...I should select the ones I considered in the DDs vs non-DDs as not necessary (to have same filtering process).
# Now I call the ASVs I filtered out from the DDs vs non-DDs comparison: 
selected18Sr <- read.csv("~/Documents/Steve Lab/Svalbard/18S_diversity_2021/Selected_esvs_18s_v2.csv")
coincidences_IDs_clean <- coincidences_IDs %>%
  mutate(ASV_ID_DDs_clean = gsub("^>", "", ASV_ID_DDs_lonng)) %>%
  separate(ASV_ID_DDs_clean, into = "ASV_ID_DDs", sep = " ", extra = "drop")
colnames(selected18Sr)
#Let's find them here: 
final_remov <- coincidences_IDs_clean %>%
  filter(ASV_ID_DDs %in% selected18Sr$ESV_ID) %>%
  dplyr::select(ASV_ID_SD_DDs,ASV_ID_DDs, sequence)
write.csv(final_remov, "ASV_to_remove_DDs_SD.csv", row.names = FALSE) #Saving it as a file
#####Continuing with the filtering out-----
selected_esvs <- read.csv('ASV_to_remove_DDs_SD.csv')
tax_table6 <- tax_table5[rownames(tax_table5) %notin% selected_esvs$ASV_ID_SD_DDs,]  
dim(tax_table6) #[1] 694   9.  #19 less than in taxtable5
#Let's see what I still have:
unique(tax_table6$Class)
#From there, the one could be a worthy to look inside is Embryophyceae...
unique(tax_table6$Genus) #Here everything is more clear
Genus_to_remove <- c("Pinus", "Morus", "Persea",
                     "Cedrus", "Helianthus", "Cucumis",
                     "Setaria", "Artemisia")
tax_table7<- tax_table6[!(tax_table6[, "Genus"] %in% Genus_to_remove), ]
dim(tax_table7) #[1] 666   9

#Now is time to check the taxa names of the remaining ASVs. Some of them could be missassigned. Blast will help.
#But before, let's see which ones are good for doing it. Let's see which ones are the top of the ASVs.

#Final Otu table:
new_otu_table <- subset(otu_table_SDD_SD, rownames(otu_table_SDD_SD) %in% rownames(tax_table7))
dim(new_otu_table) #[1] 666 25
dim(otu_table_SDD_SD)

#######Phyloseq object-----
rownames(new_otu_table)
metadata <- as.data.frame(Master_data_combined)
row.names(metadata) <-metadata$sampleID
metadata$sampleID <-NULL
taxa_table_sorted <- tax_table7[rownames(new_otu_table), ] #para ordenar
rownames(taxa_table_sorted)
rownames(new_otu_table)
class(taxa_table_sorted)
class(new_otu_table)
class(metadata)
taxa_matrix = as.matrix(taxa_table_sorted)
OTU = otu_table(new_otu_table,taxa_are_rows = TRUE)
TAX = tax_table(taxa_matrix)
DATA = sample_data(metadata)
physeq_SD_DDs<- phyloseq(OTU,TAX,DATA)
#readyyyyyyyyyyy!

#########Heatmap for matching ASVs and comparing abundances -----
#Let's get the abundances from the phyloseq object
esv_abundance <- physeq_SD_DDs %>%
  psmelt()
head(esv_abundance)
dim(esv_abundance)
colnames(esv_abundance)
# Transpose, join with metadata, then summarize
esv_wide <- esv_abundance %>%
  filter(Sample_set %in% c("DDs", "Black_1")) %>% #This only takes these two sample types from the total samples of SD_DDs
  dplyr::select(OTU, Sample, Abundance) %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)
#Create abundances matrix
mat <- esv_wide %>%
  group_by(OTU) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("OTU") %>%
  filter(rowSums(.) > 0) #When removing some samples, some ASVs will have 0 in the remining samples...so we need to filter them out.
col_annotation <- esv_abundance %>%
  filter(Sample_set %in% c("DDs", "Black_1")) %>%
  distinct(Sample, Sample_set) %>%
  column_to_rownames("Sample")
top50_OTUs <- mat %>%
  rowSums() %>%
  sort(decreasing = TRUE) %>%
  head(50) %>%
  names()

######Heat map---------
#Creating the table of the Top50:
top50 <- mat[top50_OTUs, ]
desired_order <- c(
  "BLK1","BLK2","BLK3","BLK4","BLK5","BLK6",
  "SUN0021","SUN0042","SUN0043","SUN0061","SUN0062",
  "SUN0081","SUN0082","SUN0083","SUN0101","SUN0102","SUN0103")
existing_cols <- desired_order[desired_order %in% colnames(top50)]
top50_ordered <- top50[, existing_cols]

#Now I will save my top50 as a file:
top50_sequences_final <- rep_set_ASVs %>%
  filter(ASV_ID %in% top50_OTUs) %>%
  dplyr::select(ASV_ID, ASV) 
dim(top50_sequences_final) 
write.csv(top50_sequences_final, "Top50_Heatmap_Sequences.csv", row.names = FALSE)
#Now...on that csv file I have to add new columns for the blast annotation and for the DDs analysis taxa assignation
#Now, I have the new edited file (by hand...I had to seach on Blast and add their taxa matches <97% of Per. Ident)
top50_sequences_final_TaxByBlast <- read.csv("~/Documents/Steve Lab/Svalbard/DDs_SD/18S/Top50_Heatmap_Sequences_blast.csv")
top50_sequences_final_TaxByBlast
colnames(top50_sequences_final_TaxByBlast)
#Now I have to add the column "ASV_final_tax" from this in the top50_ordered
# Unir la columna de taxonomía de Blast a tu matriz (convertida temporalmente a df)

asv_ids <- sub(":.*", "", rownames(top50_ordered))
blast_taxa <- top50_sequences_final_TaxByBlast$ASV_final_tax[match(asv_ids, top50_sequences_final_TaxByBlast$ASV_ID)]
new_labels <- paste0(asv_ids, ": ", blast_taxa)
rownames(top50_ordered) <- new_labels

#Heatmap
p18S <- pheatmap(top50_ordered,
                 cluster_cols = FALSE,     #Keeping the order
                 cluster_rows = TRUE,      #Grouping ASVs
                 annotation_col = col_annotation,
                 annotation_colors = list(
                   Sample_set = c("SD" = "#BEBFC5", "DDs" = "#ff686b")
                 ),
                 show_rownames = TRUE,
                 color = colorRampPalette(c("white", "blue", "black"))(100),
                 border_color = NA,
                 main = "A) Eukaryotes 18S: SD vs DDs")
ggsave("Black_1_SDD_comparison_HM.png", p18S, width = 6.8, height = 7, dpi = 500) #Save plot.

