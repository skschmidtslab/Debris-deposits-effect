#Get the relative abundances of phototrophs from the Euk community
#Load libraries
library(tidyverse)
library(dplyr)
library(phyloseq)
#Set directory
setwd("~/Documents/Steve Lab/Svalbard/18S_diversity_2021/phototrophs")
#Look for palletes here:
###https://colorhunt.co/palette/bfccb57c96abb7b7b7edc6b1

#------Breaking out the phyloseq in the groups samples------
#############0-215 DDs------
physeq_215DDs <- subset_samples(physeq_18, sample_type == "0-215 DDs")
#Remove ESVs which have 0 counts in all samples belonging to 0-215 DDs group/
physeq_215DDsf <- prune_taxa(taxa_sums(physeq_215DDs) > 0, physeq_215DDs)
physeq_215DDsf
otu_table(physeq_215DDsf)
sample_data(physeq_215DDsf)
tax_table(physeq_215DDsf)
#Extraction of the taxa as data frame from the phyloseq object.
taxa_df <- as.data.frame(tax_table(physeq_215DDsf))
#Converting the ESVs into ASVs.
taxa_df <- taxa_df %>%
  rownames_to_column(var = "ASV_ID") %>%      # 1. Crea la columna con los IDs
  mutate(ASV_ID = gsub("ESV", "ASV", ASV_ID)) # 2. Cambia ESV por ASV
#These are the ASVs that needs to be updated as they got selected as the most abundants.
Blast_ASVtaxa <- read.csv("~/Documents/Steve Lab/Svalbard/18S_diversity_2021/Blast_euk_ASV_correction.csv", stringsAsFactors = FALSE)

#Search ASVs and replace them
for (i in 1:nrow(Blast_ASVtaxa)) {
  asv_target <- Blast_ASVtaxa$ASV_ID[i]
  if (asv_target %in% taxa_df$ASV_ID) {
    fila_en_taxa <- which(taxa_df$ASV_ID == asv_target)
    taxa_df[fila_en_taxa, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")] <- 
      Blast_ASVtaxa[i, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
    message("success: ", asv_target)
  } else {
    warning("They are not in taxa_df: ", asv_target)
  }
}
#ASV_48 and ASV_59 are not in taxa_df, it must be bc they are not in these group of samples. That is ok!

rownames(taxa_df) <- taxa_df$ASV_ID
# taxa_df$ASV_ID <- NULL 
taxa_matrix <- as.matrix(taxa_df)
taxa_names(physeq_215DDsf) <- gsub("ESV", "ASV", taxa_names(physeq_215DDsf))
tax_table(physeq_215DDsf) <- tax_table(taxa_matrix)
physeq_215DDsf #Verifying
head(tax_table(physeq_215DDsf))

####Let's clean 
physeq_215DDsf_filtered <- prune_taxa(taxa_sums(physeq_215DDsf) > 0, physeq_215DDsf)
physeq_215DDsf_filtered2 <- prune_samples(sample_sums(physeq_215DDsf_filtered) > 0, physeq_215DDsf_filtered)
####Let's convert into relative abundances.
ps.rel_18 = transform_sample_counts(physeq_215DDsf_filtered2, function(x) x/sum(x)*100)
ps.rel_18_filtered <- prune_taxa(
  !(tax_table(ps.rel_18)[, "Genus"] %in% c("Spumella", "Pedospumella")),
  ps.rel_18)
ps.rel_18_filtered@tax_table #Spumella is a microalgae (golden algae but it is heterotroph)
####Agglomerate taxa
glom <- tax_glom(ps.rel_18_filtered, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
head(ps.melt) #I am doing this to keep the phylum and visually identify the Mosses and microalgae.
#List of phylums to keep
phylums_to_keep <- c("Cryptomonadales", "Florideophycidae", "Dinoflagellata",
                     "Streptophyta","Phragmoplastophyta", "Ochrophyta", "Chlorophyta_ph", 
                     "Chlorophyta")
ps.melt$Phylum <- as.character(ps.melt$Phylum)
ps.melt$Phylum <- ifelse(ps.melt$Phylum %in% phylums_to_keep, ps.melt$Phylum, "Other phylum")
head(ps.melt)
colnames(ps.melt)
ps.melt$Type <- ifelse(ps.melt$Phylum == "Other phylum", "Non-Phototrophs", "Phototrophs")

###Let's check % of relative abundances:
total_abundance <- sum(ps.melt$Abundance)
phototrophs_abundance <- sum(ps.melt$Abundance[ps.melt$Type == "Phototrophs"])
phototrophs_abundance
other_phylum_abundance <- total_abundance - phototrophs_abundance
other_phylum_abundance

phototroph_percentage <- (phototrophs_abundance / total_abundance) * 100
phototroph_percentage
other_phylum_percentage <- (other_phylum_abundance / total_abundance) * 100
other_phylum_percentage

#Mean and SD values:
phototroph_data <- ps.melt %>%
  filter(Type == "Phototrophs")
mean_abundance_phototroph <- mean(phototroph_data$Abundance, na.rm = TRUE) ##Mean
mean_abundance_phototroph
sd_abundance_phototroph <- sd(phototroph_data$Abundance, na.rm = TRUE) ##SD
sd_abundance_phototroph

########Data frame ------
data_215DDs <- data.frame(
  Phylum = c("Phototrophs", "Non-Phototrophs"),
  Abundance = c(56, 44)  #% of rel abundances
)
data_215DDs$angle <- cumsum(data_215DDs$Abundance) - 0.5 * data_215DDs$Abundance
colors <- c(
  "Phototrophs" = "#F1F5A8",   # azul
  "Non-Phototrophs" = "#FFD0EC"
)

########Pie plot -----------
pie_215DDs <- ggplot(data_215DDs, aes(x = 2, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  geom_text(aes(x = 2, y = Abundance, label = paste0(Abundance, "%")), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y", start = +pi/-12) +  # Girar el donut plot a 90 grados en sentido contrario a las agujas del reloj
  theme_void() +
  xlim(1, 2.5) +
  ggtitle("DDs (0-215m)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = colors) +
  theme(legend.position = "bottom")  # Ajustar la posición de la leyenda si es necesario
pie_215DDs 
ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/215DDs_Euk_pieplot.png", pie_215DDs)

########Genus at the Phylum level and at the last taxon - Phototrophs ------
ps.rel_18_filtered ###Phyloseq of 0-215 
df_physeq <- otu_table(ps.rel_18_filtered)
colnames(df_physeq)
otu_frame <- as.data.frame(df_physeq)
#Add taxonomies as new columns
taxa_df <- as.data.frame(tax_table(ps.rel_18_filtered))
colnames(taxa_df)
df_physeq <- cbind(df_physeq, taxa_df) #Join dataframes 
sample_data_df <- as.data.frame(sample_data(ps.rel_18_filtered))
phylums_to_keep <- c("Cryptomonadales", "Florideophycidae", "Dinoflagellata",
                     "Streptophyta","Phragmoplastophyta", "Ochrophyta", "Chlorophyta_ph",
                     "Chlorophyta")
df_physeq_selected <- subset(df_physeq, Phylum %in% phylums_to_keep)
dim(df_physeq)
# [1] 76 18
dim(df_physeq_selected) #Data with only the phylum we selected
# [1] 29 18 
#Abundance:
abundance<- df_physeq_selected[, grepl("^SUN", colnames(df_physeq_selected))]
#Total for each sample:
sum_sample <- rowSums(abundance)
#Relative abundance:
head(df_physeq_selected) #Let's check how it looks like
colnames(df_physeq_selected) #Good
df_physeq_selected$Relative_Abundance <- rowSums(abundance) / sum(sum_sample) * 100 #We added the rel ab column
#Replace NA by "NA"
df_physeq_selected <- df_physeq_selected %>%
  mutate(across(c(Class, Order, Family, Genus), ~replace_na(.x, "NA")))
#Let's do our rank:
df_genus_abundance <- df_physeq_selected %>%
   select(ASV_ID, Relative_Abundance) #This isolated table make us easy to see the highest abundant in the rank
df_genus_abundance$ASV_ID <- ifelse(df_physeq_selected$Relative_Abundance < 0.5, "ASVs < 0.5%", df_physeq_selected$ASV_ID)
df_genus_abundance #Now, ASV <1% will be recognize all the same as a small percentage ASV.
#Make a summary of the ASVs we are going to plot
df_genus_abundance_summed <- df_genus_abundance %>%
  group_by(ASV_ID) %>%
  mutate(Summed_Relative_Abundance = sum(Relative_Abundance)) %>%
  distinct(ASV_ID, .keep_all = TRUE) %>%
  select(ASV_ID, Summed_Relative_Abundance)
df_genus_abundance_summed #This is the list I am going to use.
#Detecting the last taxa level assigned to the ASVs
df_physeq_selected_LastTax <- df_physeq_selected %>%
  mutate(
    Taxon = case_when(
      !is.na(Genus) & Genus != "NA" ~ paste("Genus:", Genus),
      !is.na(Family) & Family != "NA" ~ paste("Family:", Family),
      !is.na(Order) & Order != "NA" ~ paste("Order:", Order),
      !is.na(Class) & Class != "NA" ~ paste("Class:", Class),
      TRUE ~ paste("Phylum:", Phylum)
    )
  )
#Adding this to the list
df_genus_abundance_summed <- df_genus_abundance_summed %>%
  left_join(
    df_physeq_selected_LastTax %>%
      select(ASV_ID, Taxon),
    by = "ASV_ID") %>%
  mutate(Taxon = ifelse(is.na(Taxon), "Others", Taxon))
df_genus_abundance_summed

library(dplyr)

df_genus_abundance_summed <- df_genus_abundance_summed %>%
  mutate(
    Description = case_when(
      Taxon == "Class: Bryopsida" ~ "Bryopsida (moss)",
      Taxon == "Genus: Atrichopsis sp." ~ "Atrichopsis (mosses)",
      Taxon == "Family: Bryaceae" ~ "Bryaceae (mosses)",
      Taxon == "Genus: Blasia sp." ~ "Blasia (liverworts)",
      Taxon == "Genus: Tortula sp. " ~ "Tortula (moss)",
      Taxon == "Genus: Chloroidium sp." ~ "Chloroidium (microalgae)",
      Taxon == "Others" ~ "Others",
      TRUE ~ Taxon ))

df_genus_abundance_summed
########pieplot ------
#Different colors for each ASV
colors <- c(
  "Bryopsida (moss)" = "#E9A89B",   # azul
  "Atrichopsis (mosses)" = "#B7C9F2",       # naranja
  "Bryaceae (mosses)" = "#EF9C66",       # verde
  "Blasia (liverworts)" = "#78ABA8",
  "Tortula (moss)" = "#FCDC94",
  "Chloroidium (microalgae)" = "#C8CFA0", 
  "Others"= "grey")

###Pie plot
pie_215DDs_18S <- ggplot(df_genus_abundance_summed, aes(x = 2, y = Summed_Relative_Abundance, fill = Description)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  #Comment here if you want to remove the labels of each portion of the pie.
  geom_text(aes(x = 2, y = Summed_Relative_Abundance, label = paste0(Summed_Relative_Abundance, "%")),
  position = position_stack(vjust = 0.5), size = 1, fontface = "bold", color = "white") +
  #Until here the comment.
  coord_polar(theta = "y", start = -pi/0.6) +  #Move angle of the pie
  theme_void() +
  scale_fill_manual(values = colors) +
  ggtitle("DDs (0-215m)") +
  theme(legend.position = "bottom")
pie_215DDs_18S
ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/Phototroph_pie_215DDs_18S.png",pie_215DDs_18S, width = 9, height = 4)

######
######
#############0-215 nonDDs------
physeq_215nonDDs <- subset_samples(physeq_18, sample_type == "0-215 non-DDs")
#Remove ESVs which have 0 counts in all samples belonging to 0-215 DDs group/
physeq_215nonDDsf <- prune_taxa(taxa_sums(physeq_215nonDDs) > 0, physeq_215nonDDs)
physeq_215nonDDsf
otu_table(physeq_215nonDDsf)
sample_data(physeq_215nonDDsf)
tax_table(physeq_215nonDDsf)
#Extraction of the taxa as data frame from the phyloseq object.
taxa_df <- as.data.frame(tax_table(physeq_215nonDDsf))
#Converting the ESVs into ASVs.
taxa_df <- taxa_df %>%
  rownames_to_column(var = "ASV_ID") %>%      # 1. Crea la columna con los IDs
  mutate(ASV_ID = gsub("ESV", "ASV", ASV_ID)) # 2. Cambia ESV por ASV
#These are the ASVs that needs to be updated as they got selected as the most abundants.
Blast_ASVtaxa <- read.csv("~/Documents/Steve Lab/Svalbard/18S_diversity_2021/Blast_euk_ASV_correction.csv", stringsAsFactors = FALSE)

#Search ASVs and replace them
for (i in 1:nrow(Blast_ASVtaxa)) {
  asv_target <- Blast_ASVtaxa$ASV_ID[i]
  if (asv_target %in% taxa_df$ASV_ID) {
    fila_en_taxa <- which(taxa_df$ASV_ID == asv_target)
    taxa_df[fila_en_taxa, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")] <- 
      Blast_ASVtaxa[i, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
    message("success: ", asv_target)
  } else {
    warning("They are not in taxa_df: ", asv_target)
  }
}
#ASV_48, ASV_1003, ASV_1037 and ASV_1193 are not in taxa_df, it must be bc they are not in these group of samples. That is ok!

rownames(taxa_df) <- taxa_df$ASV_ID
# taxa_df$ASV_ID <- NULL 
taxa_matrix <- as.matrix(taxa_df)
taxa_names(physeq_215nonDDsf) <- gsub("ESV", "ASV", taxa_names(physeq_215nonDDsf))
tax_table(physeq_215nonDDsf) <- tax_table(taxa_matrix)
physeq_215nonDDsf #Verifying
head(tax_table(physeq_215nonDDsf))

####Let's clean 
physeq_215nonDDsf_filtered <- prune_taxa(taxa_sums(physeq_215nonDDsf) > 0, physeq_215nonDDsf)
physeq_215nonDDsf_filtered2 <- prune_samples(sample_sums(physeq_215nonDDsf_filtered) > 0, physeq_215nonDDsf_filtered)
####Let's convert into relative abundances.
ps.rel_18 = transform_sample_counts(physeq_215nonDDsf_filtered2, function(x) x/sum(x)*100)
ps.rel_18_filtered <- prune_taxa(
  !(tax_table(ps.rel_18)[, "Genus"] %in% c("Spumella", "Pedospumella")),
  ps.rel_18)
ps.rel_18_filtered@tax_table #Spumella is a microalgae (golden algae but it is heterotroph)
####Agglomerate taxa
glom <- tax_glom(ps.rel_18_filtered, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
head(ps.melt) #I am doing this to keep the phylum and visually identify the Mosses and microalgae.
#List of phylums to keep
phylums_to_keep <- c("Cryptomonadales", "Florideophycidae", "Dinoflagellata",
                     "Streptophyta","Phragmoplastophyta", "Ochrophyta", "Chlorophyta_ph", 
                     "Chlorophyta")
ps.melt$Phylum <- as.character(ps.melt$Phylum)
ps.melt$Phylum <- ifelse(ps.melt$Phylum %in% phylums_to_keep, ps.melt$Phylum, "Other phylum")
head(ps.melt)
colnames(ps.melt)
ps.melt$Type <- ifelse(ps.melt$Phylum == "Other phylum", "Non-Phototrophs", "Phototrophs")

###Let's check % of relative abundances:
total_abundance <- sum(ps.melt$Abundance)
phototrophs_abundance <- sum(ps.melt$Abundance[ps.melt$Type == "Phototrophs"])
phototrophs_abundance
other_phylum_abundance <- total_abundance - phototrophs_abundance
other_phylum_abundance

phototroph_percentage <- (phototrophs_abundance / total_abundance) * 100
phototroph_percentage
other_phylum_percentage <- (other_phylum_abundance / total_abundance) * 100
other_phylum_percentage

#Mean and SD values:
phototroph_data <- ps.melt %>%
  filter(Type == "Phototrophs")
mean_abundance_phototroph <- mean(phototroph_data$Abundance, na.rm = TRUE) ##Mean
mean_abundance_phototroph
sd_abundance_phototroph <- sd(phototroph_data$Abundance, na.rm = TRUE) ##SD
sd_abundance_phototroph

########Data frame----
data_215nonDDs <- data.frame(
  Phylum = c("Phototrophs", "Non-Phototrophs"),
  Abundance = c(71, 29)  #% of rel abundances
)
data_215nonDDs$angle <- cumsum(data_215nonDDs$Abundance) - 0.5 * data_215nonDDs$Abundance
colors <- c(
  "Phototrophs" = "#F1F5A8",   # azul
  "Non-Phototrophs" = "#FFD0EC"
)

########Pie plot-----
pie_215nonDDs <- ggplot(data_215nonDDs, aes(x = 2, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  geom_text(aes(x = 2, y = Abundance, label = paste0(Abundance, "%")), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y", start = +pi/-5) +  # Girar el donut plot a 90 grados en sentido contrario a las agujas del reloj
  theme_void() +
  xlim(1, 2.5) +
  ggtitle("DDs (0-215m)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = colors) +
  theme(legend.position = "bottom")  # Ajustar la posición de la leyenda si es necesario
pie_215nonDDs 
ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/215nonDDs_Euk_pieplot.png", pie_215nonDDs)

########Genus at the Phylum level and at the last taxon - Phototrophs-----
ps.rel_18_filtered ###Phyloseq of 0-215 nonDDs 
df_physeq <- otu_table(ps.rel_18_filtered)
colnames(df_physeq)
otu_frame <- as.data.frame(df_physeq)
#Add taxonomies as new columns
taxa_df <- as.data.frame(tax_table(ps.rel_18_filtered))
colnames(taxa_df)
df_physeq <- cbind(df_physeq, taxa_df) #Join dataframes 
sample_data_df <- as.data.frame(sample_data(ps.rel_18_filtered))
phylums_to_keep <- c("Cryptomonadales", "Florideophycidae", "Dinoflagellata",
                     "Streptophyta","Phragmoplastophyta", "Ochrophyta", "Chlorophyta_ph",
                     "Chlorophyta")
df_physeq_selected <- subset(df_physeq, Phylum %in% phylums_to_keep)
dim(df_physeq)
# [1] 158  23
dim(df_physeq_selected) #Data with only the phylum we selected
# [1] 58 23
#Abundance:
abundance<- df_physeq_selected[, grepl("^SUN", colnames(df_physeq_selected))]
#Total for each sample:
sum_sample <- rowSums(abundance)
#Relative abundance:
head(df_physeq_selected) #Let's check how it looks like
colnames(df_physeq_selected) #Good
df_physeq_selected$Relative_Abundance <- rowSums(abundance) / sum(sum_sample) * 100 #We added the rel ab column
#Replace NA by "NA"
df_physeq_selected <- df_physeq_selected %>%
  mutate(across(c(Class, Order, Family, Genus), ~replace_na(.x, "NA")))
#Let's do our rank:
df_genus_abundance <- df_physeq_selected %>%
  select(ASV_ID, Relative_Abundance) #This isolated table make us easy to see the highest abundant in the rank
df_genus_abundance$ASV_ID <- ifelse(df_physeq_selected$Relative_Abundance < 0.5, "ASVs < 0.5%", df_physeq_selected$ASV_ID)
df_genus_abundance #Now, ASV <1% will be recognize all the same as a small percentage ASV.
#Make a summary of the ASVs we are going to plot
df_genus_abundance_summed <- df_genus_abundance %>%
  group_by(ASV_ID) %>%
  mutate(Summed_Relative_Abundance = sum(Relative_Abundance)) %>%
  distinct(ASV_ID, .keep_all = TRUE) %>%
  select(ASV_ID, Summed_Relative_Abundance)
df_genus_abundance_summed #This is the list I am going to use.
#Detecting the last taxa level assigned to the ASVs
df_physeq_selected_LastTax <- df_physeq_selected %>%
  mutate(
    Taxon = case_when(
      !is.na(Genus) & Genus != "NA" ~ paste("Genus:", Genus),
      !is.na(Family) & Family != "NA" ~ paste("Family:", Family),
      !is.na(Order) & Order != "NA" ~ paste("Order:", Order),
      !is.na(Class) & Class != "NA" ~ paste("Class:", Class),
      TRUE ~ paste("Phylum:", Phylum)
    )
  )
#Adding this to the list
df_genus_abundance_summed <- df_genus_abundance_summed %>%
  left_join(
    df_physeq_selected_LastTax %>%
      select(ASV_ID, Taxon),
    by = "ASV_ID") %>%
  mutate(Taxon = ifelse(is.na(Taxon), "Others", Taxon))
df_genus_abundance_summed

df_genus_abundance_summed <- df_genus_abundance_summed %>%
  mutate(
    Description = case_when(
      Taxon == "Class: Bryopsida" ~ "Bryopsida (moss)",
      Taxon == "Genus: Atrichopsis sp." ~ "Atrichopsis (mosses)",
      Taxon == "Family: Bryaceae" ~ "Bryaceae (mosses)",
      Taxon == "Genus: Blasia sp." ~ "Blasia (liverworts)",
      Taxon == "Genus: Tortula sp. " ~ "Tortula (moss)",
      Taxon == "Genus: Sanguina sp." ~ "Sanguina (microalgae)",
      Taxon == "Genus: Chloroidium sp." ~ "Chloroidium (microalgae)",
      Taxon == "Genus: Limnomonas sp." ~ "Limnomonas (microalgae)",
      Taxon == "Family: Chlamydomonadales_fa" ~ "Chlamydomonadales (microalgae)",
      Taxon == "Genus: Chloroidium" ~ "Chloroidium (microalgae)",
      Taxon == "Genus: Chlamydomonas" ~ "Chlamydomonadales (microalgae)",
      Taxon == "Others" ~ "Others",
      TRUE ~ Taxon ))
df_genus_abundance_summed


########pieplot-----
#Different colors for each ASV
colors <- c(
  "Bryopsida (moss)" = "#E9A89B",   # azul
  "Atrichopsis (mosses)" = "#B7C9F2",       # naranja
  "Bryaceae (mosses)" = "#EF9C66",       # verde
  "Blasia (liverworts)" = "#78ABA8",
  "Tortula (moss)" = "#FCDC94",
  "Chloroidium (microalgae)" = "#C8CFA0", 
  "Others"= "grey",
  "Chlamydomonadales (microalgae)" = "#BBE9FF",
  "Limnomonas (microalgae)"= "#7C96AB",
  "Sanguina (microalgae)"="#E7D4B5" )

###Pie plot
pie_215nonDDs_18S <- ggplot(df_genus_abundance_summed, aes(x = 2, y = Summed_Relative_Abundance, fill = Description)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  #Comment here if you want to remove the labels of each portion of the pie.
  geom_text(aes(x = 2, y = Summed_Relative_Abundance, label = paste0(Summed_Relative_Abundance, "%")),
            position = position_stack(vjust = 0.5), size = 1, fontface = "bold", color = "white") +
  #Until here the comment.
  coord_polar(theta = "y", start = -pi/0.6) +  #Move angle of the pie
  theme_void() +
  scale_fill_manual(values = colors) +
  ggtitle("Non-DDs (0-215m)") +
  theme(legend.position = "bottom")
pie_215nonDDs_18S

ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/Phototroph_pie_215nonDDs_18S.png",pie_215nonDDs_18S, width = 9, height = 4)

######
######
#############315-850 DDs-----
physeq_850DDs <- subset_samples(physeq_18, sample_type == "315-850 DDs")
#Remove ESVs which have 0 counts in all samples belonging to 0-215 DDs group/
physeq_850DDsf <- prune_taxa(taxa_sums(physeq_850DDs) > 0, physeq_850DDs)
physeq_850DDsf
otu_table(physeq_850DDsf)
sample_data(physeq_850DDsf)
tax_table(physeq_850DDsf)
#Extraction of the taxa as data frame from the phyloseq object.
taxa_df <- as.data.frame(tax_table(physeq_850DDsf))
#Converting the ESVs into ASVs.
taxa_df <- taxa_df %>%
  rownames_to_column(var = "ASV_ID") %>%      # 1. Crea la columna con los IDs
  mutate(ASV_ID = gsub("ESV", "ASV", ASV_ID)) # 2. Cambia ESV por ASV
#These are the ASVs that needs to be updated as they got selected as the most abundants.
Blast_ASVtaxa <- read.csv("~/Documents/Steve Lab/Svalbard/18S_diversity_2021/Blast_euk_ASV_correction.csv", stringsAsFactors = FALSE)

#Search ASVs and replace them
for (i in 1:nrow(Blast_ASVtaxa)) {
  asv_target <- Blast_ASVtaxa$ASV_ID[i]
  if (asv_target %in% taxa_df$ASV_ID) {
    fila_en_taxa <- which(taxa_df$ASV_ID == asv_target)
    taxa_df[fila_en_taxa, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")] <- 
      Blast_ASVtaxa[i, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
    message("success: ", asv_target)
  } else {
    warning("They are not in taxa_df: ", asv_target)
  }
}
#ASV_1003, ASV_1037 , ASV_1192, ASV_1193 are not in taxa_df, it must be bc they are not in these group of samples. That is ok!

rownames(taxa_df) <- taxa_df$ASV_ID
# taxa_df$ASV_ID <- NULL 
taxa_matrix <- as.matrix(taxa_df)
taxa_names(physeq_850DDsf) <- gsub("ESV", "ASV", taxa_names(physeq_850DDsf))
tax_table(physeq_850DDsf) <- tax_table(taxa_matrix)
physeq_850DDsf #Verifying
head(tax_table(physeq_850DDsf))

####Let's clean 
physeq_850DDsf_filtered <- prune_taxa(taxa_sums(physeq_850DDsf) > 0, physeq_850DDsf)
physeq_850DDsf_filtered2 <- prune_samples(sample_sums(physeq_850DDsf_filtered) > 0, physeq_850DDsf_filtered)
####Let's convert into relative abundances.
ps.rel_18 = transform_sample_counts(physeq_850DDsf_filtered2, function(x) x/sum(x)*100)
ps.rel_18_filtered <- prune_taxa(
  !(tax_table(ps.rel_18)[, "Genus"] %in% c("Spumella", "Pedospumella")),
  ps.rel_18)
ps.rel_18_filtered@tax_table #Spumella is a microalgae (golden algae but it is heterotroph)
####Agglomerate taxa
glom <- tax_glom(ps.rel_18_filtered, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
head(ps.melt) #I am doing this to keep the phylum and visually identify the Mosses and microalgae.
#List of phylums to keep
phylums_to_keep <- c("Cryptomonadales", "Florideophycidae", "Dinoflagellata",
                     "Streptophyta","Phragmoplastophyta", "Ochrophyta", "Chlorophyta_ph", 
                     "Chlorophyta")
ps.melt$Phylum <- as.character(ps.melt$Phylum)
ps.melt$Phylum <- ifelse(ps.melt$Phylum %in% phylums_to_keep, ps.melt$Phylum, "Other phylum")
head(ps.melt)
colnames(ps.melt)
ps.melt$Type <- ifelse(ps.melt$Phylum == "Other phylum", "Non-Phototrophs", "Phototrophs")

###Let's check % of relative abundances:
total_abundance <- sum(ps.melt$Abundance)
phototrophs_abundance <- sum(ps.melt$Abundance[ps.melt$Type == "Phototrophs"])
phototrophs_abundance
other_phylum_abundance <- total_abundance - phototrophs_abundance
other_phylum_abundance

phototroph_percentage <- (phototrophs_abundance / total_abundance) * 100
phototroph_percentage
other_phylum_percentage <- (other_phylum_abundance / total_abundance) * 100
other_phylum_percentage

#Mean and SD values:
phototroph_data <- ps.melt %>%
  filter(Type == "Phototrophs")
mean_abundance_phototroph <- mean(phototroph_data$Abundance, na.rm = TRUE) ##Mean
mean_abundance_phototroph
sd_abundance_phototroph <- sd(phototroph_data$Abundance, na.rm = TRUE) ##SD
sd_abundance_phototroph

########Data frame-----
data_850DDs <- data.frame(
  Phylum = c("Phototrophs", "Non-Phototrophs"),
  Abundance = c(64, 36)  #% of rel abundances
)
data_850DDs$angle <- cumsum(data_850DDs$Abundance) - 0.5 * data_850DDs$Abundance
colors <- c(
  "Phototrophs" = "#F1F5A8",   # azul
  "Non-Phototrophs" = "#FFD0EC"
)

########Pie plot -----
pie_850DDs <- ggplot(data_850DDs, aes(x = 2, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  geom_text(aes(x = 2, y = Abundance, label = paste0(Abundance, "%")), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y", start = +pi/1.18) +  # Girar el donut plot a 90 grados en sentido contrario a las agujas del reloj
  theme_void() +
  xlim(1, 2.5) +
  ggtitle("DDs (315-850m)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = colors) +
  theme(legend.position = "bottom")  # Ajustar la posición de la leyenda si es necesario
pie_850DDs 
ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/850DDs_Euk_pieplot.png", pie_850DDs)

########Genus at the Phylum level and at the last taxon - Phototrophs ----
ps.rel_18_filtered ###Phyloseq of 315-850
df_physeq <- otu_table(ps.rel_18_filtered)
colnames(df_physeq)
otu_frame <- as.data.frame(df_physeq)
#Add taxonomies as new columns
taxa_df <- as.data.frame(tax_table(ps.rel_18_filtered))
colnames(taxa_df)
df_physeq <- cbind(df_physeq, taxa_df) #Join dataframes 
sample_data_df <- as.data.frame(sample_data(ps.rel_18_filtered))
phylums_to_keep <- c("Cryptomonadales", "Florideophycidae", "Dinoflagellata",
                     "Streptophyta","Phragmoplastophyta", "Ochrophyta", "Chlorophyta_ph",
                     "Chlorophyta")
df_physeq_selected <- subset(df_physeq, Phylum %in% phylums_to_keep)
dim(df_physeq)
# [1] 246  19
dim(df_physeq_selected) #Data with only the phylum we selected
# [1] 72 19
#Abundance:
abundance<- df_physeq_selected[, grepl("^SUN", colnames(df_physeq_selected))]
#Total for each sample:
sum_sample <- rowSums(abundance)
#Relative abundance:
head(df_physeq_selected) #Let's check how it looks like
colnames(df_physeq_selected) #Good
df_physeq_selected$Relative_Abundance <- rowSums(abundance) / sum(sum_sample) * 100 #We added the rel ab column
#Replace NA by "NA"
df_physeq_selected <- df_physeq_selected %>%
  mutate(across(c(Class, Order, Family, Genus), ~replace_na(.x, "NA")))
#Let's do our rank:
df_genus_abundance <- df_physeq_selected %>%
  select(ASV_ID, Relative_Abundance) #This isolated table make us easy to see the highest abundant in the rank
df_genus_abundance$ASV_ID <- ifelse(df_physeq_selected$Relative_Abundance < 0.5, "ASVs < 0.5%", df_physeq_selected$ASV_ID)
df_genus_abundance #Now, ASV <1% will be recognize all the same as a small percentage ASV.
#Make a summary of the ASVs we are going to plot
df_genus_abundance_summed <- df_genus_abundance %>%
  group_by(ASV_ID) %>%
  mutate(Summed_Relative_Abundance = sum(Relative_Abundance)) %>%
  distinct(ASV_ID, .keep_all = TRUE) %>%
  select(ASV_ID, Summed_Relative_Abundance)
df_genus_abundance_summed #This is the list I am going to use.
#Detecting the last taxa level assigned to the ASVs
df_physeq_selected_LastTax <- df_physeq_selected %>%
  mutate(
    Taxon = case_when(
      !is.na(Genus) & Genus != "NA" ~ paste("Genus:", Genus),
      !is.na(Family) & Family != "NA" ~ paste("Family:", Family),
      !is.na(Order) & Order != "NA" ~ paste("Order:", Order),
      !is.na(Class) & Class != "NA" ~ paste("Class:", Class),
      TRUE ~ paste("Phylum:", Phylum)
    )
  )
#Adding this to the list
df_genus_abundance_summed <- df_genus_abundance_summed %>%
  left_join(
    df_physeq_selected_LastTax %>%
      select(ASV_ID, Taxon),
    by = "ASV_ID") %>%
  mutate(Taxon = ifelse(is.na(Taxon), "Others", Taxon))
df_genus_abundance_summed

df_genus_abundance_summed <- df_genus_abundance_summed %>%
  mutate(
    Description = case_when(
      Taxon == "Class: Bryopsida" ~ "Bryopsida (moss)",
      Taxon == "Genus: Atrichopsis sp." ~ "Atrichopsis (mosses)",
      Taxon == "Family: Bryaceae" ~ "Bryaceae (mosses)",
      Taxon == "Genus: Blasia sp." ~ "Blasia (liverworts)",
      Taxon == "Genus: Tortula sp. " ~ "Tortula (moss)",
      Taxon == "Genus: Sanguina sp." ~ "Sanguina (microalgae)",
      Taxon == "Genus: Chloroidium sp." ~ "Chloroidium (microalgae)",
      Taxon == "Genus: Limnomonas sp." ~ "Limnomonas (microalgae)",
      Taxon == "Family: Chlamydomonadales_fa" ~ "Chlamydomonadales (microalgae)",
      Taxon == "Genus: Chloroidium" ~ "Chloroidium (microalgae)",
      Taxon == "Genus: Chlamydomonas" ~ "Chlamydomonadales (microalgae)",
      Taxon == "Genus: Chromochloris sp." ~ "Chromochloris (microalgae)",
      Taxon == "Genus: Xanthonema" ~ "Xanthonema (microalgae)",
      Taxon == "Genus: Ochromonas" ~ "Ochromonas (microalgae)",
      Taxon == "Genus: Coccomyxa" ~ "Coccomyxa (microalgae)",
      Taxon == "Others" ~ "Others",
      TRUE ~ Taxon ))
df_genus_abundance_summed

########pieplot ----
colors <- c(
  "Bryopsida (moss)" = "#E9A89B",   # azul
  "Atrichopsis (mosses)" = "#B7C9F2",       # naranja
  "Bryaceae (mosses)" = "#EF9C66",       # verde
  "Blasia (liverworts)" = "#78ABA8",
  "Tortula (moss)" = "#FCDC94",
  "Chloroidium (microalgae)" = "#C8CFA0", 
  "Others"= "grey",
  "Chlamydomonadales (microalgae)" = "#BBE9FF",
  "Limnomonas (microalgae)"= "#7C96AB",
  "Sanguina (microalgae)"="#E7D4B5" ,
  "Xanthonema (microalgae)" = "#25671E",
  "Coccomyxa (microalgae)"= "#FE9EC7")
###Pie plot
pie_850DDs_18S <- ggplot(df_genus_abundance_summed, aes(x = 2, y = Summed_Relative_Abundance, fill = Description)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  #Comment here if you want to remove the labels of each portion of the pie.
  geom_text(aes(x = 2, y = Summed_Relative_Abundance, label = paste0(Summed_Relative_Abundance, "%")),
            position = position_stack(vjust = 0.5), size = 1, fontface = "bold", color = "white") +
  #Until here the comment.
  coord_polar(theta = "y", start = -pi/0.6) +  #Move angle of the pie
  theme_void() +
  scale_fill_manual(values = colors) +
  ggtitle("DDs (315-850m)") +
  theme(legend.position = "bottom")
pie_850DDs_18S

ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/Phototroph_pie_850DDs_18S.png",pie_850DDs_18S, width = 9, height = 4)

######
######
#############315-850 nonDDs----
physeq_850nonDDs <- subset_samples(physeq_18, sample_type == "315-850 non-DDs")
#Remove ESVs which have 0 counts in all samples belonging to 0-215 DDs group/
physeq_850nonDDsf <- prune_taxa(taxa_sums(physeq_850nonDDs) > 0, physeq_850nonDDs)
physeq_850nonDDsf
otu_table(physeq_850nonDDsf)
sample_data(physeq_850nonDDsf)
tax_table(physeq_850nonDDsf)
#Extraction of the taxa as data frame from the phyloseq object.
taxa_df <- as.data.frame(tax_table(physeq_850nonDDsf))
#Converting the ESVs into ASVs.
taxa_df <- taxa_df %>%
  rownames_to_column(var = "ASV_ID") %>%      # 1. Crea la columna con los IDs
  mutate(ASV_ID = gsub("ESV", "ASV", ASV_ID)) # 2. Cambia ESV por ASV
#These are the ASVs that needs to be updated as they got selected as the most abundants.
Blast_ASVtaxa <- read.csv("~/Documents/Steve Lab/Svalbard/18S_diversity_2021/Blast_euk_ASV_correction.csv", stringsAsFactors = FALSE)

#Search ASVs and replace them
for (i in 1:nrow(Blast_ASVtaxa)) {
  asv_target <- Blast_ASVtaxa$ASV_ID[i]
  if (asv_target %in% taxa_df$ASV_ID) {
    fila_en_taxa <- which(taxa_df$ASV_ID == asv_target)
    taxa_df[fila_en_taxa, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")] <- 
      Blast_ASVtaxa[i, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")]
    message("success: ", asv_target)
  } else {
    warning("They are not in taxa_df: ", asv_target)
  }
}
#ASV_1037, ASV_1192 and ASV_1193 are not in taxa_df, it must be bc they are not in these group of samples. That is ok!

rownames(taxa_df) <- taxa_df$ASV_ID
# taxa_df$ASV_ID <- NULL 
taxa_matrix <- as.matrix(taxa_df)
taxa_names(physeq_850nonDDsf) <- gsub("ESV", "ASV", taxa_names(physeq_850nonDDsf))
tax_table(physeq_850nonDDsf) <- tax_table(taxa_matrix)
physeq_850nonDDsf #Verifying
head(tax_table(physeq_850nonDDsf))

####Let's clean 
physeq_850nonDDsf_filtered <- prune_taxa(taxa_sums(physeq_850nonDDsf) > 0, physeq_850nonDDsf)
physeq_850nonDDsf_filtered2 <- prune_samples(sample_sums(physeq_850nonDDsf_filtered) > 0, physeq_850nonDDsf_filtered)
####Let's convert into relative abundances.
ps.rel_18 = transform_sample_counts(physeq_850nonDDsf_filtered2, function(x) x/sum(x)*100)
ps.rel_18_filtered <- prune_taxa(
  !(tax_table(ps.rel_18)[, "Genus"] %in% c("Spumella", "Pedospumella")),
  ps.rel_18)
ps.rel_18_filtered@tax_table #Spumella is a microalgae (golden algae but it is heterotroph)
####Agglomerate taxa
glom <- tax_glom(ps.rel_18_filtered, taxrank = 'Phylum', NArm = FALSE)
ps.melt <- psmelt(glom)
head(ps.melt) #I am doing this to keep the phylum and visually identify the Mosses and microalgae.
#List of phylums to keep
phylums_to_keep <- c("Cryptomonadales", "Florideophycidae", "Dinoflagellata",
                     "Streptophyta","Phragmoplastophyta", "Ochrophyta", "Chlorophyta_ph", 
                     "Chlorophyta")
ps.melt$Phylum <- as.character(ps.melt$Phylum)
ps.melt$Phylum <- ifelse(ps.melt$Phylum %in% phylums_to_keep, ps.melt$Phylum, "Other phylum")
head(ps.melt)
colnames(ps.melt)
ps.melt$Type <- ifelse(ps.melt$Phylum == "Other phylum", "Non-Phototrophs", "Phototrophs")

###Let's check % of relative abundances:
total_abundance <- sum(ps.melt$Abundance)
phototrophs_abundance <- sum(ps.melt$Abundance[ps.melt$Type == "Phototrophs"])
phototrophs_abundance
other_phylum_abundance <- total_abundance - phototrophs_abundance
other_phylum_abundance

phototroph_percentage <- (phototrophs_abundance / total_abundance) * 100
phototroph_percentage
other_phylum_percentage <- (other_phylum_abundance / total_abundance) * 100
other_phylum_percentage

#Mean and SD values:
phototroph_data <- ps.melt %>%
  filter(Type == "Phototrophs")
mean_abundance_phototroph <- mean(phototroph_data$Abundance, na.rm = TRUE) ##Mean
mean_abundance_phototroph
sd_abundance_phototroph <- sd(phototroph_data$Abundance, na.rm = TRUE) ##SD
sd_abundance_phototroph



########Data frame------
data_850nonDDs <- data.frame(
  Phylum = c("Phototrophs", "Non-Phototrophs"),
  Abundance = c(64, 36))  #% of rel abundances
data_850nonDDs$angle <- cumsum(data_850nonDDs$Abundance) - 0.5 * data_850nonDDs$Abundance
colors <- c(
  "Phototrophs" = "#F1F5A8",   # azul
  "Non-Phototrophs" = "#FFD0EC")

########Pie plot----
pie_850nonDDs <- ggplot(data_850nonDDs, aes(x = 2, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  geom_text(aes(x = 2, y = Abundance, label = paste0(Abundance, "%")), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y", start = +pi/1.18) +  # Girar el donut plot a 90 grados en sentido contrario a las agujas del reloj
  theme_void() +
  xlim(1, 2.5) +
  ggtitle("DDs (315-850m)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = colors) +
  theme(legend.position = "bottom")  # Ajustar la posición de la leyenda si es necesario
pie_850nonDDs 
ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/850nonDDs_Euk_pieplot.png", pie_850nonDDs)

########Genus at the Phylum level and at the last taxon - Phototrophs-----
ps.rel_18_filtered ###Phyloseq of 315-850
df_physeq <- otu_table(ps.rel_18_filtered)
colnames(df_physeq)
otu_frame <- as.data.frame(df_physeq)
#Add taxonomies as new columns
taxa_df <- as.data.frame(tax_table(ps.rel_18_filtered))
colnames(taxa_df)
df_physeq <- cbind(df_physeq, taxa_df) #Join dataframes 
sample_data_df <- as.data.frame(sample_data(ps.rel_18_filtered))
phylums_to_keep <- c("Cryptomonadales", "Florideophycidae", "Dinoflagellata",
                     "Streptophyta","Phragmoplastophyta", "Ochrophyta", "Chlorophyta_ph",
                     "Chlorophyta")
df_physeq_selected <- subset(df_physeq, Phylum %in% phylums_to_keep)
dim(df_physeq)
# [1] 222  18
dim(df_physeq_selected) #Data with only the phylum we selected
# [1] 71 18
#Abundance:
abundance<- df_physeq_selected[, grepl("^SUN", colnames(df_physeq_selected))]
#Total for each sample:
sum_sample <- rowSums(abundance)
#Relative abundance:
head(df_physeq_selected) #Let's check how it looks like
colnames(df_physeq_selected) #Good
df_physeq_selected$Relative_Abundance <- rowSums(abundance) / sum(sum_sample) * 100 #We added the rel ab column
#Replace NA by "NA"
df_physeq_selected <- df_physeq_selected %>%
  mutate(across(c(Class, Order, Family, Genus), ~replace_na(.x, "NA")))
#Let's do our rank:
df_genus_abundance <- df_physeq_selected %>%
  select(ASV_ID, Relative_Abundance) #This isolated table make us easy to see the highest abundant in the rank
df_genus_abundance$ASV_ID <- ifelse(df_physeq_selected$Relative_Abundance < 0.5, "ASVs < 0.5%", df_physeq_selected$ASV_ID)
df_genus_abundance #Now, ASV <1% will be recognize all the same as a small percentage ASV.
#Make a summary of the ASVs we are going to plot
df_genus_abundance_summed <- df_genus_abundance %>%
  group_by(ASV_ID) %>%
  mutate(Summed_Relative_Abundance = sum(Relative_Abundance)) %>%
  distinct(ASV_ID, .keep_all = TRUE) %>%
  select(ASV_ID, Summed_Relative_Abundance)
df_genus_abundance_summed #This is the list I am going to use.
#Detecting the last taxa level assigned to the ASVs
df_physeq_selected_LastTax <- df_physeq_selected %>%
  mutate(
    Taxon = case_when(
      !is.na(Genus) & Genus != "NA" ~ paste("Genus:", Genus),
      !is.na(Family) & Family != "NA" ~ paste("Family:", Family),
      !is.na(Order) & Order != "NA" ~ paste("Order:", Order),
      !is.na(Class) & Class != "NA" ~ paste("Class:", Class),
      TRUE ~ paste("Phylum:", Phylum)
    )
  )
#Adding this to the list
df_genus_abundance_summed <- df_genus_abundance_summed %>%
  left_join(
    df_physeq_selected_LastTax %>%
      select(ASV_ID, Taxon),
    by = "ASV_ID") %>%
  mutate(Taxon = ifelse(is.na(Taxon), "Others", Taxon))
df_genus_abundance_summed

df_genus_abundance_summed <- df_genus_abundance_summed %>%
  mutate(
    Description = case_when(
      Taxon == "Class: Bryopsida" ~ "Bryopsida (moss)",
      Taxon == "Genus: Atrichopsis sp." ~ "Atrichopsis (mosses)",
      Taxon == "Family: Bryaceae" ~ "Bryaceae (mosses)",
      Taxon == "Genus: Blasia sp." ~ "Blasia (liverworts)",
      Taxon == "Genus: Tortula sp. " ~ "Tortula (moss)",
      Taxon == "Genus: Sanguina sp." ~ "Sanguina (microalgae)",
      Taxon == "Genus: Chloroidium sp." ~ "Chloroidium (microalgae)",
      Taxon == "Genus: Limnomonas sp." ~ "Limnomonas (microalgae)",
      Taxon == "Family: Chlamydomonadales_fa" ~ "Chlamydomonadales (microalgae)",
      Taxon == "Genus: Chloroidium" ~ "Chloroidium (microalgae)",
      Taxon == "Genus: Chlamydomonas" ~ "Chlamydomonadales (microalgae)",
      Taxon == "Genus: Chromochloris sp." ~ "Chromochloris (microalgae)",
      Taxon == "Genus: Xanthonema" ~ "Xanthonema (microalgae)",
      Taxon == "Genus: Ochromonas" ~ "Ochromonas (microalgae)",
      Taxon == "Genus: Coccomyxa" ~ "Coccomyxa (microalgae)",
      Taxon == "Others" ~ "Others",
      TRUE ~ Taxon ))
df_genus_abundance_summed

########Pieplot-----
colors <- c(
  "Bryopsida (moss)" = "#E9A89B",   # azul
  "Atrichopsis (mosses)" = "#B7C9F2",       # naranja
  "Bryaceae (mosses)" = "#EF9C66",       # verde
  "Blasia (liverworts)" = "#78ABA8",
  "Tortula (moss)" = "#FCDC94",
  "Chloroidium (microalgae)" = "#C8CFA0", 
  "Others"= "grey",
  "Chlamydomonadales (microalgae)" = "#BBE9FF",
  "Limnomonas (microalgae)"= "#7C96AB",
  "Sanguina (microalgae)"="#E7D4B5" ,
  "Xanthonema (microalgae)" = "#25671E",
  "Coccomyxa (microalgae)"= "#FE9EC7")

###Pie plot
pie_850nonDDs_18S <- ggplot(df_genus_abundance_summed, aes(x = 2, y = Summed_Relative_Abundance, fill = Description)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  #Comment here if you want to remove the labels of each portion of the pie.
  # geom_text(aes(x = 2, y = Summed_Relative_Abundance, label = paste0(Summed_Relative_Abundance, "%")),
  #           position = position_stack(vjust = 0.5), size = 1, fontface = "bold", color = "white") +
  #Until here the comment.
  coord_polar(theta = "y", start = -pi/0.6) +  #Move angle of the pie
  theme_void() +
  scale_fill_manual(values = colors) +
  ggtitle("Non-DDs (315-850m)") +
  theme(legend.position = "bottom")
pie_850nonDDs_18S
ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/Phototroph_pie_850nonDDs_18S.png",pie_850nonDDs_18S, width = 9, height = 4)
