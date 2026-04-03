#Load libraries
library(dplyr)
# install.packages("ggbreak")
library(ggbreak)
library(patchwork)

#Lets get the abundance:
esv_abundance <- physeq_18S %>%
  psmelt()
#format
head(esv_abundance)
dim(esv_abundance)  #[1] 41965    12
#Filtering just samples 0-250 on stripe:
head(esv_abundance)
colnames(esv_abundance)
esv_abundance$OTU <- gsub("^ESV_", "ASV_", esv_abundance$OTU)


###ON Stripe - 0-250m-------

#Selecting samples that just belong to 0-250 stripe
subset_on250 <- esv_abundance[esv_abundance$sample_type == "0-215 DDs", ]
subset_on250 <- as.data.frame(subset_on250)
#Now, removing the ESVs which doesn't belong to this sample type:
esv_unique_on250 <- subset_on250 %>%
  group_by(OTU) %>%
  filter(sum(Abundance) > 0) %>%
  ungroup()

head(esv_unique_on250)
#Relative abundance per sample:
rel_ab_on250 <- esv_unique_on250 %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = (Abundance / sum(Abundance)) * 100)
head(rel_ab_on250)

write.csv(rel_ab_on250, "rel_ab_on250.csv", row.names = FALSE)


#EL SUN0102 esta con NAN, no se porq.
#Remove other columns we don't care:
rel_ab_on250 <- rel_ab_on250 %>%
  dplyr::select(Relative_Abundance, OTU, Sample, Genus)
head(rel_ab_on250)

#In order to probe the rel abun, let's take a sample, sum up all the ESVs and the total should be 100:
sample_SUN0022 <- subset(rel_ab_on250, Sample == "SUN0022")
total_abundance_SUN0022 <- sum(sample_SUN0022$Relative_Abundance)
total_abundance_SUN0022 #it is 100%! it is ok!

#deleting NA 
#(I have noticed that I have some NA in my rela abundance, 
#but it should be a number, so I removed it)
#The explanation should be that those ESVs doesn't belong to that sample group.
rel_ab_on250 <- rel_ab_on250 %>%
  filter(!is.na(Relative_Abundance))


############################## OUTLIERS######
##### Remove outliers based in boxplots (It is optional, 
#sometimes it is better to keep the outliers) #Better decide based on a specific analysis.
## Start here
cleaned_dataframe <- rel_ab_on250 %>%
  group_by(OTU) %>%
  mutate(
    lower_threshold = quantile(Relative_Abundance, 0.25) - 1.5 * IQR(Relative_Abundance, na.rm = TRUE),
    upper_threshold = quantile(Relative_Abundance, 0.75) + 1.5 * IQR(Relative_Abundance, na.rm = TRUE)
  ) %>%
  filter(Relative_Abundance >= lower_threshold & Relative_Abundance <= upper_threshold) %>%
  ungroup()
rel_ab_on250 <- cleaned_dataframe
### end here
##############################end of outliers

###Order by mean:--------------------------------
summary_on250_mean <- aggregate(Relative_Abundance ~ OTU, data = rel_ab_on250, FUN = mean) %>%
  arrange(desc(Relative_Abundance))
#rename column
colnames(summary_on250_mean) <- c("OTU", "mean_abundance")

head(summary_on250_mean)
#getting the top 10
top10_on250_mean <- summary_on250_mean %>%
  arrange(desc(mean_abundance)) %>%
  head(10)
top10_on250_mean
#Let's rename our data where all the poinst or observations per sample are present
#All data points belong to the top we referred prevously.
subset_on250_mean <- rel_ab_on250 %>%
  filter(OTU %in% top10_on250_mean$OTU)
subset_on250_mean

#Create a blast vector.
blast_vec <- c(
  "ASV_1"  = "Bryopsida (Class)",
  "ASV_7" = "Blasia sp.",
  "ASV_3" = "Bryaceae (Family)",
  # "ASV_9" = "Tortula sp.",
  # "ASV_11"= "Chloroidium sp.",
  "ASV_10" = "Pragmopora sp.",
  "ASV_1003"= "Mucor sp.",
  "ASV_1037"= "Poteriospumella sp.",
  "ASV_105"= "Basidiomycota (Phylum)",
  "ASV_102"= "Sphaeropleales (Order)",
  "ASV_121" = "Basidiomycota (Phylum)",
  "ASV_139" = "Malassezia sp.")

#I went to the rep set file and took the ASV sequences and blasted them
subset_on250_mean$Blast_taxa <- blast_vec[subset_on250_mean$OTU]

##Box plot:----
#color:
genus_colors <- c("Bryopsida (Class)" = "#FF9800", 
                  "Blasia sp." = "#D1BB9E", 
                  "Bryaceae (Family)" = "#DC0083", 
                  "Tortula sp." = "#7E8EF1", 
                  "Chloroidium sp." = "#D6D46D",
                  "Pragmopora sp." = "#A94438",
                  "Mucor sp."= "#0F67B1",
                  "Poteriospumella sp." = "#03346E",
                  "Xanthonema sp." = "#F075AA",
                  "Tremella sp." = "grey",
                  "Atrichopsis sp."= "#6EACDA",
                  "Limnomonas sp."= "#36BA98",
                  "Rhogostoma sp." = "#FFA38F",
                  "Tetracladium sp." = "#FFC7ED",
                  "Rhogostomidae (Family)" = "#179BAE",
                  "Sanguina sp." = "#e3e656",
                  "Exophiala sp." = "#4C3BCF",
                  "Chromochloris sp."= "#FF204E",
                  "Chaetothyriales (Order)"= "purple",
                  "NA" = "#3B3C36",
                  "Malassezia sp." = "#F1B42F",
                  "Basidiomycota (Phylum)" = "#81613E", 
                  "Sphaeropleales (Order)"="#9FCB98", 
                  "Pottia sp." = "#C9BEFF")


# Crear el gráfico
boxplot_on250_mean <- ggplot(subset_on250_mean, aes(x = factor(OTU, levels = top10_on250_mean$OTU), 
                                                    y = Relative_Abundance, 
                                                    color = Blast_taxa)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_point(data = top10_on250_mean, aes(x = factor(OTU), y = mean_abundance), 
             shape = 18, size = 3, color = "red") +
  scale_y_continuous(limits = c(0, 101), breaks = c(0, 10, 20, 30, 40, 100)) + 
  scale_color_manual(values = genus_colors, name = "Taxa") +
  labs(x = "ASV", y = "Relative Abundance %") +
  ggtitle("A) DDs (0-215m)")+
  scale_y_break(c(41, 95), scales = 0.18) + 
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    panel.border = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(colour = "white"),
    axis.line.y = element_line(colour = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right")
print(boxplot_on250_mean)
#Save as png
ggsave("Relative Abudance_ranked_median and mean/Rel_Ab_ESV/On250_ESV18S.png", boxplot_on250_mean, width = 6, height = 4, dpi = 500)

###OFF Stripe - 0-250m-------

#Selecting samples that just belong to 0-250 stripe
subset_off250 <- esv_abundance[esv_abundance$sample_type == "0-215 non-DDs", ]
subset_off250 <- as.data.frame(subset_off250)
#Now, removing the ESVs which doesn't belong to this sample type:
esv_unique_off250 <- subset_off250 %>%
  group_by(OTU) %>%
  filter(sum(Abundance) > 0) %>%
  ungroup()
#Relative abundance per sample:
rel_ab_off250 <- esv_unique_off250 %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = (Abundance / sum(Abundance)) * 100)
head(rel_ab_off250)
#Remove other columns we don't care:
rel_ab_off250 <- rel_ab_off250 %>%
  dplyr::select(Relative_Abundance, OTU, Sample, Genus)
head(rel_ab_off250)
#In order to probe the rel abun, let's take a sample, sum up all the ESVs and the total should be 100:
sample_SUN0072 <- subset(rel_ab_off250, Sample == "SUN0072")
total_abundance_SUN0072 <- sum(sample_SUN0072$Relative_Abundance)
total_abundance_SUN0072 #it is 100%! it is ok!

#deleting NA 
#(I have noticed that I have some NA in my rela abundance, 
#but it should be a number, so I removed it)
#The explanation should be that those ESVs doesn't belong to that sample group.
rel_ab_off250 <- rel_ab_off250 %>%
  filter(!is.na(Relative_Abundance))

############################## OUTLIERS######
##### Remove outliers based in boxplots (It is optional, 
#sometimes it is better to keep the outliers) #Better decide based on a specific analysis.
## Start here
cleaned_dataframe <- rel_ab_off250 %>%
  group_by(OTU) %>%
  mutate(
    lower_threshold = quantile(Relative_Abundance, 0.25) - 1.5 * IQR(Relative_Abundance, na.rm = TRUE),
    upper_threshold = quantile(Relative_Abundance, 0.75) + 1.5 * IQR(Relative_Abundance, na.rm = TRUE)
  ) %>%
  filter(Relative_Abundance >= lower_threshold & Relative_Abundance <= upper_threshold) %>%
  ungroup()
rel_ab_off250 <- cleaned_dataframe
### end here
##############################end of outliers


###Order by mean:--------------------------------
summary_off250_mean <- aggregate(Relative_Abundance ~ OTU, data = rel_ab_off250, FUN = mean) %>%
  arrange(desc(Relative_Abundance))
#rename column
colnames(summary_off250_mean) <- c("OTU", "mean_abundance")

#getting the top 10
top10_off250_mean <- summary_off250_mean %>%
  arrange(desc(mean_abundance)) %>%
  head(10)
top10_off250_mean #See the taxa names of the ASVs...in order to verify the blast assignation. 

#Let's rename our data where all the poinst or observations per sample are present
#All data points belong to the top we referred prevously.
subset_off250_mean <- rel_ab_off250 %>%
  filter(OTU %in% top10_off250_mean$OTU)
subset_off250_mean

blast_vec <- c(
  "ASV_3" = "Bryaceae (Family)",
  "ASV_1"  = "Bryopsida (Class)",
  "ASV_4" = "Atrichopsis sp.",
  "ASV_2" = "Atrichopsis sp.",
  "ASV_13" = "Tetracladium sp.",
  "ASV_53"= "Limnomonas sp.",
  "ASV_59" = "Rhogostomidae (Family)",
  "ASV_44"= "Sanguina sp.",
  "ASV_7"= "Blasia sp.",
  "ASV_16"= "Rhogostoma sp.",
  "ASV_9" = "Pottia sp.")

subset_off250_mean$Blast_taxa <- blast_vec[subset_off250_mean$OTU]

##Box plot:----
#color
genus_colors <- c("Bryopsida (Class)" = "#FF9800", 
                  "Blasia sp." = "#D1BB9E", 
                  "Bryaceae (Family)" = "#DC0083", 
                  "Tortula sp." = "#7E8EF1", 
                  "Chloroidium sp." = "#D6D46D",
                  "Pragmopora sp." = "#A94438",
                  "Mucor sp."= "#0F67B1",
                  "Poteriospumella sp." = "#03346E",
                  "Xanthonema sp." = "#F075AA",
                  "Tremella sp." = "grey",
                  "Atrichopsis sp."= "#6EACDA",
                  "Limnomonas sp."= "#36BA98",
                  "Rhogostoma sp." = "#FFA38F",
                  "Tetracladium sp." = "#FFC7ED",
                  "Rhogostomidae (Family)" = "#179BAE",
                  "Sanguina sp." = "#e3e656",
                  "Exophiala sp." = "#4C3BCF",
                  "Chromochloris sp."= "#FF204E",
                  "Chaetothyriales (Order)"= "purple",
                  "NA" = "#3B3C36",
                  "Malassezia sp." = "#F1B42F",
                  "Basidiomycota (Phylum)" = "#81613E", 
                  "Sphaeropleales (Order)"="#9FCB98", 
                  "Pottia sp." = "#C9BEFF")

boxplot_off250_mean <- ggplot(subset_off250_mean, aes(x = factor(OTU, levels = top10_off250_mean$OTU), y = Relative_Abundance, color = Blast_taxa)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5)  +
  ggtitle("B) non-DDs (0-215m)") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = genus_colors, name = "Taxa") +  # Asigna los colores definidos manualmente
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black")) +
  theme(axis.title.x = element_text(angle = 0, size = 12)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 13, color = "black")) +
  theme(axis.title.y = element_text(size = 14)) +
  labs(x = "ASV", y = "Relative Abundance %") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(face = "bold", size = 16))  # Título en negrita
boxplot_off250_mean <- boxplot_off250_mean +
  geom_point(data = top10_off250_mean, aes(x = factor(OTU), y = mean_abundance), shape = 18, size = 3, color = "red")
boxplot_off250_mean
#Save as png
ggsave("Relative Abudance_ranked_median and mean/Rel_Ab_ESV/Off250_ESV18S.png", boxplot_off250_mean, width = 6, height = 4, dpi = 500)

###ON Stripe - 300-850m-------

#Selecting samples that just belong to 0-250 stripe
subset_on850 <- esv_abundance[esv_abundance$sample_type == "315-850 DDs", ]
subset_on850 <- as.data.frame(subset_on850)
#Now, removing the ESVs which doesn't belong to this sample type:
esv_unique_on850 <- subset_on850 %>%
  group_by(OTU) %>%
  filter(sum(Abundance) > 0) %>%
  ungroup()
#Relative abundance per sample:
rel_ab_on850 <- esv_unique_on850 %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = (Abundance / sum(Abundance)) * 100)
head(rel_ab_on850)
#Remove other columns we don't care:
rel_ab_on850 <- rel_ab_on850 %>%
  dplyr::select(Relative_Abundance, OTU, Sample, Genus)
head(rel_ab_on850)
#In order to probe the rel abun, let's take a sample, sum up all the ESVs and the total should be 100:
sample_SUN0181 <- subset(rel_ab_on850, Sample == "SUN0181")
total_abundance_SUN0181 <- sum(sample_SUN0181$Relative_Abundance)
total_abundance_SUN0181 #it is 100%! it is ok!

#deleting NA 
#(I have noticed that I have some NA in my rela abundance, 
#but it should be a number, so I removed it)
#The explanation should be that those ESVs doesn't belong to that sample group.
rel_ab_on850 <- rel_ab_on850 %>%
  filter(!is.na(Relative_Abundance))

############################## OUTLIERS######
##### Remove outliers based in boxplots (It is optional, 
#sometimes it is better to keep the outliers) #Better decide based on a specific analysis.
## Start here
cleaned_dataframe <- rel_ab_on850 %>%
  group_by(OTU) %>%
  mutate(
    lower_threshold = quantile(Relative_Abundance, 0.25) - 1.5 * IQR(Relative_Abundance, na.rm = TRUE),
    upper_threshold = quantile(Relative_Abundance, 0.75) + 1.5 * IQR(Relative_Abundance, na.rm = TRUE)
  ) %>%
  filter(Relative_Abundance >= lower_threshold & Relative_Abundance <= upper_threshold) %>%
  ungroup()
rel_ab_on850 <- cleaned_dataframe
### end here
##############################end of outliers


###Order by mean:--------------------------------
summary_on850_mean <- aggregate(Relative_Abundance ~ OTU, data = rel_ab_on850, FUN = mean) %>%
  arrange(desc(Relative_Abundance))
#rename column
colnames(summary_on850_mean) <- c("OTU", "mean_abundance")

#getting the top 10
top10_on850_mean <- summary_on850_mean %>%
  arrange(desc(mean_abundance)) %>%
  head(10)

#Let's rename our data where all the poinst or observations per sample are present
#All data points belong to the top we referred prevously.
subset_on850_mean <- rel_ab_on850 %>%
  filter(OTU %in% top10_on850_mean$OTU)
subset_on850_mean

blast_vec <- c(
  "ASV_2" = "Atrichopsis sp.",
  "ASV_5" = "Exophiala sp.",
  "ASV_1"  = "Bryopsida (Class)",
  "ASV_3" = "Bryaceae (Family)",
  "ASV_4" = "Atrichopsis sp.",
  "ASV_16"= "Rhogostoma sp.",
  "ASV_11" = "Chloroidium sp.",
  "ASV_10" = "Pragmopora sp.",
  "ASV_48"= "Chromochloris sp.",
  "ASV_7"= "Blasia sp.", 
  "ASV_8" = "Chaetothyriales (Order)",
  "ASV_25" = "Chaetothyriales (Order)"
)

subset_on850_mean$Blast_taxa <- blast_vec[subset_on850_mean$OTU]

##Box plot:----
#colors

genus_colors <- c("Bryopsida (Class)" = "#FF9800", 
                  "Blasia sp." = "#D1BB9E", 
                  "Bryaceae (Family)" = "#DC0083", 
                  "Tortula sp." = "#7E8EF1", 
                  "Chloroidium sp." = "#D6D46D",
                  "Pragmopora sp." = "#A94438",
                  "Mucor sp."= "#0F67B1",
                  "Poteriospumella sp." = "#03346E",
                  "Xanthonema sp." = "#F075AA",
                  "Tremella sp." = "grey",
                  "Atrichopsis sp."= "#6EACDA",
                  "Limnomonas sp."= "#36BA98",
                  "Rhogostoma sp." = "#FFA38F",
                  "Tetracladium sp." = "#FFC7ED",
                  "Rhogostomidae (Family)" = "#179BAE",
                  "Sanguina sp." = "#e3e656",
                  "Exophiala sp." = "#4C3BCF",
                  "Chromochloris sp."= "#FF204E",
                  "Chaetothyriales (Order)"= "purple",
                  "NA" = "#3B3C36",
                  "Malassezia sp." = "#F1B42F",
                  "Basidiomycota (Phylum)" = "#81613E", 
                  "Sphaeropleales (Order)"="#9FCB98", 
                  "Pottia sp." = "#C9BEFF")

boxplot_on850_mean <- ggplot(subset_on850_mean, aes(x = factor(OTU, levels = top10_on850_mean$OTU), y = Relative_Abundance, color = Blast_taxa)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5)  +
  ggtitle("C) DDs (315-850m)") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = genus_colors, name = "Taxa") +  # Asigna los colores definidos manualmente
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black")) +
  theme(axis.title.x = element_text(angle = 0, size = 12)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 13, color = "black" )) +
  theme(axis.title.y = element_text(size = 14)) +
  labs(x = "ASV", y = "Relative Abundance %") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(face = "bold", size = 16))  # Título en negrita
boxplot_on850_mean <- boxplot_on850_mean +
  geom_point(data = top10_on850_mean, aes(x = factor(OTU), y = mean_abundance), shape = 18, size = 3, color = "red")
boxplot_on850_mean
#Save as png
ggsave("Relative Abudance_ranked_median and mean/Rel_Ab_ESV/On850_ESV18S.png", boxplot_on850_mean, width = 6, height = 4, dpi = 500)

###Off Stripe - 300-850m-------

#Selecting samples that just belong to 0-250 stripe
subset_off850 <- esv_abundance[esv_abundance$sample_type == "315-850 non-DDs", ]
subset_off850 <- as.data.frame(subset_off850)
#Now, removing the ESVs which doesn't belong to this sample type:
esv_unique_off850 <- subset_off850 %>%
  group_by(OTU) %>%
  filter(sum(Abundance) > 0) %>%
  ungroup()
#Relative abundance per sample:
rel_ab_off850 <- esv_unique_off850 %>%
  group_by(Sample) %>%
  mutate(Relative_Abundance = (Abundance / sum(Abundance)) * 100)
head(rel_ab_off850)
#Remove other columns we don't care:
rel_ab_off850 <- rel_ab_off850 %>%
  dplyr::select(Relative_Abundance, OTU, Sample, Genus)
head(rel_ab_off850)
#In order to probe the rel abun, let's take a sample, sum up all the ESVs and the total should be 100:
sample_SUN0113 <- subset(rel_ab_off850, Sample == "SUN0113")
total_abundance_SUN0113 <- sum(sample_SUN0113$Relative_Abundance)
total_abundance_SUN0113 #it is 100%! it is ok!

#deleting NA 
#(I have noticed that I have some NA in my rela abundance, 
#but it should be a number, so I removed it)
#The explanation should be that those ESVs doesn't belong to that sample group.
rel_ab_off850 <- rel_ab_off850 %>%
  filter(!is.na(Relative_Abundance))

############################## OUTLIERS######
##### Remove outliers based in boxplots (It is optional, 
#sometimes it is better to keep the outliers) #Better decide based on a specific analysis.
## Start here
cleaned_dataframe <- rel_ab_off850 %>%
  group_by(OTU) %>%
  mutate(
    lower_threshold = quantile(Relative_Abundance, 0.25) - 1.5 * IQR(Relative_Abundance, na.rm = TRUE),
    upper_threshold = quantile(Relative_Abundance, 0.75) + 1.5 * IQR(Relative_Abundance, na.rm = TRUE)
  ) %>%
  filter(Relative_Abundance >= lower_threshold & Relative_Abundance <= upper_threshold) %>%
  ungroup()
rel_ab_off850 <- cleaned_dataframe
### end here
##############################end of outliers

###Order by mean:--------------------------------
summary_off850_mean <- aggregate(Relative_Abundance ~ OTU, data = rel_ab_off850, FUN = mean) %>%
  arrange(desc(Relative_Abundance))
#rename column
colnames(summary_off850_mean) <- c("OTU", "mean_abundance")
total_meanabundance <- sum(summary_off850_mean$mean_abundance) #It is 100

#getting the top 10
top10_off850_mean <- summary_off850_mean %>%
  arrange(desc(mean_abundance)) %>%
  head(10)

#Let's rename our data where all the poinst or observations per sample are present
#All data points belong to the top we referred prevously.
subset_off850_mean <- rel_ab_off850 %>%
  filter(OTU %in% top10_off850_mean$OTU)
subset_off850_mean

blast_vec <- c(
  "ASV_2" = "Atrichopsis sp.",
  "ASV_1"  = "Bryopsida (Class)",
  "ASV_5" = "Exophiala sp.",
  "ASV_10" = "Pragmopora sp.",
  "ASV_4" = "Atrichopsis sp.",
  "ASV_7"= "Blasia sp.",
  "ASV_3" = "Bryaceae (Family)",
  "ASV_13" = "Tetracladium sp.",
  "ASV_16"= "Rhogostoma sp.",
  "ASV_33"= "Mortierella sp.",
  "ASV_8" = "Chaetothyriales (Order)",
  "ASV_25" = "Chaetothyriales (Order)"
)

subset_off850_mean$Blast_taxa <- blast_vec[subset_off850_mean$OTU]

##Box plot:----
#COLOR:

genus_colors <- c("Bryopsida (Class)" = "#FF9800", 
                  "Blasia sp." = "#D1BB9E", 
                  "Bryaceae (Family)" = "#DC0083", 
                  "Tortula sp." = "#7E8EF1", 
                  "Chloroidium sp." = "#D6D46D",
                  "Pragmopora sp." = "#A94438",
                  "Mucor sp."= "#0F67B1",
                  "Poteriospumella sp." = "#03346E",
                  "Xanthonema sp." = "#F075AA",
                  "Tremella sp." = "grey",
                  "Atrichopsis sp."= "#6EACDA",
                  "Limnomonas sp."= "#36BA98",
                  "Rhogostoma sp." = "#FFA38F",
                  "Tetracladium sp." = "#FFC7ED",
                  "Rhogostomidae (Family)" = "#179BAE",
                  "Sanguina sp." = "#e3e656",
                  "Exophiala sp." = "#4C3BCF",
                  "Chromochloris sp."= "#FF204E",
                  "Chaetothyriales (Order)"= "purple",
                  "NA" = "#3B3C36",
                  "Malassezia sp." = "#F1B42F",
                  "Basidiomycota (Phylum)" = "#81613E", 
                  "Sphaeropleales (Order)"="#9FCB98", 
                  "Pottia sp." = "#C9BEFF")

boxplot_off850_mean <- ggplot(subset_off850_mean, aes(x = factor(OTU, levels = top10_off850_mean$OTU), y = Relative_Abundance, color = Blast_taxa)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5)  +
  ggtitle("D) non-DDs (315-850m)") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_color_manual(values = genus_colors, name = "Taxa") +  # Asigna los colores definidos manualmente
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11, color = "black")) +
  theme(axis.title.x = element_text(angle = 0, size = 12)) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 13, color = "black")) +
  theme(axis.title.y = element_text(size = 14)) +
  labs(x = "ASV", y = "Relative Abundance %") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(face = "bold", size = 16))  # Título en negrita
boxplot_off850_mean <- boxplot_off850_mean +
  geom_point(data = top10_off850_mean, aes(x = factor(OTU), y = mean_abundance), shape = 18, size = 3, color = "red")
# Mostrar el gráfico
print(boxplot_off850_mean)

#Save as png
ggsave("Relative Abudance_ranked_median and mean/Rel_Ab_ESV/Off850_ESV18S.png", boxplot_off850_mean, width = 6, height = 4, dpi = 500)

#Combining plots:
final_plot <- 
  (boxplot_on250_mean | boxplot_off250_mean) /
  (boxplot_on850_mean | boxplot_off850_mean)
final_plot
ggsave("~/Documents/Steve Lab/Svalbard/Paper figures/RelAb_18S_4GROUPS_Composited.png", 
       final_plot, width = 11, height = 8, dpi = 800)
