#Creating a PCoA in 3D
#Load libraries:
library(plotly)
library(magrittr)
library(htmlwidgets)

#Eukaryotes (18S)------
#We need ordBC and pseq.rel
physeq_18_F <- prune_samples(sample_sums(physeq_18S) > 0, physeq_18S)
pseq18.rel <- microbiome::transform(physeq_18_F, "compositional")
#Verifying the sums per sample
sample_sums(pseq18.rel)
#Distances for PCoA:
DistBC18 = phyloseq::distance(pseq18.rel, method = "bray")
ordBC <- ordinate(pseq18.rel, method = "PCoA", distance = DistBC18)

#Extract the first 3 axes from your PCoA results
pcoa_3d_coords <- as.data.frame(ordBC$vectors[, 1:3])
colnames(pcoa_3d_coords) <- c("Axis1", "Axis2", "Axis3")
#Adding the metadata
pcoa_3d_coords$sample_type <- sample_data(pseq18.rel)$sample_type
set.seed(123) # to keep the results consistent
pcoa_3d_coords$Axis1 <- jitter(pcoa_3d_coords$Axis1, amount = 0.005)
pcoa_3d_coords$Axis2 <- jitter(pcoa_3d_coords$Axis2, amount = 0.005)
pcoa_3d_coords$Axis3 <- jitter(pcoa_3d_coords$Axis3, amount = 0.005)
# ---------------------------------------------------------

#Extract variance percentages for the labels
variation <- ordBC$values$Relative_eig * 100
pc1 <- round(variation[1], 1)
pc2 <- round(variation[2], 1)
pc3 <- round(variation[3], 1)

#Create the interactive 3D plot with Percentages in Axis labels
pcoa_3d_18S <- plot_ly(
  data = pcoa_3d_coords, 
  x = ~Axis1, y = ~Axis2, z = ~Axis3, 
  color = ~sample_type, 
  colors = c("#ff686b", "#4cc9f0", "#b30000", "#005f99"), 
  type = 'scatter3d', 
  mode = 'markers',
  marker = list(size = 5, line = list(color = 'black', width = 1))
) %>%
  layout(
    title = "3D PCoA-Eukaryotes",
    scene = list(
      xaxis = list(title = paste0("PCoA Axis 1 (", pc1, "%)")),
      yaxis = list(title = paste0("PCoA Axis 2 (", pc2, "%)")),
      zaxis = list(title = paste0("PCoA Axis 3 (", pc3, "%)"))
    ))
pcoa_3d_18S
#Save it as .html file:
saveWidget(pcoa_3d_18S, "~/Documents/Steve Lab/Svalbard/Paper figures/PCoA_3D_Euk.html", selfcontained = TRUE)


#Prokaryotes (16S)------

#We need ordBC and pseq.rel
physeq_16_filtered <- prune_taxa(taxa_sums(physeq) > 0, physeq) #First let's remove the ASVs that have 0 in at least one sample of the phyloseq
physeq <- physeq_16_filtered 
pseq16.rel <- microbiome::transform(physeq, "compositional")
#Verifying the sums per sample
sample_sums(pseq16.rel)
#Distances for PCoA:
DistBC18 = phyloseq::distance(pseq16.rel, method = "bray")
ordBC <- ordinate(pseq16.rel, method = "PCoA", distance = DistBC18)

#Extract the first 3 axes from your PCoA results
pcoa_3d_coords <- as.data.frame(ordBC$vectors[, 1:3])
colnames(pcoa_3d_coords) <- c("Axis1", "Axis2", "Axis3")
#Adding the metadata
pcoa_3d_coords$sample_type <- sample_data(pseq16.rel)$sample_type
set.seed(123) # to keep the results consistent
pcoa_3d_coords$Axis1 <- jitter(pcoa_3d_coords$Axis1, amount = 0.005)
pcoa_3d_coords$Axis2 <- jitter(pcoa_3d_coords$Axis2, amount = 0.005)
pcoa_3d_coords$Axis3 <- jitter(pcoa_3d_coords$Axis3, amount = 0.005)
# ---------------------------------------------------------

#Extract variance percentages for the labels
variation <- ordBC$values$Relative_eig * 100
pc1 <- round(variation[1], 1)
pc2 <- round(variation[2], 1)
pc3 <- round(variation[3], 1)

#Create the interactive 3D plot with Percentages in Axis labels
pcoa_3d_16S <- plot_ly(
  data = pcoa_3d_coords, 
  x = ~Axis1, y = ~Axis2, z = ~Axis3, 
  color = ~sample_type, 
  colors = c("#ff686b", "#4cc9f0", "#b30000", "#005f99"), 
  type = 'scatter3d', 
  mode = 'markers',
  marker = list(size = 5, line = list(color = 'black', width = 1))
) %>%
  layout(
    title = "3D PCoA-Prokaryotes",
    scene = list(
      xaxis = list(title = paste0("PCoA Axis 1 (", pc1, "%)")),
      yaxis = list(title = paste0("PCoA Axis 2 (", pc2, "%)")),
      zaxis = list(title = paste0("PCoA Axis 3 (", pc3, "%)"))
    ))
pcoa_3d_16S

#save it
saveWidget(pcoa_3d_16S, "~/Documents/Steve Lab/Svalbard/Paper figures/PCoA_3D_Prok.html", selfcontained = TRUE)
