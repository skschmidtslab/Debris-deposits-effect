#Core microbiome
#Set working directory first 
#Packages:
library(VennDiagram)

####16S----
table(meta(physeq)$sample_type, useNA = "always") #This shows # of samples in each group.
physeq.rel<- microbiome::transform(physeq, "compositional") #Make the abundances compositional
sample_type <- unique(as.character(meta(physeq.rel)$sample_type)) 
print(sample_type)
list_core <- c() # an empty object to store information
#Function for getting the core microbiome with 60% of prevalence:
for (n in sample_type){ # for each variable n in DiseaseState
  print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(physeq.rel, sample_type == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.60)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  print(list_core)
}

#Ven diagram:
venn.plot16 <- venn.diagram(
  x = list(
    "0-215 non-DDs" = list_core[["0-215 non-DDs"]],
    "0-215 DDs"     = list_core[["0-215 DDs"]],
    "315-850 non-DDs" = list_core[["315-850 non-DDs"]],
    "315-850 DDs"   = list_core[["315-850 DDs"]]
  ),
  filename = NULL,
  fill = c("#4cc9f0", "#ff686b", "#005f99", "#b30000"),  # Colors
  alpha = 0.6,
  cex = 2.5,                 #Size of numbers
  fontface = "bold",
  cat.cex = 2.5,             #Size of letters
  cat.fontface = "bold",
  cat.col = c("#4cc9f0", "#ff686b", "#005f99", "#b30000"),  # Colors 
  cat.pos = c(-17, 8, -1.1, -1.1),  #Supperposition adjustment
  cat.dist = c(0.22, 0.22, 0.1, 0.1),#distance of the letter
  margin = 0.05,
  main.cex = 2,
  main.fontface = "bold",
  main.pos = c(0.5, 1)
)

#Show it!
grid.newpage()
pushViewport(viewport(width = 0.8, height = 0.8))  #Size control (0–1)
grid.draw(venn.plot16)
popViewport()

ggsave("Core Microbiome/Core_Microbiome_16S.png",
       plot = venn.plot16,width = 13, height = 7.3, units = "in", dpi = 700)

####18S----
table(meta(physeq_18S)$sample_type, useNA = "always") #This shows # of samples in each group.
physeq18S.rel<- microbiome::transform(physeq_18S, "compositional") #Make the abundances compositional
sample_type <- unique(as.character(meta(physeq18S.rel)$sample_type)) 
print(sample_type)
list_core <- c() # an empty object to store information
#Function for getting the core microbiome with 60% of prevalence:
for (n in sample_type){ # for each variable n in DiseaseState
  print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(physeq18S.rel, sample_type == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.60)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  print(list_core)
}

#Ven diagram:
venn.plot18 <- venn.diagram(
  x = list(
    "0-215 non-DDs" = list_core[["0-215 non-DDs"]],
    "0-215 DDs"     = list_core[["0-215 DDs"]],
    "315-850 non-DDs" = list_core[["315-850 non-DDs"]],
    "315-850 DDs"   = list_core[["315-850 DDs"]]
  ),
  filename = NULL,
  fill = c("#4cc9f0", "#ff686b", "#005f99", "#b30000"),  # Colors
  alpha = 0.6,
  cex = 2.5,                 #Size of numbers
  fontface = "bold",
  cat.cex = 2.5,             #Size of letters
  cat.fontface = "bold",
  cat.col = c("#4cc9f0", "#ff686b", "#005f99", "#b30000"),  # Colors 
  cat.pos = c(-17, 8, -1.1, -1.1),  #Supperposition adjustment
  cat.dist = c(0.22, 0.22, 0.1, 0.1),#distance of the letter
  margin = 0.05,
  main.cex = 2,
  main.fontface = "bold",
  main.pos = c(0.5, 1)
)

#Show it!
grid.newpage()
pushViewport(viewport(width = 0.8, height = 0.8))  #Size control (0–1)
grid.draw(venn.plot18)
popViewport()

ggsave("Core Microbiome/Core_Microbiome_18S.png",
       plot = venn.plot18,width = 13, height = 7.3, units = "in", dpi = 700)

