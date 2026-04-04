# Supraglacial debris deposits influences the rate and pattern of succession along a High-Arctic glacial foreland. 
## Project description
This study evaluates how debris deposits influence microbial succession along the Midtre Lovénbreen glacier forefield by comparing a debris-free transect with a debris-deposit transect.

## Repository content and data:
1. 16S ASV file: Raw output from the DADA2 pipeline for prokaryotes.
2. 18S ASV file: Raw output from the DADA2 pipeline for eukaryotes.
3. Metadata: Sample and correspondent environmental data.
4. Physeq_diversity_16S.R: Script for phyloseq object construction using 16S ASVs.
5. Physeq_diversity_18S.R: Script for phyloseq object construction using 18S ASVs.
7. Rel_abundances_18S: Scripts for getting the relative abundances using the physeq object created before.
8. Core_microbiome: Scripts to identify core taxa per group of samples for 16S and 18S (with a prevalence of 60%), and to represent them in a venn diagram. 
9. DDs_SD_18S.R: Script to visualize total microbial abundances in supraglacial debris and their subsequent projection onto the glacier forefield using heatmaps.
10. Phototroph_18S.R: Focused analysis of phototrophic communities. It generates pie charts representing relative abundances of phototrophs such as mosses and microalgae (18S).
   
Note: Some plots were improve in color, format (arrows and lines) and resolution in Inkscape.
