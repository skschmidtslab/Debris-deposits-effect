# Supraglacial debris deposits influences the rate and pattern of succession along a High-Arctic glacial foreland. 
## Project description
This study evaluates how debris deposits influence microbial succession along the Midtre Lovénbreen glacier forefield by comparing a debris-free transect with a debris-deposit transect.

## Repository content and data:
1. 16S ASV file: Raw output from the DADA2 pipeline for prokaryotes.
2. 18S ASV file: Raw output from the DADA2 pipeline for eukaryotes.
3. Metadata: Sample and correspondent environmental data.
4. Phyloseq_RelAb_Prok.R: Script for data cleaning and phyloseq object construction using 16S data. It includes workflows for calculating relative abundances, generating bar plots, and identifying core taxa (60%) visualized through Venn diagrams.
5. Phyloseq_RelAb_Euk.R: Script for phyloseq object construction using 18S data. It follows the same processing steps as the prokaryotic workflow for consistency.
6. Diversity_Prok.R: Statistical analysis for prokaryotes, including Alpha diversity (Richness), Beta diversity, PCoA visualizations, and testing via PERMANOVA and PERMDISPER.
7. Diversity_Euk.R: Statistical analysis for eukaryotes, covering Alpha and Beta diversity, PCoA plots, and robustness testing (PERMANOVA/PERMDISPER).
8. SupraglacialDebris_DebrisDeposits_Heatmap.R: Script to visualize total microbial abundances in supraglacial debris and their subsequent projection onto the glacier forefield using heatmaps.
9. Biogeo_Env_factors.R: Analysis of correlations between environmental variables, ASV distributions, and diversity indices. This script also includes comparative analyses of DNA concentration and pH across sample groups.
10. Phototroph_18S.R: Focused analysis of phototrophic communities. It generates pie charts representing relative abundances of phototrophs such as mosses and microalgae (18S).
