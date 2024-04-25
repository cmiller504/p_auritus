# p_auritus
This repository includes input data for analyses of population genomic structure, quantifying the relative impact of isolation of environment on genomic differentiation, and mapping genomic turnover of a tropical forest frog (doi: 10.3389/fcosc.2024.1366248)

Central African rainforests are predicted to be disproportionately affected by Q6 future climate change. How species will cope with these changes is unclear, but rapid environmental changes will likely impose strong selection pressures. Here we examined environmental drivers of genomic variation in the central African puddle frog (Phrynobatrachus auritus) to identify areas of elevated environmentally-associated turnover. We also compared current and future climate models to pinpoint areas of high genomic vulnerability where allele frequencies will have to shift the most in order to keep pace with future climate change. Neither physical landscape barriers nor the effects of past Pleistocene refugia influenced genomic differentiation. Alternatively, geographic distance and seasonal aspects of precipitation are the most important drivers of SNP allele frequency variation. Patterns of genomic differentiation coincided with key ecological gradients across the forest-savanna ecotone, montane areas, and a coastal to interior rainfall gradient. Areas of greatest vulnerability were found in the lower Sanaga basin, the southeastern region of Cameroon, and southwest Gabon. In contrast with past conservation efforts that have focused on hotspots of species richness or endemism, our findings highlight the importance of maintaining environmentally heterogeneous landscapes to preserve genomic variation and ongoing evolutionary processes in the face of climate change.

Description of the data and file structure

Sample and SNP data:

pauritus_sampleID_coordinates.csv This file contains a sample table with sample ID, site name, and GPS coordinate

pauritus_m3M5n4_maf_rSNP.vcf  This is the filtered VCF file used for downstream analyses

Scripts and input data for analyses:

Generalized Dissimilarity modeling (GDM)

GDM_Pauritus.R script

GDM_input_table_maf_rSNP_fst.csv This is the pairwise matrix of Fst values

bio01.asc, bio04.asc, bio12.asc, bio15.asc, bio19.asc Environmental grids for BioClim variables related to temperature and precipitation

ElevationRivers100.csv, ElevationRiversScaled.csv, ccsm.csv, miroc.csv Resistance distances (CIRCUITSCAPE outputs) for isolation by landscape barriers and Pleistocene refugia

Gradient Forest analyses

GradientForest_Pauritus.R script

frog_frqs.txt minor allele frequencies calculated by population

env_GF.csv values for latitude and longitude and the 5 environmental variables for each sampled site

random_grid_EVs_100000_below7.csv random grids of 100,000 sample points of the environmental variables

Additional: Creating sample bias files for hindcast species distribution model (SDM) of P. auritus using Maxent

create_sampleBias_file.R script

Pauritus_GPS.csv table of GBIF localities

elev.asc ascii file for elevation

ne_10m_land.zip shapefile for land

Sharing/Access information

Sample metadata and sequencing reads are also available on NCBI.

http://...will update when accession numbers are available

