# --------------------------------------------------------
# R Script: SpikeIN scale factor calculation for Dekkera_bruxellensis
# Author: [Mitra]
# Date: [2024 28 02]
# --------------------------------------------------------

# Load necessary libraries
# Load necessary libraries
suppressMessages(library(phyloseq, dplyr, microbiome))

# Set the working directory
setwd("C:/Users/mghotbi/Desktop/R results16s_21/collection/UHM/Ch_Al 24 Feb R_ITS-20240217T150223Z-001/Ch_Al 24 Feb R_ITS")

# Read metadata
Metadata_ITS <- read.csv("ITS_METADATA_origin.csv", header = TRUE, sep = ",", row.names = 1)

# Filter Spiked_ITS
Spiked_metadata <- Metadata_ITS %>%
  filter(plate_ID %in% c("BAD_plate", "AlexRurik_MAD_plate1", "22UHMwf_JDp3_MAD2", "22UHMwf_JDp1"))


# Import data
physeqITS=readRDS("physeqITS.rds")

# Define spiked cell counts

physeqITS_collapsed <- physeqITS %>%
  tax_glom(taxrank = "Species")%>%
  psmelt()%>%
  as.data.frame()

# Read saved phyloseq object
physeqITS <- readRDS("physeqITS.rds")

# Define spiked cell counts
spiked_cells_s_Dekkera_bruxellensis <- 733
spiked_cells_s_Tetragenococcus_halophilus <- 1847

# Process spiked data
{
#No_spiked <- subset_samples(physeqITS, !(sample.id %in% c( "UHM891-20480_S162", "UHM102-20043_S98", "UHM102-20044_S99",  # Add  sample IDs ..)))
# define spiked data
#spiked_ITS <- subset_samples(physeqITS, !(sample.id %in% sample_names(No_spiked)))
  }
  
#before factor calculation
#collapsing redundant

#s_Dekkera_bruxellensis
Dekkera_bruxellensis <- subset_taxa(physeqITS, Species == "s__Dekkera_bruxellensis") 
Dekkera.only = merge_taxa(Dekkera_bruxellensis, taxa_names(Dekkera_bruxellensis))
Dekkera.only
#Dekkera.mergedSP <- merge_taxa(Dekkera_bruxellensis, taxa_names(Dekkera_bruxellensis)[1:11])

#No_Dekkera_bruxellensis
No_Dekkera_bruxellensis <- subset_taxa(physeqITS, Species != "s__Dekkera_bruxellensis")
MergedITS<-merge_phyloseq(No_Dekkera_bruxellensis,Dekkera.only)

#So this file contains one ASV identified as s__Dekkera_bruxellensis

MergedITS<-spiked_ITS

#factors extra####
{
  Sample_or_Control<-as.factor(sample_data(MergedITS)$Sample_or_Control)
  Plate_ID<-as.factor(sample_data(MergedITS)$Plate_ID)
  Env_broad_scale<-as.factor(sample_data(MergedITS)$Env_broad_scale)
  quant_reading<-as.factor(sample_data(MergedITS)$quant_reading)
  Family<-as.factor(sample_data(MergedITS)$Family)
  host_genus<-as.factor(sample_data(MergedITS)$host_genus)
  host_taxon<-as.factor(sample_data(MergedITS)$host_taxon)
  Metamorphosis<-as.factor(sample_data(MergedITS)$Metamorphosis)
  animal_ecomode<-as.factor(sample_data(MergedITS)$animal_ecomode)
  Diet_Detailed<-as.factor(sample_data(MergedITS)$Diet_Detailed)
  Habitat<-as.factor(sample_data(MergedITS)$Habitat)
  site<-as.factor(sample_data(MergedITS)$site)
  ecomode<-as.factor(sample_data(MergedITS)$ecomode)
  Reproduction<-as.factor(sample_data(MergedITS)$Reproduction)
  ecoregion_IV<-as.factor(sample_data(MergedITS)$ecoregion_IV)
  Clade_Order<-as.factor(sample_data(MergedITS)$Clade_Order)
  animal_type<-as.factor(sample_data(MergedITS)$animal_type)
  ecoregion_III<-as.factor(sample_data(MergedITS)$ecoregion_III)
}


#spiked-> ready for the scaling factor calc

spiked_ITS %>%
  sample_data() %>%
  as.data.frame() %>%
  head()

#In case you have not removed "No-spiked" yet then use this script####
{#define No_spiked and exclude it from the file
# Subset phyloseq object to samples with spike-ins
#spiked_ITS <- subset_samples(MergedITS, !(sample.id %in% sample_names(No_spiked)))
#saveRDS(spiked_ITS,file = "spiked_ITS.rds")
}

# Calculate scaling factors using spiked_ITS file
# Create a data frame with sample names and their corresponding total read counts
spiked_ITS_total_reads <- data.frame(Sample = rownames(sample_data(spiked_ITS)), Total_Reads = sample_sums(otu_table(spiked_ITS)))
spiked_ITS_total_reads

#Dekkera read count
spiked_Dekkera_bruxellensis_reads <- data.frame(Sample = rownames(sample_data(spiked_ITS)), Total_Reads = sample_sums(otu_table(spiked_ITS)))

# Calculate scaling factor for Dekkera_bruxellensis (spiked cells/spiked reads)
denominator_Dekkera_bruxellensis <- spiked_Dekkera_bruxellensis_reads$Total_Reads
scaling_factor_Dekkera_bruxellensis <- ifelse(denominator_Dekkera_bruxellensis != 0, spiked_cells_s_Dekkera_bruxellensis / denominator_Dekkera_bruxellensis, 0)

head(scaling_factor_Dekkera_bruxellensis)
write.csv(scaling_factor_Dekkera_bruxellensis, file = "scaling_factor_Dekkera_bruxellensis.csv", row.names = FALSE)

#now use the scalling factor to calculate the Absolute abundance
# Scale the ASV counts by the scaling factor and round the result
physeqITS_scaled_ASV <- round(otu_table(spiked_ITS) * scaling_factor_Dekkera_bruxellensis)

write.csv(physeqITS_scaled_ASV,file = "Spiked_physeqITS_Absolute.csv")

#Replace the rel abundance you have with the Absolute abundance in you phyloseq obj
#make a new phyloseq object with the Absolute abundance 
library(phyloseq)

# Load OTU table
otu <- read.csv("Spiked_physeqITS_Absolute.csv", header = TRUE, sep = ",",row.names = 1)
otumat <- as.matrix(physeqITS_scaled)
library("ape")

random_tree = rtree(ntaxa(spiked_ITS), rooted=TRUE, tip.label=taxa_names(spiked_ITS))
adj_scaled_spiked_ITS <- phyloseq(otu_table(otumat, taxa_are_rows = TRUE),
                                  tax_table(spiked_ITS),
                                  phy_tree(random_tree),
                                  sample_data(spiked_ITS))

saveRDS(adj_scaled_spiked_ITS,file = "Last_adj_scaled_spiked_ITS.rds")


print("Congratulations, you have your absolute abundance!")


