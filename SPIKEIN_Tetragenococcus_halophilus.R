# --------------------------------------------------------
# R Script: SpikeIN scale factor calculation for Tetragenococcus_halophilus
# Author: [Mitra]
# Date: [2024 28 02]
# --------------------------------------------------------

# Load necessary libraries
suppressMessages({
  library(ape)
  library(phyloseq)
  library(dplyr)
  library(microbiome)
})

# working directory
setwd("C:/Users/mghotbi/Desktop/R results16s_21/collection/UHM/Cl_Al_16S_Feb24-20240217T150211Z-001/Cl_Al_16S_Feb24")

# Load and tidy
physeq16S <- readRDS("physeq_pure_16S.rds")

set.seed(123)
physeq16S <- tidy_phyloseq(physeq16S)

# Exclude Tetragenococcus_halophilus and save it
physeq16S_NO_Tetragenococcus_halophilus <- physeq16S %>%
  subset_taxa(Species != "Tetragenococcus_halophilus")
saveRDS(physeq16S_NO_Tetragenococcus_halophilus, file = "physeq16S_NO_Tetragenococcus_halophilus.rds")


# Define spiked cell counts
spiked_cells_s_Dekkera_bruxellensis <- 733
spiked_cells_s_Tetragenococcus_halophilus <- 1847


#physeq16S@sam_data$plate.ID

# Exclude irrelevant samples and keep those spiked 
spiked_16S <- subset_samples(physeq16S, plate.ID %in% c("BAD_plate", "AlexRurik_MAD_plate1", "22UHMwf_JDp3_MAD2", "22UHMwf_JDp1"))

#Normalize based on the spiked volume of the samples
summary(otu_table(spiked_16S))  # Summary statistics of ASV counts before normalization

if (class(spiked_16S@otu_table) == "matrix") {
  spiked_16S@otu_table <- spiked_16S@otu_table / 3
} else {
  spiked_16S@otu_table <- as.matrix(spiked_16S@otu_table) / 3
}

summary(otu_table(spiked_16S))  # Summary statistics of ASV counts after normalization

physeq16S_adj_scaled <- spiked_16S

# scaling factor calculation starts from here
spiked_16S_total_reads <- data.frame(Sample = rownames(sample_data(physeq16S_adj_scaled)), Total_Reads = sample_sums(otu_table(physeq16S_adj_scaled)))
write.csv(spiked_16S_total_reads, file = "spiked_16S_total_reads.csv")

spiked_Tetragenococcus_halophilus_reads <- data.frame(Sample = rownames(sample_data(physeq16S_adj_scaled)), Total_Reads = sample_sums(otu_table(physeq16S_adj_scaled)))
write.csv(spiked_Tetragenococcus_halophilus_reads, file = "spiked_Tetragenococcus_halophilus_reads.csv")

denominator_spiked_Tetragenococcus_halophilus_reads <- spiked_Tetragenococcus_halophilus_reads$Total_Reads
scaling_factor_Tetragenococcus_halophilus <- ifelse(denominator_spiked_Tetragenococcus_halophilus_reads != 0, spiked_cells_s_Tetragenococcus_halophilus / denominator_spiked_Tetragenococcus_halophilus_reads, 0)
write.csv(scaling_factor_Tetragenococcus_halophilus, file = "scaling_factor_Tetragenococcus_halophilus.csv", row.names = FALSE)

# Scale the ASV counts by crossing ASVs to the scalling factor and round the result
physeq16S_adj_scaled_count <- round(otu_table(physeq16S_adj_scaled) * scaling_factor_Tetragenococcus_halophilus)
write.csv(physeq16S_adj_scaled_count, file = "physeq16S_adj_scaled_AbsoluteCount.csv")

# Create a new phyloseq object with the calculated absolute abundance
otu <- read.csv("physeq16S_adj_scaled_AbsoluteCount.csv", header = TRUE, sep = ",", row.names = 1)
otumat <- as.matrix(otu)
random_tree = rtree(ntaxa(physeq16S_adj_scaled), rooted=TRUE, tip.label=taxa_names(physeq16S_adj_scaled))

adj_scaled_spiked_16S <- phyloseq(otu_table(otumat, taxa_are_rows = TRUE),
                           tax_table(physeq16S_adj_scaled),
                           phy_tree(random_tree),
                           sample_data(physeq16S_adj_scaled))    

#Alternatively you could save your time by directly crossing spiked-in calling factor to the ASV table of your new phyloseq object

random_tree <- rtree(ntaxa(physeq16S_adj_scaled), rooted = TRUE, tip.label = taxa_names(physeq16S_adj_scaled))
physeq16S_adj_scaled_absolute_abundance <- phyloseq(
  otu_table(otu_table(physeq16S_adj_scaled) * scaling_factor_Tetragenococcus_halophilus, taxa_are_rows = TRUE),
  tax_table(spiked_16S),
  phy_tree(random_tree),
  sample_data(spiked_16S)
)

print("You did it, Congratulations!")
saveRDS(physeq16S_adj_scaled_absolute_abundance,"physeq16S_adj_scaled_absolute_abundance240319.rds")            
