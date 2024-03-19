# --------------------------------------------------------
# R Script: SpikeIN scale factor calculation for Tetragenococcus_halophilus
# Author: [Mitra]
# Date: [2024 28 02]
# --------------------------------------------------------

# Load necessary libraries
suppressMessages(library(phyloseq, dplyr, microbiome,ape))

# working directory
setwd("C:/Users/mghotbi/Desktop/R results16s_21/collection/UHM/Cl_Al_16S_Feb24-20240217T150211Z-001/Cl_Al_16S_Feb24")

# Load and tidy
physeq16S <- readRDS("physeq16S.rds")
set.seed(123)
physeq16S <- tidy_phyloseq(physeq16S)

# Exclude Tetragenococcus_halophilus and save it
physeq16S_NO_Tetragenococcus_halophilus <- physeq16S %>%
  subset_taxa(Species != "Tetragenococcus_halophilus")
saveRDS(physeq16S_NO_Tetragenococcus_halophilus, file = "physeq16S_NO_Tetragenococcus_halophilus.rds")
meta <- sample_data(physeq16S)

# Define spiked cell counts
spiked_cells_s_Dekkera_bruxellensis <- 733
spiked_cells_s_Tetragenococcus_halophilus <- 1847

# Exclude irrelevant samples and spike-in normalization
spiked_16S <- subset_samples(physeq16S, !(sample.id %in% c(
  "UHM746-21478_S117", "UHM747-21477_S106", "UHM748-21467_S170", "UHM748-21487_S129", 
  "UHM749-21479_S128", "UHM759-21466_S159", "UHM759-21486_S118", "UHM775-21485_S107", "UHM776-21482_S161", "UHM777-21484_S183", "UHM779-21468_S181", "UHM779-21488_S140", "UHM782-21480_S139", "UHM810-21472_S138", "UHM811-21471_S127", "UHM813-21481_S150", "UHM818-21469_S105", "UHM818-21489_S151", "UHM819-21473_S149", "UHM820-21470_S116",
  "UHM820-21490_S162", "UHM827-21474_S160", "UHM829-21476_S182", "UHM832-21483_S172", "UHM905-20598_S45"
)))

#Normalize based on the spiked volume of the samples
if (class(spiked_16S@otu_table) == "matrix") {
  spiked_16S@otu_table <- spiked_16S@otu_table / 3
} else {
  spiked_16S@otu_table <- as.matrix(spiked_16S@otu_table) / 3
}
physeq16S_adj_scaled <- spiked_16S

# scaling factor
spiked_16S_total_reads <- data.frame(Sample = rownames(sample_data(physeq16S_adj_scaled)), Total_Reads = sample_sums(otu_table(physeq16S_adj_scaled)))
write.csv(spiked_16S_total_reads, file = "spiked_16S_total_reads.csv")

spiked_Tetragenococcus_halophilus_reads <- data.frame(Sample = rownames(sample_data(physeq16S_adj_scaled)), Total_Reads = sample_sums(otu_table(physeq16S_adj_scaled)))
write.csv(spiked_Tetragenococcus_halophilus_reads, file = "spiked_Tetragenococcus_halophilus_reads.csv")

denominator_spiked_Tetragenococcus_halophilus_reads <- spiked_Tetragenococcus_halophilus_reads$Total_Reads
scaling_factor_Tetragenococcus_halophilus <- ifelse(denominator_spiked_Tetragenococcus_halophilus_reads != 0, spiked_cells_s_Tetragenococcus_halophilus / denominator_spiked_Tetragenococcus_halophilus_reads, 0)
write.csv(scaling_factor_Tetragenococcus_halophilus, file = "scaling_factor_Tetragenococcus_halophilus.csv", row.names = FALSE)

# Scale the ASV counts by the scaling factor and round the result
physeq16S_adj_scaled_count <- round(otu_table(physeq16S_adj_scaled) * scaling_factor_Tetragenococcus_halophilus)
write.csv(physeq16S_adj_scaled_count, file = "physeq16S_adj_scaled_AbsoluteCount.csv")

# Create a new phyloseq object with the calculated absolute counts
otu <- read.csv("physeq16S_adj_scaled_AbsoluteCount.csv", header = TRUE, sep = ",", row.names = 1)
otumat <- as.matrix(otu)
random_tree = rtree(ntaxa(adj_spiked_16S), rooted=TRUE, tip.label=taxa_names(adj_spiked_16S))

adj_scaled_spiked_16S <- phyloseq(otu_table(otumat, taxa_are_rows = TRUE),
                           tax_table(physeq16S_adj_scaled),
                           phy_tree(random_tree),
                           sample_data(physeq16S_adj_scaled))    

#Alternatively you could save time by directly crossing spiked-in calling factor to the ASV table while you are making you new phyloseq object

# Create new phyloseq object with scaled ASV x spiked-in calling factor 
random_tree <- rtree(ntaxa(spiked_16S), rooted = TRUE, tip.label = taxa_names(spiked_16S))
physeq16S_adj_scaled <- phyloseq(
  otu_table(otu_table(spiked_16S) * scaling_factor_Tetragenococcus_halophilus, taxa_are_rows = TRUE),
  tax_table(spiked_16S),
  phy_tree(random_tree),
  sample_data(spiked_16S)
)

print("You did it, Congratulations!")

                 
