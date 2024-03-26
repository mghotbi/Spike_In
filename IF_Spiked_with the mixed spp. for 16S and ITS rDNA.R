# --------------------------------------------------------
# R Script: SpikeIN scale factor calculation for Tetragenococcus_halophilus
# Author: [Mitra]
# Date: [2024 26 03]
# --------------------------------------------------------


#all in one for Spike_In
getwd()
# package list:
({
  
reqpkg = c("edgeR", "PoiClaClu","microbiome", "phyloseq","DESeq2", "foreach", "doParallel", 
           "plyr", "reshape2", "ggplot2", "grid","vegan", "scales", "cluster","ape", "dplyr",
           "igraph", "ggnet", "microbiomeutilities", "intergraph", "network", "SpiecEasi",
           "data.table", "decontam", "ggtext", "devtools", "dada2", "ggplot2", "phyloseq",
           "vegan", "ggpubr", "agridat", "lme4", "rstatix", "emmeans", "microbiomeMarker", 
           "lmerTest")
inpkg = installed.packages()[, "Package"]
neededpkg = reqpkg[!reqpkg %in% inpkg]

if (length(neededpkg) > 0) {
  cat("The following package(s) need to be installed:", paste(neededpkg, collapse = ", "), "\n")
  cat("Installing package(s)...\n")
  install.packages(neededpkg)
  cat("Package(s) installed successfully.\n")
} else {
  cat("All required packages are already installed.\n")
}
cat("Loading required packages...\n")
invisible(lapply(reqpkg, library, character.only = TRUE))
cat("Hey there ;) Packages loaded successfully\n :) ")

})

setwd("C:/Users/mghotbi/Desktop/Spike_In factor Calc Workshop/Tetra")
otu<-read.csv("asv2.csv",header = T,row.names = 1,sep = ",")
tax<-read.csv("tax2.csv",header = T,row.names = 1,sep = ",")
meta<-read.csv("Metadata_16s.csv",row.names = NULL,header = T,sep = ",")

#metadata processing
({
# biosample is a common column between animal metadata and sequence metadata
# the sample id should be changed to biosample which is common between animal metadata, 16S and ITS
result <- meta %>%
  rename(biosample = X16s_biosample) %>%  # Renaming X16s_biosample to biosample # we need a common column*
  left_join(seq_meta, by = "biosample") %>% 
  mutate(sample.names = sample.id) %>%  
  rename(sample.id = biosample,  # again renaming biosample to sample.id we do that for ASV and Tax files too
         sample_codes = sample.id)

n <- which(names(result) == "sample_codes")  #where is it
result <- result[,-n]  # Removing the extra sample.id

tibble(result)
write.csv(result, file = "meta.csv", row.names = FALSE) 
meta<-read.csv("meta.csv",row.names = NULL,header = T,sep = ",")

})

# Create phyloseq object with the new metadata consists of spike_volume
cat("prepare the new metadata, including spike volumes :D 👍 ~<:<<<<<> 😊 ")

({
  
otumat <- as.matrix(otu)
taxmat<-as.matrix(tax)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
sampledata = sample_data(data.frame(meta, row.names=sample_names(physeq), stringsAsFactors=FALSE))
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
#MyTree16 <- read.tree("tree16.nwk")
physeq16S = merge_phyloseq(physeq, sampledata)
set.seed(123)
physeq16S <- prune_samples(sample_sums(physeq16S) > 0, physeq16S)
#taxa_names(physeq16S) <- paste0("ASV", seq(ntaxa(physeq16S)))
colnames(tax_table(physeq16S)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
physeq16S <- subset_taxa(physeq16S, Class != "Chloroplast")
physeq16S <- subset_taxa(physeq16S, Order != "Mitochondria")
physeq16S@otu_table@.Data
row.names(otu_table(physeq16S))

})

#clean it up # this part is not mandatory
cat(" ¯\\_(ツ)_/¯ ")

({
  
tidy_phyloseq <- function(my_phyloseq){
  colnames(tax_table(my_phyloseq)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_table(my_phyloseq)[,colnames(tax_table(my_phyloseq))] <- gsub(tax_table(my_phyloseq)[,colnames(tax_table(my_phyloseq))],pattern="[a-z]__",replacement="")
  tax_table(my_phyloseq)[tax_table(my_phyloseq)[,"Phylum"]==NA,"Phylum"] <- "Unidentified"
  return(my_phyloseq)
}

})

physeq16S<-tidy_phyloseq(physeq16S)

#if you like to see the stat before calculating absolute abundance
({

#summary based on ASV_ID 
summary_phyloseq <- function(physeq) {
  otu_table <- otu_table(physeq)
  summary_stats <- data.frame(
    ASV_ID = rownames(otu_table),
    Mean = apply(otu_table, 1, mean),
    Median = apply(otu_table, 1, median),
    SD = apply(otu_table, 1, sd),
    SE = apply(otu_table, 1, function(x) sd(x) / sqrt(length(x)))
  )
  return(summary_stats)
}


#summary based on sample_id
summ_phyloseq <- function(physeq) {
  otu_table <- otu_table(physeq)
  summary_stats <- data.frame(
    Sample_ID = colnames(otu_table),
    Mean = apply(otu_table, 2, mean),
    Median = apply(otu_table, 2, median),
    SD = apply(otu_table, 2, sd),
    SE = apply(otu_table, 2, function(x) sd(x) / sqrt(length(x)))
  )
    return(summary_stats)
}

})

summary(physeq16S@tax_table)
#for each ASV 
summary_phyloseq(physeq16S)
#for each sample
summ_phyloseq(physeq16S)

#keep spiked samples -> 264 samples are spiked # 
spiked_16S <- subset_samples(physeq16S, spiked_volume %in% c("2", "1"))
spiked_16S@tax_table@.Data


#So in case of normalizing 16S/3
#Normalize based on the volume of DNA added to 16S rRNA amplicon PCR #### so each 16S ASV divided to 3
({
  
#(this differs from ITS rDNA (1Micro litter))
if (class(spiked_16S@otu_table) == "matrix") {
  spiked_16S@otu_table <- spiked_16S@otu_table / 3
} else {
  spiked_16S@otu_table <- as.matrix(spiked_16S@otu_table) / 3
}

physeq16S_adj_scaled <- spiked_16S

})

#If you decide to go forward without normalization
#or  non-normalized skip the last step####and keep going with spiked_16S
({
spiked_16S<-spiked_16S

})


#If you like to resample
#or resample x3? from 16S####

({
  
set.seed(20140206)
savedatafile = "simulation"
# Number of ASVs in simulation.  This is the maximum.  For less-diverse
nOTUs = 2000L
# Minimum number of reads
minobs = 1L
# Samples per simulation
J = 40L
# Effect Size. 
mixfacs = c(1, 1.15, 1.25, 1.5, 1.75, 2, 2.5, 3.5)
ns = c(1000, 2000, 5000, 10000)
# replicate numb to repeat for each comb of simulation
reps = 1:5
# delimiter for command 
comdelim = "_"
# Number of cores I guess I dont have more.
Ncores = 8
# rarefying power params
rarefy_fracs = c(0, seq(0.05, 0.25, 0.05), 0.4)
rarereps = 1:2
# The combinations of simulation parameters, used later in the script.
simparams = apply(expand.grid(ns, reps, mixfacs), 1, paste0, collapse = comdelim)
# Define date-stamp for file names
datestamp = gsub(":", "_", gsub("[[:space:]]+", "-", date()), fixed = TRUE)
print(datestamp)

sampsums = sample_sums(spiked_16S)
# Trim ASVs that do not appear in very many samples (prevalence)
samobs = apply(otu_table(spiked_16S), 1, function(x, m) sum(x > m), m = minobs)
otudf = data.frame(prev = samobs, sums = taxa_sums(spiked_16S))
otudf = otudf[order(-otudf$prev, -otudf$sums), ]

template = prune_taxa(rownames(otudf)[1:nOTUs], spiked_16S)
template 

})

({
  
  proportion = function(physeq) {
# Normalize total seq
  normf = function(x, tot = max(sample_sums(physeq))) {
    tot * x/sum(x)  }
  physeq = transform_sample_counts(physeq, normf)
  return(physeq)
}

  
  randomsubsample = function(physeq, smalltrim = 0.15, replace = TRUE) {
    require("phyloseq")
    samplemin = sort(sample_sums(physeq))[-(1:floor(smalltrim * nsamples(physeq)))][1]
    physeqr = rarefy_even_depth(physeq, samplemin, rngseed = FALSE, replace = replace, 
                                trimOTUs = TRUE)
    return(physeqr)
  }
  

  simpletrim = function(physeq, J) {
    if (taxa_are_rows(physeq)) {
      physeq = t(physeq)
    }
     prevalence = apply(as(otu_table(physeq), "matrix"), 2, function(x, minobs) {
      sum(x > minobs)
    }, minobs = 2 * J)/(2 * J)
    # Will only keep OTUs that appear in more than 2% of samples
    keepOTUs = prevalence > 0.05
    # Keep only those OTUs with total reads greater than 0.5x the number of samples
    keepOTUs = keepOTUs & taxa_sums(physeq) > (0.5 * J)
    return(prune_taxa(keepOTUs, physeq))
  }
  

#last purification  
  deseq_varstab = function(physeq, sampleConditions = rep("A", nsamples(physeq)), 
                           ...) {
    require("DESeq2")
    # Enforce orientation.
    if (!taxa_are_rows(physeq)) {
      physeq <- t(physeq)
    }
    x = as(otu_table(physeq), "matrix")
    # The same tweak as for edgeR to avoid NaN problems that cause the workflow
    # to stall/crash.
    x = x + 1
    # Create annotated data.frame with the taxonomy table, in case it is useful
    # later
    taxADF = as(data.frame(as(tax_table(physeq), "matrix"), stringsAsFactors = FALSE), 
                "AnnotatedDataFrame")
    cds = newCountDataSet(x, sampleConditions, featureData = taxADF)
    # First estimate library size factors
    cds = estimateSizeFactors(cds)
    # Variance estimation, passing along additional options
    cds = estimateDispersions(cds, ...)
    # Determine which column(s) have the dispersion estimates
    dispcol = grep("disp\\_", colnames(fData(cds)))
    # Enforce that there are no infinite values in the dispersion estimates
    if (any(!is.finite(fData(cds)[, dispcol]))) {
      fData(cds)[which(!is.finite(fData(cds)[, dispcol])), dispcol] <- 0
    }
    vsmat = exprs(varianceStabilizingTransformation(cds))
    otu_table(physeq) <- otu_table(vsmat, taxa_are_rows = TRUE)
    return(physeq)
  }
  
  
spiked_16S<-simpletrim(spiked_16S, J = 0.01)
spiked_16S_resampled<-randomsubsample(spiked_16S)
#to continue with resampled phylseq use spiked_16S_resampled<-physeq16S_adj_scaled
spiked_16S_resampled<-physeq16S_adj_scaled 

})


# before starting the scaling factor calculation
# we merge 3 Tetragenococcus_halophilus ASVs rooted from spikeIn sp and the rest will stay in the file
# e49935179f23c00fbaf37e529b14aecd
# 2ddb215ff668b6a24a9ccf2bb076a453
# 65ab824f29da7101040d5c7ab82b6441
# 45c10f23ebd7376ca0b353e85bbee992

hashcodes <- c(
  "e49935179f23c00fbaf37e529b14aecd",
  "2ddb215ff668b6a24a9ccf2bb076a453",
  "65ab824f29da7101040d5c7ab82b6441",
  "45c10f23ebd7376ca0b353e85bbee992")

#only Tetra
spiked_Tetra <- subset_taxa(spiked_16S, row.names(tax_table(spiked_16S)) %in% hashcodes)
saveRDS(spiked_Tetra,"spiked_tetra.rds")
#then merge them to one
Tetragenococcus_halophilus.spiked.only = merge_taxa(spiked_Tetra, taxa_names(spiked_Tetra))
saveRDS(Tetragenococcus_halophilus.spiked.only,"Tetragenococcus_halophilus.spiked.only.rds")
#file without  Tetra
No_Tetragenococcus_halophilus <- subset_taxa(spiked_16S, Species != "Tetragenococcus_halophilus" & Species != " Tetragenococcus_sp.")
saveRDS(No_Tetragenococcus_halophilus,"No_Tetragenococcus_halophilus.rds")
#lastly return Tetra as one ASV to the file for the rest of calc
spiked_16S_no_duplicate <-merge_phyloseq(No_Tetragenococcus_halophilus,Tetragenococcus_halophilus.spiked.only)
saveRDS(spiked_16S_no_duplicate,"spiked_16S_no_duplicate.rds")

spiked_16S_no_duplicate@tax_table@.Data

spiked_16S_no_duplicate<-physeq16S_adj_scaled

#file is ready for spike_In calc
#####Here we start with spikeIn factor calculation
#lets keep going with the 16S/3 
# Define spiked cell counts
spiked_cells_s_Dekkera_bruxellensis <- 733
spiked_cells_s_Tetragenococcus_halophilus <- 1847


# Exclude Tetragenococcus_halophilus and save it
physeq16S_NO_Tetragenococcus_halophilus <- physeq16S_adj_scaled %>%
  subset_taxa(Species != "Tetragenococcus_halophilus")
saveRDS(physeq16S_NO_Tetragenococcus_halophilus, file = "physeq16S_NO_Tetragenococcus_halophilus.rds")

# Calculate total reads and Tetragenococcus_halophilus reads
spiked_16S_total_reads <- data.frame(Sample = rownames(sample_data(physeq16S_adj_scaled)), Total_Reads = sample_sums(otu_table(physeq16S_adj_scaled)))
write.csv(spiked_16S_total_reads, file = "spiked_16S_total_reads.csv")

spiked_Tetragenococcus_halophilus_reads <- data.frame(Sample = rownames(sample_data(physeq16S_adj_scaled)), Total_Reads = sample_sums(otu_table(physeq16S_adj_scaled)))
write.csv(spiked_Tetragenococcus_halophilus_reads, file = "spiked_Tetragenococcus_halophilus_reads.csv")

# Calculate scaling factor for Tetragenococcus_halophilus
#so we divide spike cells to number of spiked reads while considering that some plates have a full spike volume (2 µL) and others have a half spike volume (1 µL)
denominator_spiked_Tetragenococcus_halophilus_reads <- spiked_Tetragenococcus_halophilus_reads$Total_Reads
metadata <- sample_data(physeq16S_adj_scaled)
# Calculate scaling factor based on spike volume in metadata (either 2=full or 1=half, 0 is removed)
scaling_factor_Tetragenococcus_halophilus <- ifelse(denominator_spiked_Tetragenococcus_halophilus_reads != 0, {
  ifelse(metadata$spiked_volume == 1, {
    spiked_cells_s_Tetragenococcus_halophilus / 2 / denominator_spiked_Tetragenococcus_halophilus_reads
  }, ifelse(metadata$spiked_volume == 2, {
    spiked_cells_s_Tetragenococcus_halophilus / denominator_spiked_Tetragenococcus_halophilus_reads
  }, 0))
}, 0)

write.csv(scaling_factor_Tetragenococcus_halophilus, file = "scaling_factor_Tetragenococcus_halophilus.csv", row.names = FALSE)
print(scaling_factor_Tetragenococcus_halophilus)
cat("cool!")


# Turn ASV counts to absolute count by crossing ASVs to the scaling factor and round it (funnily we dont have half reads ;) )
physeq16S_adj_scaled_count <- round(otu_table(physeq16S_adj_scaled) * scaling_factor_Tetragenococcus_halophilus)
write.csv(physeq16S_adj_scaled_count, file = "physeq16S_adj_scaled_AbsoluteCount.csv")

# Create a new phyloseq object with the absolute abundance
#otu <- read.csv("physeq16S_adj_scaled_AbsoluteCount.csv", header = TRUE, sep = ",", row.names = 1)
#otumat <- as.matrix(otu)
#random_tree = rtree(ntaxa(physeq16S_adj_scaled), rooted=TRUE, tip.label=taxa_names(physeq16S_adj_scaled))
#adj_scaled_spiked_16S <- phyloseq(otu_table(otumat, taxa_are_rows = TRUE),tax_table(physeq16S_adj_scaled),phy_tree(random_tree),sample_data(physeq16S_adj_scaled))    

# cross the spiked-in scalling factor to the ASV while making a new phyloseq obj 

random_tree <- rtree(ntaxa(physeq16S_adj_scaled), rooted = TRUE, tip.label = taxa_names(physeq16S_adj_scaled))
physeq16S_adj_scaled_absolute_abundance <- phyloseq(
  otu_table(otu_table(physeq16S_adj_scaled) * scaling_factor_Tetragenococcus_halophilus, taxa_are_rows = TRUE),
  tax_table(physeq16S_adj_scaled),
  phy_tree(random_tree),
  sample_data(physeq16S_adj_scaled))


#check the stat after the absolute calc
summary(physeq16S_adj_scaled_absolute_abundance@tax_table)
#for each ASV 
summary_phyloseq(physeq16S_adj_scaled_absolute_abundance)
#for each sample
summ_phyloseq(physeq16S_adj_scaled_absolute_abundance)

sentence <- "🐢🐢🐍🦎🐍🦎🐸🐸🐸 🥳🥳🥳 Congratulations! 🪻🌼🌻🌹🥀🌺🌸Your efforts have produced amazing results, hooray 🥳🥳🥳🍄🍄🍄🍄 👍 👌👌 😊"
for (i in 1:100) {
  cat(sentence, "\n")
}
