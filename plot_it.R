mytheme <- theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = 'white', color = "#e1deda"),
    panel.border = element_rect(linetype = "solid", fill = NA),
    axis.line.x = element_line(colour = 'black', size = 0.6),
    axis.line.y = element_line(colour = 'black', size = 0.6),
    axis.ticks = element_line(colour = 'black', size = 0.35),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 13, color = "black", face = "bold"), 
    legend.key.size = unit(0.8, 'cm'),
    axis.title.x = element_text(family = "Times New Roman", size = 12, color = "black"), 
    axis.title.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"), 
    axis.text.x = element_text(family = "Times New Roman", size = 13, angle = 0, color = "black", face = "bold"), # Making x-axis text bold
    axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
    plot.title = element_text(color = "black", size = 12, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )



myshape= c(19,3,1,2,9,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
MG=c("black","firebrick4","olivedrab3", "tomato2" ,"#66a182","#2e4057","#8d96a0","gold", "#0072B2", "#D55E00","#999999","#2e4057","#8d96a0","#0e669b","#00798c","dodgerblue4", "steelblue2","#00AFBB", "#E7B800", "#FC4E07","lightskyblue4","#82cfd0","#b2e0e4","honeydew3","#8d96a3","lavender","#CC6686","lavenderblush2","mistyrose3","#e1deda","darkgoldenrod","burlywood","papayawhip","wheat4","cornsilk3","khaki2","beige","gray60","gray80","gray96","cadetblue4","#82cfd0","#b2e0e4","honeydew2","mintcream","#0e668b","#a3c4dc", "lightskyblue1","aliceblue","blue","red","green", "lavender","lavenderblush2","mistyrose3","#e1deda","darkgoldenrod","burlywood","papayawhip","wheat4","cornsilk3","khaki2","beige","gray60")


plot_it <- function(physeq) {
  # Rarefy and tidy phyloseq object
  physeq <- rarefy_even_depth(physeq, rngseed = 500, 
                              sample.size = max(10, 0.9 * min(sample_sums(physeq))), 
                              replace = TRUE)
  physeq <- tidy_phyloseq(physeq)
  
  # Exclude specific taxa
  physeq <- subset_taxa(physeq, Species != "Tetragenococcus_halophilus" & Species != "Tetragenococcus_sp.")
  
  # Glom and prune taxa
  glom.1 <- tax_glom(physeq, taxrank = "Genus")
  glom.1 <- prune_taxa(taxa_sums(glom.1) > 0, glom.1)
  glom.2 <- tax_glom(glom.1, taxrank = "Genus")
  glom.3 <- prune_taxa(!grepl("uncultured", tax_table(glom.2)[, "Genus"]), glom.2)
  Genus10 <- names(sort(taxa_sums(glom.3), decreasing = TRUE)[1:20])
  gen10 <- prune_taxa(Genus10, glom.3)
  
  # Convert to data frame
  ps <- psmelt(gen10)
  
  # Filter and preprocess data
  dat <- ps %>%
    filter(!is.na(Habitat)) %>%
    filter(!is.na(Ecoregion_III)) %>%
    filter(!is.na(Host_genus)) %>%
    filter(!is.na(Animal_type)) %>%
    filter(Animal_type != "Toad") %>%
    filter(!is.na(Wild.vs.Captive)) %>%
    filter(Genus != " uncultured" & Genus != " Tetragenococcus")
  
  # Write filtered data to a CSV file
  write.csv(dat, file = "datn.csv")
  
  # Plot
  dat %>%
    group_by(Animal_type) %>%
    mutate(Proportion = Abundance / sum(Abundance, na.rm = TRUE)) %>%
    filter(Proportion > 0 & !is.na(Phylum)) %>%
    ggplot(aes(y = Phylum, x = log10(Proportion), fill = Phylum)) +
    geom_density_ridges2(scale = 1, alpha = 0.8, show.legend = FALSE) +
    ggtitle("Distribution of Relative Abundances by Phylum") +
    labs(x = "Log10(Proportion)", y = NULL) +
    mytheme+scale_fill_manual(values = MG)
  
}
