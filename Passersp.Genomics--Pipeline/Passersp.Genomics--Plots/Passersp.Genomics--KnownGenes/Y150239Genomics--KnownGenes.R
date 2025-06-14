### The BEGINNING ~~~~~
##
# Y150239Genomics--TWISST by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads packages ~
pacman::p_load(tidyverse, scales, reshape2, ggh4x, lemon, patchwork, GenomicRanges, txdbmaker,
               rtracklayer, GOstats, GSEABase, outliers, clusterProfiler,BSgenome.Ggallus.UCSC.galGal4, VariantAnnotation, vcfR)


# Imports the House Sparrow annotation ~
HouseGFF <- import("house_sparrow.gff")
HouseGFF_dff <- as.data.frame(HouseGFF)


# Filter for your two genes
genes_of_interest <- c("IV00_00042043", "IV00_00042044", "IV00_00042045", "IV00_00042115")
gene_coords <- HouseGFF_dff %>%
               filter(type == "gene", ID %in% genes_of_interest) %>%
               dplyr::select(seqnames, start, end, strand, ID)


blast_hit <- data.frame(seqnames = "chr10",
                      start = 16156242,
                      end = 16158934,
                      strand = "-",  
                      ID = "rs14949856_blast_hit")

gene_coords <- bind_rows(gene_coords, blast_hit)


# Loads VCF data ~
VCF_auto <- read.vcfR("../../../../LargeFiles/Passersp.Genomics--TriangularR/AllSamples_bcftools.raw.vcf.Autosomes.TriangularR.Focal.ALL.vcf", verbose = TRUE)


# Further process the data ~
genotypes <- extract.gt(VCF_auto)
fix_info <- getFIX(VCF_auto)


# Example: subset for first gene
gene1 <- gene_coords[1, ]
gene2 <- gene_coords[2, ]
gene3 <- gene_coords[3, ]
gene4 <- gene_coords[4, ]
gene5 <- gene_coords[5, ]


# Filter SNPs that fall within that gene region
snp_in_gene1 <- fix_info %>%
  as.data.frame() %>%
  filter(CHROM == as.character(gene1$seqnames),
         POS >= gene1$start,
         POS <= gene1$end)

snp_in_gene2 <- fix_info %>%
  as.data.frame() %>%
  filter(CHROM == as.character(gene2$seqnames),
         POS >= gene2$start,
         POS <= gene2$end)

snp_in_gene3 <- fix_info %>%
  as.data.frame() %>%
  filter(CHROM == as.character(gene3$seqnames),
         POS >= gene3$start,
         POS <= gene3$end)

snp_in_gene4 <- fix_info %>%
  as.data.frame() %>%
  filter(CHROM == as.character(gene4$seqnames),
         POS >= gene4$start,
         POS <= gene4$end)

snp_in_gene5 <- fix_info %>%
                as.data.frame() %>%
                filter(CHROM == as.character(gene5$seqnames),
                POS >= gene5$start,
                POS <= gene5$end)


snp_pos_gene1 <- snp_in_gene1$POS
rows_gene1 <- match(snp_pos_gene1, as.numeric(fix_info[, "POS"]))

snp_pos_gene2 <- snp_in_gene2$POS
rows_gene2 <- match(snp_pos_gene2, as.numeric(fix_info[, "POS"]))

snp_pos_gene3 <- snp_in_gene3$POS
rows_gene3 <- match(snp_pos_gene3, as.numeric(fix_info[, "POS"]))

snp_pos_gene4 <- snp_in_gene4$POS
rows_gene4 <- match(snp_pos_gene4, as.numeric(fix_info[, "POS"]))

snp_pos_gene5 <- snp_in_gene5$POS
rows_gene5 <- match(snp_pos_gene5, as.numeric(fix_info[, "POS"]))


# Sets individual ~
individual <- "PI22NLD0001M_SAMPLE"


individual_geno_gene1 <- genotypes[rows_gene1, individual]
individual_geno_gene2 <- genotypes[rows_gene2, individual]
individual_geno_gene3 <- genotypes[rows_gene3, individual]
individual_geno_gene4 <- genotypes[rows_gene4, individual]
individual_geno_gene5 <- genotypes[rows_gene5, individual]


print(individual_geno_gene5)


geno_df <- cbind(getFIX(VCF_auto), extract.gt(VCF_auto)[, "PI22NLD0001M_SAMPLE"])
colnames(geno_df)[ncol(geno_df)] <- "genotype"
geno_df <- as.data.frame(geno_df)
geno_df$POS <- as.numeric(geno_df$POS)


snp_sub_gene1 <- geno_df %>%
                 filter(CHROM == as.character(gene1$seqnames),
                 POS >= gene1$start,
                 POS <= gene1$end)

snp_sub_gene2 <- geno_df %>%
                 filter(CHROM == as.character(gene2$seqnames),
                 POS >= gene2$start,
                 POS <= gene2$end)

snp_sub_gene3 <- geno_df %>%
                 filter(CHROM == as.character(gene3$seqnames),
                 POS >= gene3$start,
                 POS <= gene3$end)

snp_sub_gene4 <- geno_df %>%
                 filter(CHROM == as.character(gene4$seqnames),
                 POS >= gene4$start,
                 POS <= gene4$end)

snp_sub_gene5 <- geno_df %>%
                 filter(CHROM == as.character(gene5$seqnames),
                 POS >= gene5$start,
                 POS <= gene5$end)

n_snps_gene1 <- nrow(snp_sub_gene1)
n_snps_gene2 <- nrow(snp_sub_gene2)
n_snps_gene3 <- nrow(snp_sub_gene3)
n_snps_gene4 <- nrow(snp_sub_gene4)
n_snps_gene5 <- nrow(snp_sub_gene5)


print(paste("Number of SNPs in", gene1$ID, ":", n_snps_gene1))
print(paste("Number of SNPs in", gene2$ID, ":", n_snps_gene2))
print(paste("Number of SNPs in", gene3$ID, ":", n_snps_gene3))
print(paste("Number of SNPs in", gene4$ID, ":", n_snps_gene4))
print(paste("Number of SNPs in", gene5$ID, ":", n_snps_gene5))


Ggallus <- BSgenome.Ggallus.UCSC.galGal4


# Get 1,000 bp upstream and downstream of SNP position
start_pos <- 15038239 - 500
end_pos   <- 15038239 + 500
query_seq <- as.character(subseq(seq_chr10, start = start_pos, end = end_pos))


# Output as FASTA
cat(">rs14949856_window\n", query_seq, "\n", file = "rs14949856_500bp.fasta")


seq_chr10 <- Ggallus$chr10
base <- subseq(seq_chr10, start = 15038239, end = 15038239)






# Optional: extract gene name
HouseGFF$GeneNames <- stringr::str_extract(HouseGFF$attributes, "(?<=Name=)[^;]+")


# Edits GOTerms lightly ~
GOTerms <- GOTerms |> 
           mutate(GO_Term = GO_Term |> 
           str_replace_all("\\(InterPro\\)", "") |> 
           str_replace_all("\\|", ", "),
           Gene_ID = str_replace_all(Gene_ID, "-.*", ""))


# Creates GOTermsWithData  ~
GOTermsWithData <- GOTerms %>%
                   filter(GO_Term != "-") %>%
                   separate_rows(GO_Term, sep = ", ") %>%
                   mutate(Evidence = "IEA") %>%
                   dplyr::select(GO_Term, Evidence, Gene_ID)


# Creates the GoFrame ~
goFrame <- GOFrame(as.data.frame(GOTermsWithData, organism = "Passerd"))
goAllFrame <- GOAllFrame(goFrame)
GSC <- GeneSetCollection(goAllFrame, setType = GOCollection())


# Sets Gene Universe ~
GenesUniverse <- (unique(GOTerms$Gene_ID))


# Defines categories
categories <- c("upper", "lower", "outliers")


# Initializes lists ~ 
Regions_GR_list <- list()
GeneOverlaps_list <- list()
GenesInRange_df_list <- list()
GenesInRangeOnlyGenes_df_list <- list()
FocalGenes_list <- list()
GOTermsOrtho_Edited_list <- list()
GO_Params_list <- list()
GO_Enrich_list <- list()
GO_Enrich_Top50_list <- list()


for (cat in categories) {df_name <- paste0("filtered_positions_", cat, "_df")
                         Regions_GR_list[[cat]] <- GRanges(seqnames = get(df_name)$CHR,
                                                   ranges = IRanges(start = as.numeric(get(df_name)$Start),
                                                   end = as.numeric(get(df_name)$End)))
  
  # Finds overlaps ~
  GeneOverlaps_list[[cat]] <- findOverlaps(HouseGFF, Regions_GR_list[[cat]])
  
  
  # Extracts genes ~
  GenesInRange <- HouseGFF[queryHits(GeneOverlaps_list[[cat]])]
  GenesInRange_df_list[[cat]] <- data.frame(
                                 Chromosome = as.character(seqnames(GenesInRange)),
                                 Start = start(GenesInRange),
                                 End = end(GenesInRange),
                                 Gene_ID = mcols(GenesInRange)$ID,
                                 Gene_Name = mcols(GenesInRange)$Name,
                                 Type = mcols(GenesInRange)$type)
  GenesInRange_df_list[[cat]]$Gene_ID <- sub("-.*", "", GenesInRange_df_list[[cat]]$Gene_ID)
  
  
  # Filters for genes only ~
  GenesInRangeOnlyGenes_df_list[[cat]] <- GenesInRange_df_list[[cat]] %>%
                                          filter(Type == "gene") %>%
                                          unique()
  
  
  # Extracts focal genes ~
  FocalGenes_list[[cat]] <- as.data.frame(GenesInRangeOnlyGenes_df_list[[cat]]$Gene_ID)
  colnames(FocalGenes_list[[cat]]) <- "Gene_ID"
  
  
  # Saves the lists of Focal Genes ~
  write.table(FocalGenes_list[[cat]], file = paste0("Y150239Genomics--GOAnalysis--FocalGenes_", cat, ".csv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  
  # Edits GOTermsOrtho ~ 
  GOTermsOrtho_Edited_list[[cat]] <- GOTermsOrtho %>%
                                     filter(Gene_ID_Zebra != "" & Gene_ID_House != "") %>%
                                     dplyr::select(Orthogroup, Gene_ID_House, Gene_ID_Zebra) %>%
                                     mutate(Gene_ID_House = str_replace_all(Gene_ID_House, "-.*", "")) %>%
                                     mutate(Gene_ID_House = strsplit(Gene_ID_House, ", ")) %>%
                                     unnest(Gene_ID_House) %>%
                                     filter(str_detect(Gene_ID_House, paste(FocalGenes_list[[cat]]$Gene_ID, collapse = "|"))) %>%
                                     separate_rows(Gene_ID_Zebra, sep = ", ") %>%
                                     dplyr::select(Gene_ID_Zebra)
  
  
  # Saves GOTermsOrtho ~
  write.table(GOTermsOrtho_Edited_list[[cat]], file = paste0("Y150239Genomics--GOAnalysis--GOTermsOrtho_Edited_", cat, ".csv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  
  # Sets GO Analysis parameters ~
GO_Params_list[[cat]] <- GSEAGOHyperGParams(name = paste0("Passerd GO Enrich - ", cat),
                                            geneSetCollection = GSC,
                                            geneIds = FocalGenes_list[[cat]]$Gene_ID,
                                            universeGeneIds = GenesUniverse,
                                            ontology = "BP",
                                            pvalueCutoff = .05,
                                            conditional = FALSE,
                                            testDirection = "over")
  
  
  # Runs GO analysis ~
  Over <- hyperGTest(GO_Params_list[[cat]])
  
  
  # Stores GO enrichment results ~
  GO_Enrich_list[[cat]] <- as.data.frame(summary(Over))
  
  
  # Saves full GO enrichment table ~
  GO_Enrich_list[[cat]] %>% arrange(Pvalue) %>%
                            write.csv(file = paste0("Y150239Genomics--GOAnalysis_", cat, ".csv"))
  
  
  # Get top 50 enriched terms
  GO_Enrich_Top50_list[[cat]] <- GO_Enrich_list[[cat]] %>%
                                 arrange(Pvalue) %>%
                                 head(50)
  
  
  # Defines capitalization function ~
  capitalise_words <- function(text) {
    words <- str_split(text, " ")[[1]]
    exclude_patterns <- c("of", "to", "and", "the", "in")
    patterns_map <- setNames(exclude_patterns, tolower(exclude_patterns))
    process_hyphenated <- function(word) {
      parts <- str_split(word, "-", simplify = TRUE)
      parts <- sapply(seq_along(parts), function(i) {
        part <- parts[i]
        if (i > 1) tolower(part) else str_to_title(part)})
      str_c(parts, collapse = "-")}
    words <- sapply(words, function(word) {
      word_lower <- tolower(word)
      if (word_lower %in% names(patterns_map)) {
        patterns_map[[word_lower]]
      } else if (str_detect(word, "-")) {
        process_hyphenated(word)
      } else {
        str_to_title(word)}})
    str_c(words, collapse = " ")}
  GO_Enrich_Top50_list[[cat]]$Term <- sapply(GO_Enrich_Top50_list[[cat]]$Term, capitalise_words)
  GO_Enrich_Top50_list[[cat]]$Category <- cat}


# Combines GOEnrich_Top50 results ~
GOEnrich_Top50_combined <- do.call(rbind, GO_Enrich_Top50_list)


# Identifies native rows in GOEnrich_Top50_combined ~
original_rows <- GOEnrich_Top50_combined %>%
                 mutate(Count = as.character(Count), 
                 Size = as.character(Size),
                 is_expanded = FALSE)


# Gets expands rows ~
expanded_rows <- GOEnrich_Top50_combined %>%
                 distinct(Term) %>%
                 crossing(Category = unique(GOEnrich_Top50_combined$Category)) %>%
                 anti_join(original_rows, by = c("Term", "Category")) %>%
                 mutate(Pvalue = NA, Count = "", Size = "", is_expanded = TRUE)


# Combines native and expanded rows ~
expanded_df <- bind_rows(original_rows, expanded_rows) %>%
               select(-is_expanded)


# Ensures column order matches original ~
expanded_df <- expanded_df %>% 
               select(names(GOEnrich_Top50_combined)) %>%
               arrange(desc(Term))


# Reorders Category ~
expanded_df$Category <- factor(expanded_df$Category, ordered = TRUE,
                               levels = c("lower",
                                          "upper",
                                          "outliers"))

# Saves file ~
write.csv(expanded_df, file = "Y150239Genomics--GOAnalysis_Top50_Combined.csv", row.names = FALSE)


# Defines y-strip facet labels ~
y_strip_labels <- c("outliers" = "All Outliers",
                    "upper" = "Upper Fence",
                    "lower" = "Lower Fence")


# Creates GoAnalysis plot ~
GOAnalysis_Plot <-
ggplot(expanded_df, aes(x = 1, y = Term)) +
  geom_tile(aes(fill = -log(Pvalue), width = 4), colour = "#000000") +
  geom_text(aes(label = ifelse(Count != "" & Size != "", paste(Count, "/", Size), "")), color = "#000000", size = 3) +
  coord_fixed() +
  scale_fill_gradient(low = "#e5f5f9", high = "#238b45", na.value = "#FAFAFA",) +
  facet_nested(. ~ Category, labeller = labeller(Category = y_strip_labels),
               strip = strip_nested(text_x = elem_list_text(size = 12, family = "Optima", face = "bold", angle = 90),
                                    background_x = elem_list_rect(fill = "#FAFAFA", colour = "#000000", linewidth = .2),
                                    by_layer_x = TRUE)) +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete(expand = c(0, 0)) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(.25, "lines"),
        legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Optima", color = "#000000", size = 16, face = "bold"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(family = "Optima", color = "#000000", size = 8.5, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "#000000", linewidth = .2)) +
  guides(fill = guide_colourbar(title = "Pvalue (-log)", title.theme = element_text(size = 12, family =  "Optima", face = "bold"),
                                label.theme = element_text(size = 10, family =  "Optima", face = "bold"), label.position = "right",
                                barwidth = 1.25, barheight = 12, order = 1, frame.linetype = 1, frame.colour = NA,
                                ticks.colour = "#000000", direction = "vertical", even.steps = TRUE,
                                draw.ulim = TRUE, draw.llim = TRUE))


# Saves plot ~
ggsave(GOAnalysis_Plot, file = "Y150239Genomics--TWISST_GOAnalysis_NEW.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 12, height = 15, scale = 1, dpi = 600)
ggsave(GOAnalysis_Plot, file = "Y150239Genomics--TWISST_GOAnalysis.jpeg",
       limitsize = FALSE, width = 12, height = 15, scale = 1, dpi = 600)


###################################################################################################################################################################################
###################################################################################################################################################################################


# Creates plot ~
Weights_Area_Real_Plot <-
  ggplot() +
  geom_area(data = fulldfUp, aes(x = as.numeric(Mid), y = as.numeric(Value), fill = Weight, group = Weight),
            position = "identity", colour = "#000000", alpha = .3, linetype = 1, linewidth = .2) +
  facet_grid(CHR ~., scales = "free_x", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     #limits = c(0, 147697000),
                     expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                     #breaks = c(50000, 100000, 150000, 200000), 
                     #labels = c("50K", "100K", "150K", "200K"),
                     #limits = c(0, 201000),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Area_Real_Plot, file = "Y150239Genomics--TWISST_AreaReal.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 40, height = 30, scale = 1, dpi = 600)
#ggsave(Weights_Area_Real_Plot, file = "Meerkerkgenomics--TWISST_AreaReal.jpeg",
#       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 600)

# Creates plot ~
Weights_Area_Real_Smooth_Plot <-
  fulldfUp %>%
  ggplot(aes(x = as.numeric(Mid), y = as.numeric(Value), fill = Weight, group = Weight)) +
  stat_smooth(geom = 'area', method = 'loess', position = "identity", linetype = 1, linewidth = .2, colour = "#000000", span = .15, alpha = .3) +
  facet_grid(CHR ~., scales = "free_y", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     limits = c(0, 147697000),
                     expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                     #breaks = c(50000, 100000, 150000, 200000), 
                     #labels = c("50K", "100K", "150K", "200K"),
                     #limits = c(0, 201000),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Area_Real_Smooth_Plot, file = "Meerkerkgenomics--TWISST_WithMeerkerk_AreaRealSmooth.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 40, height = 30, scale = 1, dpi = 600)
#ggsave(Weights_Area_Real_Smooth_Plot, file = "Meerkerkgenomics--TWISST_AreaRealSmooth.jpeg",
#       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 600)


# Creates plot ~
Weights_Line_Plot <-
  ggplot() +
  geom_line(data = fulldfUp, aes(x = as.numeric(Mid), y = as.numeric(Value), colour = Weight, group = Weight),
            position = "identity", linetype = 1, linewidth = .2) +
  facet_rep_grid(CHR ~., scales = "free_y", labeller = labeller(CHR = y_strip_labels)) +
  scale_x_continuous("Genomic Position",
                     breaks = c(25000000, 50000000, 75000000, 100000000, 125000000), 
                     labels = c("25Mb", "50Mb", "75Mb", "100Mb", "125Mb"),
                     limits = c(0, 147697000),
                     expand = c(0, 0)) +
  scale_y_continuous("Weights Across Chrmosomes",
                     #breaks = c(50000, 100000, 150000, 200000), 
                     #labels = c("50K", "100K", "150K", "200K"),
                     #limits = c(0, 201000),
                     expand = c(0, 0)) +
  scale_colour_manual(values = c("#c2a5cf", "#3288bd", "#d53e4f")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_line(color = "#ededed", linetype = "dashed", linewidth = .00005),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "#000000", linewidth = .3),
        axis.title.x = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 30, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, face = "bold", color = "#000000", margin = margin(t = 0, r = 30, b = 0, l = 0)),
        axis.text = element_text(colour = "#000000", size = 8, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        strip.text = element_text(colour = "#000000", size = 11.5, face = "bold"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 30, b = 40, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background = element_blank()) +
  guides(colour = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                               label.theme = element_text(size = 19), override.aes = list(linewidth = 1.75, linetype = 1)),
         fill = guide_legend(title = "Topologies", title.theme = element_text(size = 21, face = "bold"),
                             label.theme = element_text(size = 19), override.aes = list(linewidth = .3, linetype = 1)))


# Saves plot ~
ggsave(Weights_Line_Plot, file = "Meerkerkgenomics--TWISST_WithMeerkerk_LineReal.pdf",
       device = cairo_pdf, limitsize = FALSE, width = 40, height = 30, scale = 1, dpi = 600)
#ggsave(Weights_Line_Plot, file = "Meerkerkgenomics--TWISST_LineReal.jpeg",
#       limitsize = FALSE, width = 30, height = 35, scale = 1, dpi = 600)


#
##
### The END ~~~~~