### The BEGINNING ~~~~~
##
# Passer sp. Genomics -- KnownGenes | Written by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads packages ~
pacman::p_load(tidyverse, scales, reshape2, ggh4x, lemon, patchwork, GenomicRanges, GenomicFeatures, txdbmaker,
               rtracklayer, GOstats, GSEABase, outliers, clusterProfiler, biomaRt, AnnotationDbi,
               BSgenome.Ggallus.UCSC.galGal4, TxDb.Ggallus.UCSC.galGal4.refGene,
               VariantAnnotation, vcfR, Biostrings)

# Loads VCF data ~
VCF_auto <- read.vcfR("../../../../LargeFiles/Passersp.Genomics--TriangularR/AllSamples_bcftools.raw.vcf.Autosomes.TriangularR.ALL.vcf", verbose = TRUE)
VCF_allo <- read.vcfR("../../../../LargeFiles/Passersp.Genomics--TriangularR/AllSamples_bcftools.raw.vcf.Allosome.TriangularR.ALL.vcf", verbose = TRUE)


vcf_auto_df <- as.data.frame(VCF_auto@fix)
vcf_allo_df <- as.data.frame(VCF_allo@fix)

vcf_auto_df$POS <- as.numeric(vcf_auto_df$POS)
vcf_allo_df$POS <- as.numeric(vcf_allo_df$POS)


# Imports the House Sparrow annotation ~
HouseGFF <- import("house_sparrow.gff")
HouseGFF_dff <- as.data.frame(HouseGFF)


# Filter for your two genes
genes_of_interest_Melissah <- c("IV00_00022782", "IV00_00008837", "IV00_00008839", "IV00_00051423", "IV00_00016654",
                                "IV00_00016726", "IV00_00015896", "IV00_00038225", "IV00_00038223")
gene_coords <- HouseGFF_dff %>%
               filter(type == "gene", ID %in% genes_of_interest_Melissah) %>%
               dplyr::select(seqnames, start, end, strand, ID) %>%
               dplyr::mutate(start = as.numeric(as.character(start)),
                             end = as.numeric(as.character(end)),
                             start = start - 15000,
                             end = end + 15000) %>%
               dplyr::rename(CHR = seqnames, Start = start, End = end, Strand = strand, GeneID = ID) %>%
               dplyr::select(-Strand)


# Loads Zebra Finch cDNA .fasta file ~  
ZebraFinch_cDNA <- "Taeniopygia_guttata.bTaeGut1_v1.p.cdna.all.fa.gz"
ZebraFinch_AllcDNA <- readDNAStringSet(ZebraFinch_cDNA)


# Sets genes of interest ~  
genes_of_interest <- c("ASIP", "OCA2", "MLANA", "PMEL")


# Genes gene information ~
ZebraFinch_GeneList <- lapply(genes_of_interest, function(g) grep(g, names(ZebraFinch_AllcDNA)))


# Combines all gene information ~
ZebraFinch_GeneList <- unique(unlist(ZebraFinch_GeneList))


# Extracts gene sequences ~
ZebraFinch_GeneList <- ZebraFinch_AllcDNA[ZebraFinch_GeneList]


# Corrects headers ~
names(ZebraFinch_GeneList) <- sub(" .*", "", names(ZebraFinch_GeneList))


# Writes gene sequences into a .fasta file ~
writeXStringSet(ZebraFinch_GeneList, filepath = "ZebraFinch_SelectedGenes.fasta")


# Reads Zebra Finch BLAST results ~
ZebraFinch_BLAST <- read.table("ZebraFinch_BLASTresults.txt", sep = "\t", header = FALSE, col.names = c("GeneID", "CHR", "Score", "Length", "Mismatch",
                                                                                                        "GapOpen", "qStart", "qEnd", "Start", "End", "E-value", "Bitscore"))

# Gets ZebraFinch_BestHits ~
ZebraFinch_BestHits <- ZebraFinch_BLAST %>%
                       group_by(GeneID) %>%
                       slice_max(Bitscore, n = 1, with_ties = FALSE) %>%
                       dplyr::select(CHR, Start, End, GeneID) %>%
                       dplyr::mutate(Start = as.numeric(as.character(Start)),
                                     End = as.numeric(as.character(End)),
                                     Start = pmax(Start - 15000, 0),
                                     End = End + 15000) %>%
                       group_by(CHR, Start, End) %>%
                       summarise(GeneID = paste(unique(GeneID), collapse = ", "), .groups = "drop")


# Combines coordinates ~
fulldf <- rbind(gene_coords, ZebraFinch_BestHits)


# Function to extract genotypes from VCF given chromosome, start, and end
extract_genotypes <- function(vcf_obj, chrom, start, end) {
                              region_mask <- getCHROM(vcf_obj) == chrom & getPOS(vcf_obj) >= start & getPOS(vcf_obj) <= end
                              vcf_region <- vcf_obj[region_mask, ]
                              return(extract.gt(vcf_region))}

geno_list <- list()

for (i in seq_len(nrow(fulldf))) {row <- fulldf[i, ]
                                         chrom <- as.character(row$CHR)
                                         start <- row$Start
                                         end <- row$End
                                         gene_id <- row$GeneID
  
  vcf_to_use <- if (chrom %in% getCHROM(VCF_auto)) VCF_auto else VCF_allo
  gt_matrix <- extract_genotypes(vcf_to_use, chrom, start, end)
  geno_list[[gene_id]] <- gt_matrix}


# Flatten the list into a tidy data frame
geno_df <- purrr::imap_dfr(geno_list, function(gt_matrix, gene_id) {
  if (is.null(gt_matrix) || nrow(gt_matrix) == 0) return(NULL)
  gt_df <- as.data.frame(gt_matrix) %>%
    mutate(SNP_ID = rownames(gt_matrix)) %>%
    pivot_longer(-SNP_ID, names_to = "Sample_ID", values_to = "Genotype") %>%
    mutate(Gene_ID = gene_id)
  return(gt_df)})


# Expands PCA_Annot by adding Population ~
geno_df$Population <- ifelse(grepl("PI22NLD0001M_SAMPLE", geno_df$Sample_ID), "Focal Ind.",
                      ifelse(grepl("PD22NLD0146F_SAMPLE", geno_df$Sample_ID), "Garderen_01",
                      ifelse(grepl("PD22NLD0147F_SAMPLE", geno_df$Sample_ID), "Garderen_02",
                      ifelse(grepl("PDOM2022NLD0077M_SAMPLE", geno_df$Sample_ID), "Meerkerk_01",
                      ifelse(grepl("FR0", geno_df$Sample_ID), "Sales",
                      ifelse(grepl("KAZ", geno_df$Sample_ID), "Chokpak",
                      ifelse(grepl("Lesina", geno_df$Sample_ID), "Lesina",
                      ifelse(grepl("Crotone", geno_df$Sample_ID), "Crotone",
                      ifelse(grepl("Guglionesi", geno_df$Sample_ID), "Guglionesi",
                      ifelse(grepl("PDOM2022NLD0", geno_df$Sample_ID), "Utrecht", "Error"))))))))))


# Expands PCA_Annot by adding Species ~
geno_df$Species <- ifelse(geno_df$Population %in% c("Utrecht", "Sales"), "House",
                   ifelse(geno_df$Population %in% c("Chokpak", "Lesina"), "Spanish",
                   ifelse(geno_df$Population %in% c("Crotone", "Guglionesi"), "Italian",
                   ifelse(geno_df$Population %in% c("Garderen_01", "Garderen_02", "Meerkerk_01"), "Control",
                   ifelse(geno_df$Population %in% c("Focal Ind."), "Focal Ind.", "Error")))))


# Finds interesting SNPs ~
interesting_snps <- geno_df %>%
                    group_by(SNP_ID) %>%
                    filter((all(Genotype[Species %in% c("House", "Control")] %in% c(NA, "0/0", "0/1")) &&
                    all(Genotype[Species == "Focal Ind."] == "1/1")) |
                    (all(Genotype[Species %in% c("House", "Control")] == "1/1") &&
                    all(Genotype[Species == "Focal Ind."] %in% c("0/0", "0/1")))) %>%
                    ungroup()


geno_df <- geno_df %>% mutate(Condition = case_when(Gene_ID %in% interesting_snps$Gene_ID & 
                              SNP_ID %in% interesting_snps$SNP_ID &
                               all(Genotype[Species %in% c("House", "Control")] %in% c(NA, "0/0", "0/1")) &
                               all(Genotype[Species == "Focal Ind."] == "1/1") ~ "Cond1",
                              Gene_ID %in% interesting_snps$Gene_ID & 
                              SNP_ID %in% interesting_snps$SNP_ID &
                               all(Genotype[Species %in% c("House", "Control")] == "1/1") &
                               all(Genotype[Species == "Focal Ind."] %in% c("0/0", "0/1")) ~ "Cond2", TRUE ~ NA_character_))


snp_summary <- geno_df %>% group_by(Gene_ID, SNP_ID) %>%
                           summarise(Cond1 = all(Genotype[Species %in% c("House", "Control")] %in% c(NA, "0/0", "0/1")) &&
                                             all(Genotype[Species == "Focal Ind."] == "1/1"),
                                     Cond2 = all(Genotype[Species %in% c("House", "Control")] == "1/1") &&
                                             all(Genotype[Species == "Focal Ind."] %in% c("0/0", "0/1")), .groups = "drop") %>%
                           group_by(Gene_ID) %>%
                           summarise("Number of SNPs" = n(),
                                     "Private Homo Alternative" = sum(Cond1),
                                     "Private Homo Ref. or Hetero" = sum(Cond2), .groups = "drop")


Oi <- fulldf %>%
      left_join(snp_summary, by = c("GeneID" = "Gene_ID")) %>%
      replace_na(list("Number of SNPs" = 0,
                      "Private Homo Alternative" = 0,
                      "Private Homo Ref. or Hetero" = 0))


Oi <- Oi %>%
      mutate(`Region Length` = as.numeric(End) - as.numeric(Start) + 1,
      Coverage = round((`Number of SNPs` / `Region Length`) * 100, 4))


Oi <- Oi %>% 
          dplyr::select(GeneID, CHR, Start, End, `Region Length`, "Number of SNPs", Coverage, "Private Homo Alternative", "Private Homo Ref. or Hetero")


# Saves file ~
write.table(Oi, file = "Passersp.Genomics--KnownGenes.txt", row.names = FALSE, quote = FALSE, sep = "\t")





# Extract full headers
full_headers <- names(ZebraFinch_AllcDNA)[unique(unlist(lapply(genes_of_interest, function(g) grep(g, names(ZebraFinch_AllcDNA)))))]  

# Get GeneID (first word in header)
gene_ids <- sub(" .*", "", full_headers)

# Extract GeneName using regex (look for 'gene_symbol:NAME')
gene_names <- sub(".*gene_symbol:([A-Za-z0-9_-]+).*", "\\1", full_headers)

# Create the named vector: GeneID -> GeneName
gene_key <- as.data.frame(setNames(gene_names, gene_ids))





gene_key <- c("ENSTGUT00000003914.2" = "ASIP",
              "ENSTGUT00000010894.2" = "OCA2",
              "ENSTGUT00000025229.1" = "OCA2",
              "ENSTGUT00000034149.1" = "MLANA",
              "ENSTGUT00000036307.1" = "MLANA",
              "ENSTGUT00000045010.1" = "PMEL",
              "IV00_00022782" = "TYRP1",
              "IV00_00008837" = "HERC2",
              "IV00_00008839" = "HERC2",
              "IV00_00051423" = "MFSD12",
              "IV00_00016654" = "HAND2", 
              "IV00_00016726" = "DCT",
              "IV00_00015896" = "HPGDS",
              "IV00_00038225" = "WNT7A",
              "IV00_00038223" = "WNT7A")
  

Oi <- Oi %>%
      mutate(GeneName = gene_names[GeneID]) %>% 
      relocate(GeneName, .before = 1)


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


# Get 1,000 bp upstream and downstream of SNP position
start_pos <- 15038239 - 500
end_pos   <- 15038239 + 500
query_seq <- as.character(subseq(seq_chr10, start = start_pos, end = end_pos))


# Output as FASTA
cat(">rs14949856_window\n", query_seq, "\n", file = "rs14949856_500bp.fasta")


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