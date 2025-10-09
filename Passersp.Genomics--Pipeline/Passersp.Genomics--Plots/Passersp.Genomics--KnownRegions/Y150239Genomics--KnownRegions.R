### The BEGINNING ~~~~~
##
# Passer sp. Genomics -- KnownRegions | Written by George Pacheco ~


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads packages ~
pacman::p_load(tidyverse, scales, reshape2, ggh4x, lemon, patchwork, GenomicRanges, GenomicFeatures, txdbmaker,
               rtracklayer, GOstats, GSEABase, outliers, clusterProfiler, biomaRt, AnnotationDbi,
               BSgenome.Ggallus.UCSC.galGal4, TxDb.Ggallus.UCSC.galGal4.refGene,
               VariantAnnotation, vcfR, Rsamtools, Biostrings, fuzzyjoin)


# Loads VCF data ~
VCF_auto <- read.vcfR("../../../../LargeFiles/Passersp.Genomics--TriangularR/AllSamples_bcftools.raw.vcf.Autosomes.TriangularR.ALL.vcf", verbose = TRUE)
VCF_allo <- read.vcfR("../../../../LargeFiles/Passersp.Genomics--TriangularR/AllSamples_bcftools.raw.vcf.Allosome.TriangularR.ALL.vcf", verbose = TRUE)


# Gets data frames with the genotypes ~
vcf_auto_df <- as.data.frame(VCF_auto@fix)
vcf_allo_df <- as.data.frame(VCF_allo@fix)
vcf_auto_df$POS <- as.numeric(vcf_auto_df$POS)
vcf_allo_df$POS <- as.numeric(vcf_allo_df$POS)


# Imports the House Sparrow annotation ~
HouseGFF <- import("../../../../LargeFiles/Passersp.Genomics--KnownGenes/house_sparrow.gff")
HouseGFF_dff <- as.data.frame(HouseGFF)


# Imports the Red Junglefowl assembly ~
RedJunglefowlCHR10 <- readDNAStringSet("../../../../LargeFiles/Passersp.Genomics--KnownGenes/Gallus_gallus.Galgal4.dna.chromosome.10.fa.gz")


# Defines region of the trans-eQTL as reported by DOI: 10.1038/s41598-020-57710-7 ~
RedJunglefowl_POS <- 15038239
RedJunglefowl_StartPOS <- RedJunglefowl_POS - 1000
RedJunglefowl_EndPOS   <- RedJunglefowl_POS + 1000


# Extracts the region ~
RedJunglefowl_QuerySeq <- subseq(RedJunglefowlCHR10, start = RedJunglefowl_StartPOS, end = RedJunglefowl_EndPOS)


# Converts to string ~
RedJunglefowl_QuerySeqString <- as.character(RedJunglefowl_QuerySeq)


# Create a new DNAStringSet object ~
RedJunglefowl_Set <- DNAStringSet(RedJunglefowl_QuerySeqString)


# Name it accordingly ~ 
names(RedJunglefowl_Set) <- "RedJunglefowl_Chr10_Pos15038239_±1Kb"


# Loads Zebra Finch cDNA .fasta file ~  
ZebraFinch_cDNA <- "../../../../LargeFiles/Passersp.Genomics--KnownGenes/Taeniopygia_guttata.bTaeGut1_v1.p.cdna.all.fa.gz"
ZebraFinch_AllcDNA <- readDNAStringSet(ZebraFinch_cDNA)


# Sets the pigmentation genes of interest ~
pigmentation_genes_of_interest <- c("IV00_00022782", "IV00_00008837", "IV00_00008839", "IV00_00051423", "IV00_00016654",
                                    "IV00_00016726", "IV00_00015896", "IV00_00038225", "IV00_00038223", "IV00_00042043",
                                    "IV00_00042044", "IV00_00042045", "IV00_00042115")


# Gets coordinates for the pigmentation genes of interest ~
pigmentation_genes_of_interest_coords <- HouseGFF_dff %>%
                                         filter(type == "gene", ID %in% pigmentation_genes_of_interest) %>%
                                         dplyr::select(seqnames, start, end, strand, ID) %>%
                                         dplyr::mutate(start = as.numeric(as.character(start)),
                                                       end = as.numeric(as.character(end)),
                                                       start = start - 15000,
                                                       end = end + 15000) %>%
                                         dplyr::rename(CHR = seqnames, Start = start, End = end, Strand = strand, GeneID = ID) %>%
                                         dplyr::select(-Strand)


# Sets pigmentation genes of interest found in the Zebra Finch ~  
pigmentation_genes_of_interest_ZebraFinch <- c("ASIP", "OCA2", "MLANA", "PMEL")


# Genes gene information ~
ZebraFinch_GeneList <- lapply(pigmentation_genes_of_interest_ZebraFinch, function(x) grep(x, names(ZebraFinch_AllcDNA)))


# Combines all Zebra Finch gene information ~
ZebraFinch_GeneList <- unique(unlist(ZebraFinch_GeneList))


# Extracts gene sequences ~
ZebraFinch_GeneList <- ZebraFinch_AllcDNA[ZebraFinch_GeneList]


# Corrects headers ~
names(ZebraFinch_GeneList) <- sub(" .*", "", names(ZebraFinch_GeneList))


CombinedGeneSet <- c(ZebraFinch_GeneList, RedJunglefowl_Set)


# Writes Zebra Finch gene sequences into a .fasta file ~
writeXStringSet(CombinedGeneSet, filepath = "RegionsOfInterest_BLAST.fasta")


# The BLAST was performed on Saga with the command below ~
# module load BLAST+/2.14.1-gompi-2023a
# blastn -query ZebraFinch_SelectedGenes.fasta -db housesparrow_db -out ZebraFinch_BLASTresults.txt -outfmt 6 -num_threads 6


# Loads Zebra Finch BLAST results ~
ZebraFinch_BLAST <- read.table("RegionsOfInterest_BLASTresults.txt", sep = "\t", header = FALSE,
                               col.names = c("GeneID", "CHR", "Score", "Length", "Mismatch",
                                             "GapOpen", "qStart", "qEnd", "Start", "End", "E-value", "Bitscore"))


# Gets ZebraFinch_BestHits while adding a +- 15K flanking region ~
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


# Combines all gene coordinates ~
fulldf <- rbind(pigmentation_genes_of_interest_coords, ZebraFinch_BestHits)


# Sets function to extract genotypes from VCF ~
extract_genotypes <- function(vcf_obj, chrom, start, end) {region_mask <- getCHROM(vcf_obj) == chrom & getPOS(vcf_obj) >= start & getPOS(vcf_obj) <= end
                                                           vcf_region <- vcf_obj[region_mask, ]
                     return(extract.gt(vcf_region))}


# Runs the extract_genotypes function  ~
geno_list <- list()
for (x in seq_len(nrow(fulldf))) {row <- fulldf[x, ]
                                         chrom <- as.character(row$CHR)
                                         start <- row$Start
                                         end <- row$End
                                         gene_id <- row$GeneID
  
  vcf_to_use <- if (chrom %in% getCHROM(VCF_auto)) VCF_auto else VCF_allo
  gt_matrix <- extract_genotypes(vcf_to_use, chrom, start, end)
  geno_list[[gene_id]] <- gt_matrix}


# Transforms the list into a data frame ~
geno_df <- purrr::imap_dfr(geno_list, function(gt_matrix, gene_id) {if (is.null(gt_matrix) || nrow(gt_matrix) == 0) return(NULL)
                           gt_df <- as.data.frame(gt_matrix) %>%
                           mutate(SNP_ID = rownames(gt_matrix)) %>%
                           pivot_longer(-SNP_ID, names_to = "Sample_ID", values_to = "Genotype") %>%
                           mutate(Gene_ID = gene_id)
                           return(gt_df)})


# Expands geno_df by adding Population ~
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


# Expands geno_df by adding Species ~
geno_df$Species <- ifelse(geno_df$Population %in% c("Utrecht", "Sales"), "House",
                   ifelse(geno_df$Population %in% c("Chokpak", "Lesina"), "Spanish",
                   ifelse(geno_df$Population %in% c("Crotone", "Guglionesi"), "Italian",
                   ifelse(geno_df$Population %in% c("Meerkerk_01"), "Control",
                   ifelse(geno_df$Population %in% c("Focal Ind."), "Focal Ind.", "Error")))))


# Identify condition per SNP and summarise by gene ~
genotype_summary <- geno_df %>%
                    group_by(Gene_ID, SNP_ID) %>%
                    summarise(Cond1 = all(Genotype[Species %in% c("House", "Control")] %in% c(NA, "0/0", "0/1")) &&
                                      all(Genotype[Species == "Focal Ind."] == "1/1"),
                              Cond2 = all(Genotype[Species %in% c("House", "Control")] == "1/1") &&
                                      all(Genotype[Species == "Focal Ind."] %in% c("0/0", "0/1")), .groups = "drop") %>%
                    group_by(Gene_ID) %>%
                    summarise(`Number of SNPs` = n(),
                              `Private Homo Alternative` = sum(Cond1),
                              `Private Homo Ref. or Hetero` = sum(Cond2), .groups = "drop")

# Keep SNP-level summary
#snp_summary <- geno_df %>%
#               group_by(Gene_ID, SNP_ID) %>%
#               summarise(Cond1 = all(Genotype[Species %in% c("House", "Control")] %in% c(NA, "0/0", "0/1")) &&
#                                 all(Genotype[Species == "Focal Ind."] == "1/1"),
#                         Cond2 = all(Genotype[Species %in% c("House", "Control")] == "1/1") &&
#                                 all(Genotype[Species == "Focal Ind."] %in% c("0/0", "0/1")), .groups = "drop")

# Filter to SNPs that obey at least one condition
#snp_summary_hits <- snp_summary %>%
#                    filter(Cond1 | Cond2)


# Merge summary with coordinate table ~
fulldfUp <- fulldf %>%
            left_join(genotype_summary, by = c("GeneID" = "Gene_ID")) %>%
                      replace_na(list(`Number of SNPs` = 0,
                                      `Private Homo Alternative` = 0,
                                      `Private Homo Ref. or Hetero` = 0)) %>%
                      mutate(`Region Length` = as.numeric(End) - as.numeric(Start) + 1,
                             Coverage = round((`Number of SNPs` / `Region Length`) * 100, 4)) %>%
             dplyr::select(GeneID, CHR, Start, End, `Region Length`, "Number of SNPs", Coverage, "Private Homo Alternative", "Private Homo Ref. or Hetero")


# Sets gene key ~
gene_key <- c("ENSTGUT00000003914.2" = "ASIP",
              "ENSTGUT00000010894.2" = "OCA2",
              "ENSTGUT00000025229.1" = "OCA2",
              "ENSTGUT00000034149.1" = "MLANA",
              "ENSTGUT00000036307.1" = "MLANA",
              "ENSTGUT00000045010.1" = "PMEL",
              "IV00_00042043" = "CREBBP",
              "IV00_00042044" = "CREBBP",
              "IV00_00042045" = "CREBBP",
              "IV00_00042115" = "WDR24",
              "RedJunglefowl_Chr10_Pos15038239_±1Kb" = "Trans-eQTL Affecting CREBBP & WDR24",
              "IV00_00022782" = "TYRP1",
              "IV00_00008837" = "HERC2",
              "IV00_00008839" = "HERC2",
              "IV00_00051423" = "MFSD12",
              "IV00_00016654" = "HAND2", 
              "IV00_00016726" = "DCT",
              "IV00_00015896" = "HPGDS",
              "IV00_00038225" = "WNT7A",
              "IV00_00038223" = "WNT7A")


# Expands fulldfUp by adding GeneName ~
fulldfUp <- fulldfUp %>%
            rowwise() %>%
            mutate(GeneName = paste(unique(na.omit(gene_key[trimws(unlist(strsplit(GeneID, ",")))])),
            collapse = ", ")) %>%
            ungroup() %>%
            relocate(GeneName, .before = 1) %>%
            mutate(CHR = as.character(CHR),
                   Start = as.numeric(as.character(Start)),
                   End = as.numeric(as.character(End)))


# Saves table ~
write.table(fulldfUp, file = "Passersp.Genomics--KnownRegions.txt", row.names = FALSE, quote = FALSE, sep = "\t")


# Creats .bed file ~
bed_df <- fulldfUp %>%
          dplyr::select(CHR, Start, End) %>%
          dplyr::mutate(Start = pmax(Start - 1, 0))


# Saves .bed file ~
write.table(bed_df, file = "regions.bed", sep = "\t",
            col.names = FALSE, row.names = FALSE, quote = FALSE)


# The VCF annotation was performed on Saga with the command below ~
#for query in Autosomes Allosome
#do
#java -Xmx4g -jar /cluster/projects/nn10082k/Pacheco/Software/snpEff/snpEff.jar \
#-v Passer -stats AllSamples_bcftools.raw.vcf.${query}.TriangularR_SNPEffSummary.html \
#/cluster/work/users/georgep/TriangularR/AllSamples_bcftools.raw.vcf.${query}.TriangularR.All.vcf \
#> /cluster/work/users/georgep/TriangularR/AllSamples_bcftools.raw.vcf.${query}.TriangularR.All.Annotated.vcf
#done


# Loads .fasta file ~
HS_fasta_path <- "../../../../LargeFiles/Passersp.Genomics--KnownGenes/house_sparrow_genome_assembly-18-11-14_masked.fasta"
HS_fasta <- FaFile(HS_fasta_path)
HS_SeqInfoObj <- seqinfo(HS_fasta)


# Defines function to process the annotated VCFs ~
process_annotated_vcf <- function(vcf_path, seqinfo_obj) {vcf <- readVcf(vcf_path, genome = seqinfo_obj)
                                                          gr <- rowRanges(vcf)
                                                          ann <- info(vcf)$ANN
                                                          df <- as.data.frame(gr)
                                                          df$ANN <- as.list(ann)
                                                          df %>%
                                                          unnest_longer(ANN) %>%
                                                          separate(ANN,
                                                          into = c("Allele", "Annotation", "Impact", "Gene_Name", "Gene_ID",
                                                                  "Feature_Type", "Feature_ID", "Transcript_Biotype", 
                                                                  "Rank", "HGVS.c", "HGVS.p", "cDNA.pos", 
                                                                  "CDS.pos", "AA.pos", "Distance", "ERRORS"),
                                                          sep = "\\|", fill = "right")}

# Applies function ~
VCF_auto_unnested <- process_annotated_vcf("../../../../LargeFiles/Passersp.Genomics--KnownGenes/AllSamples_bcftools.raw.vcf.Autosomes.TriangularR.All.Annotated.Subset.vcf.gz",
                                          HS_SeqInfoObj)
VCF_allo_unnested <- process_annotated_vcf("../../../../LargeFiles/Passersp.Genomics--KnownGenes/AllSamples_bcftools.raw.vcf.Allosome.TriangularR.All.Annotated.Subset.vcf.gz",
                                          HS_SeqInfoObj)


# Combines the two data frames ~
VCF_AllCHRs_unnested <- bind_rows(VCF_auto_unnested, VCF_allo_unnested)


# Filters for high-impact variants ~
high_impact_variants <- VCF_AllCHRs_unnested %>%
                        filter(Impact %in% c("HIGH", "MODERATE")) %>%
                        rename(CHR = seqnames, 
                               Start = start,
                               End = end) %>%
                        mutate(CHR = as.character(CHR),
                               Start = as.numeric(as.character(Start)),
                               End = as.numeric(as.character(End)))


# Incorporates genotype information ~
filtered_variants <- high_impact_variants %>%
                     fuzzy_inner_join(fulldfUp, by = c("CHR" = "CHR", "Start" = "Start", "End" = "End"),
                     match_fun = list(`==`, `>=`, `<=`))


#
##
### The END ~~~~~