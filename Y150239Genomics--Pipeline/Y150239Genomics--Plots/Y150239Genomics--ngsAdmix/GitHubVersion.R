### The BEGINNING ~~~~~
##
# ~ Plots Y150239Genomics--evalAdmix | Written by George Pacheco with help from Jose Samaniego.


# Cleans the environment ~ 
rm(list=ls())


# Loads required packages ~
pacman::p_load(tidyverse, reshape2, scales, optparse, plyr, RColorBrewer, extrafont)
source("visFuns.R")


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Imports data while incorporating annotation ~
corres <- list()
annot <- list()
corres_files <- dir(pattern = ".corres")
annot_files <- dir(pattern = ".labels")
for (k in seq_along(annot_files)) {
  annot[[k]] <- read.table(annot_files[k], sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(annot[[k]]) <- c("Sample_ID")
  annot[[k]]$Population <- ifelse(grepl("FR0", annot[[k]]$Sample_ID), "Sales",
                           ifelse(grepl("KAZ", annot[[k]]$Sample_ID), "Chokpak",
                           ifelse(grepl("Lesina", annot[[k]]$Sample_ID), "Lesina",
                           ifelse(grepl("Crotone", annot[[k]]$Sample_ID), "Crotone",
                           ifelse(grepl("Guglionesi", annot[[k]]$Sample_ID), "Guglionesi",
                           ifelse(grepl("PI22NLD0001M", annot[[k]]$Sample_ID), "Y150239",
                           ifelse(grepl("PD22NLD0146F", annot[[k]]$Sample_ID), "Garderen",
                           ifelse(grepl("PD22NLD0147F", annot[[k]]$Sample_ID), "Garderen",
                           ifelse(grepl("PDOM2022NLD0077M", annot[[k]]$Sample_ID), "Meerkerk",
                           ifelse(grepl("PDOM2022NLD0", annot[[k]]$Sample_ID), "Utrecht", "Error"))))))))))
  corres_df <- as.data.frame(read.table(corres_files[k]))
  rownames(corres_df) <- annot[[k]]$Sample_ID
  colnames(corres_df) <- annot[[k]]$Sample_ID
  corres_df$CHRType <- str_extract(corres_files[k], "(Allosome|Autosomes)")
  corres_df$K <- str_extract(corres_files[k], "(K2|K3|K4|K5|K6|K7)")
  corres[[k]] <- corres_df}


# Perform individual computation ~
compute_on_matrices <- function(corres_ind_list, pop_list, ord_list) {
  result <- lapply(1:length(corres_ind_list), function(i) {
    corres_ind <- corres_ind_list[[i]]
    pop <- pop_list[[i]]
    ord <- ord_list[[i]]
    corres_sub <- corres_ind
    
    # Create annotation data frame
    annotation_df <- data.frame(Sample_ID = rownames(corres_sub),
                                CHRType = rep(NA, nrow(corres_sub)),
                                K = rep(NA, nrow(corres_sub)))
    
    # Match population and add annotation columns
    for (j in 1:length(pop)) {
      annotation_df[pop == names(pop[j]), c("CHRType", "K")] <- unique(ord[pop == names(pop[j])])
    }
    
    # Combine annotation and individual matrix
    indiv_matrix <- cbind(annotation_df, corres_sub)
    
    return(indiv_matrix)
  })
  return(result)}


compute_mean_correlation <- function(corres_ind, pop, exclude_cols) {
  N <- length(pop)
  mean_cors <- matrix(NA, ncol = length(unique(pop)), nrow = length(unique(pop)))
  colnames(mean_cors) <- unique(pop)
  rownames(mean_cors) <- unique(pop)
  
  for (i1 in 1:(length(unique(pop)))) {
    for (i2 in 1:(length(unique(pop)))) {
      p1 <- unique(pop)[i1]
      p2 <- unique(pop)[i2]
      
      # Filter out annotation columns if corres_sub has columns
      corres_sub <- corres_ind[which(pop == p1), which(pop == p2)]
      
      # Handle NAs in corres_sub
      corres_sub[is.na(corres_sub)] <- 0
      
      # Check if corres_sub is not empty
      if (!is.null(corres_sub) && nrow(corres_sub) > 0 && ncol(corres_sub) > 0) {
        # Check if there's only one individual in one or both populations
        if (nrow(corres_sub) <= 1 || ncol(corres_sub) <= 1) {
          mean_cors[i1, i2] <- 0
        } else {
          mean_cors[i1, i2] <- mean(corres_sub, na.rm = TRUE)
        }
      } else {
        mean_cors[i1, i2] <- NA
      }
    }
  }
  return(mean_cors)}


# Applies functions ~
main_list <- compute_mean_correlation(corres_ind_list = corres,
                                 pop_list = lapply(annot, function(ann) ann$Population),
                                 ord_list = lapply(corres, orderInds))

exclude_cols <- c(70, 71)

# Run the function
mean_correlation <- compute_mean_correlation(corres, annot$Population, exclude_cols)





IND <- corres_results[[1]][["individual"]]
MEAN <- corres_results[[1]][["mean_correlation"]]


# /////////////////////// Matrices /////////////////////// # 










# Loads the data ~
Corres_2 <- as.matrix(read.table("AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.corres1"))


# Reads the annotation file ~
ids <- read.table("AllSamples_bcftools.raw.vcf.Filtered.Autosomes.NoKinship.NoTreeSparrow.MAFfiltered.Pruned.K2.labels1", stringsAsFactors = FALSE, sep = "\t", header = FALSE)


# Adds column ids names ~
colnames(ids) <- c("Sample_ID")


# Expands ids by adding Population ~
ids$Population <- ifelse(grepl("FR0", ids$Sample_ID), "Sales",
                         ifelse(grepl("KAZ", ids$Sample_ID), "Chokpak",
                                ifelse(grepl("Lesina", ids$Sample_ID), "Lesina",
                                       ifelse(grepl("Crotone", ids$Sample_ID), "Crotone",
                                              ifelse(grepl("Guglionesi", ids$Sample_ID), "Guglionesi",
                                                     ifelse(grepl("PI22NLD0001M", ids$Sample_ID), "Y150239",
                                                            ifelse(grepl("PD22NLD0146F", ids$Sample_ID), "Garderen",
                                                                   ifelse(grepl("PD22NLD0147F", ids$Sample_ID), "Garderen",
                                                                          ifelse(grepl("PDOM2022NLD0077M", ids$Sample_ID), "Meerkerk",
                                                                                 ifelse(grepl("PDOM2022NLD0", ids$Sample_ID), "Utrecht", "Error"))))))))))


row.names(Corres_2) <- ids$Sample_ID
colnames(Corres_2) <- ids$Sample_ID


cor_mat = Corres_2
pop = as.vector(ids[, 2])
ord = NULL
superpop=NULL
title = "Evaluation of 1000G admixture proportions with K = 3" 
min_z = -.1
max_z = .1
cex.main = 1.5
cex.lab = 1.5
cex.legend = 1.5
color_palette=c("#001260", "#EAEDE9", "#601200")
pop_labels = c(TRUE, TRUE)
plot_legend = TRUE
adjlab = .1
rotatelabpop = 0
rotatelabsuperpop = 0
lineswidth = 1
lineswidthsuperpop = 2
adjlabsuperpop = .16
cex.lab.2 = 1.5


N <- dim(cor_mat)[1]


if(is.null(ord)&!is.null(pop)) ord <- order(pop)
if(is.null(ord)&is.null(pop)) ord <- 1:nrow(cor_mat)


if(is.null(pop)){
  pop <- rep(" ", nrow(cor_mat))
  lineswidth <- 0}


pop <- pop[ord]


N_pop <- vapply(unique(pop[ord]), function(x) sum(pop == x), 1)


cor_mat <- cor_mat[ord, ord]


## Set lower part of matrix as population mean correlation
mean_cors <- matrix(ncol=length(unique(pop)), nrow = length(unique(pop)))
colnames(mean_cors) <- unique(pop)
rownames(mean_cors) <- unique(pop)

for(i1 in 1:(length(unique(pop)))){
  for(i2 in 1:(length(unique(pop)))){
    p1 <- unique(pop)[i1]
    p2 <- unique(pop)[i2]
    mean_cors[i1,i2]<- mean(cor_mat[which(pop == p1),
                                    which(pop == p2)][!is.na(cor_mat[which(pop==p1),
                                                                     which(pop==p2)])])}}

for(i1 in 1:(N-1)){
  for(i2 in (i1+1):N){
    cor_mat[i1, i2] <- mean_cors[pop[i2], pop[i1]]}}

z_lims <- c(min_z, max_z)

if(all(is.na(z_lims))) z_lims <- c(-max(abs(cor_mat[!is.na(cor_mat)])),
                                   max(abs(cor_mat[!is.na(cor_mat)])))


if(any(is.na(z_lims))) z_lims <- c(-z_lims[!is.na(z_lims)], z_lims[!is.na(z_lims)])

min_z <- z_lims[1]
max_z <- z_lims[2]

diag(cor_mat) <- 10
nHalf <- 10

# make sure col palette is centered on 0
Min <- min_z
Max <- max_z
Thresh <- 0

## Make vector of colors for values below threshold
rc1 <- colorRampPalette(colors = color_palette[1:2], space = "Lab")(nHalf)

## Make vector of colors for values above threshold
rc2 <- colorRampPalette(colors = color_palette[2:3], space="Lab")(nHalf)
rampcols <- c(rc1, rc2)

rampcols[c(nHalf, nHalf+1)] <- rgb(t(col2rgb(color_palette[2])), maxColorValue=256) 

rb1 <- seq(Min, Thresh, length.out=nHalf+1)
rb2 <- seq(Thresh, Max, length.out=nHalf+1)[-1]
rampbreaks <- c(rb1, rb2)

rlegend <- as.raster(matrix(rampcols, ncol=1)[length(rampcols):1,])
if(plot_legend){
  layout(matrix(1:2, ncol = 2), width = c(4, 1), height = c(1, 1))
  par(mar=c(5, 4, 4, 0),oma = c(1, 4.5, 2, 0))
}else
  par(mar = c(5, 4, 4, 5), oma = c(1, 4.5, 2, 0))
image(t(cor_mat), col = rampcols, breaks = rampbreaks,
      yaxt = "n", xaxt = "n", zlim = c(min_z, max_z), useRaster = TRUE,
      main = title, 
      oldstyle = TRUE, cex.main = cex.main, xpd = NA)
image(ifelse(t(cor_mat>max_z), 1, NA), col = "darkred", add = TRUE)
if(min(cor_mat) < min_z) image(ifelse(t(cor_mat < min_z), 1, NA), col = "darkslateblue", add = TRUE)
image(ifelse(t(cor_mat == 10), 1, NA), col = "black", add = TRUE)


# Puts pop info ~
if(pop_labels[2])
  text(sort(tapply(1:length(pop), pop, mean) / length(pop)), -adjlab, unique(pop), xpd = NA, cex = cex.lab, srt = rotatelabpop)
if(pop_labels[1])
  text(-adjlab,sort(tapply(1:length(pop), pop, mean) / length(pop)), unique(pop), xpd = NA, cex = cex.lab, srt = 90-rotatelabpop)
abline(v = grconvertX(cumsum(sapply(unique(pop), function(x){sum(pop == x)})) / N, "npc", "user"),
       col = 1, lwd = lineswidth, xpd = FALSE)
abline(h = grconvertY(cumsum(sapply(unique(pop), function(x){sum(pop == x)})) / N, "npc", "user"),
       col = 1, lwd = lineswidth, xpd = FALSE)


# Put superpop if not null ~
if(!is.null(superpop)){
  superpop <- superpop[ord]
  if(pop_labels[2])
    text(sort(tapply(1:length(superpop),superpop,mean)/length(superpop)),-adjlabsuperpop,unique(superpop),xpd = NA,cex=cex.lab.2, srt=rotatelabsuperpop, font=2)
  if(pop_labels[1])
    text(-adjlabsuperpop,sort(tapply(1:length(superpop),superpop,mean)/length(superpop)),unique(superpop),xpd = NA, cex=cex.lab.2,srt=90-rotatelabsuperpop,font=2)
  abline(v=grconvertX(cumsum(sapply(unique(superpop),function(x){sum(superpop == x)}))/N,"npc","user"),
         col = 1,lwd = lineswidthsuperpop,xpd=F)
  abline(h=grconvertY(cumsum(sapply(unique(superpop),function(x){sum(superpop == x)}))/N, "npc", "user"),
         col = 1, lwd = lineswidthsuperpop,xpd=F)}

if(plot_legend){
  par(mar = c(5, .5, 4, 2))
  plot(c(0, 1), c(0, 1), type = 'n', axes = F, xlab = '', ylab = '', main = '')    
  
  rasterImage(rlegend, 0, .25, .4, .75)
  text(x = .8, y = c(.25, .5, .75),
       labels = c(-max(abs(min_z),abs(max_z)), 0, max(abs(min_z), abs(max_z))),
       cex=cex.legend, xpd = NA)
  
  if(plot_legend){
    par(mar = c(5, .5, 4, 2))
    plot(c(0, 1), c(0, 1), type = 'n', axes = F, xlab = '', ylab = '', main = '')    
    
    rasterImage(rlegend, 0, 0.25, 0.4,0.75)
    text(x=0.8, y = c(0.25, 0.5, 0.75),
         labels = c(-max(abs(min_z),abs(max_z)), 0, max(abs(min_z),abs(max_z))),
         cex=cex.legend,xpd = NA)}}


orderInds <- function(q=NULL, pop=NULL, popord=NULL){
  # Function to order individuals for admixture and evalAdmix plots. 
  # recommended is to use pop, then if q is given it will order within pop by admixture proporiton. poporder allows to pre-specify order of populations
  # if only q is given will group individuals by main cluster they are assigned
  
  ordpop <- function(x, pop, q){
    idx <- which(pop==x)
    main_k <- which.max(apply(as.matrix(q[idx,]),2,mean))
    ord <- order(q[idx,main_k])
    idx[ord]
  } 
  
  if(!is.null(pop)){
    
    if(is.null(popord)) popord <- unique(pop)
    
    if(!is.null(q)){ 
      
      ord <- unlist(sapply(popord, ordpop, pop=pop, q=q))
      
    } else if (is.null(q)) {
      
      ord <- unlist(sapply(popord, function(x) which(pop==x)))
      
    }
  } else if (is.null(pop)&!is.null(q)) {
    
    # get index of k with max value per individual
    main_k <- apply(q, 1, which.max)
    
    # get max q per indivdiual
    main_q <- q[cbind(1:nrow(q),main_k)]
    
    ord <- order(main_k, main_q)
    
  } else {stop("Need at least an argument to order.")}
  
  return(ord)}


orderK <- function(q, refinds= NULL,refpops = NULL, pop=NULL){
  # Function to order ancestral populations, useful to keep cluster colors in admix plot the same when comparing results across different k values
  # if you give refinds will use maximum Q value of each individual to define clusters
  # if you give refpops (must also give pops) will use maximum mean admixture proportions within inds from pop to define clusters
  # if any refpops or refinds have same cluster as maximum, the admixture plot will look really bad (you will lose a cluster and another will be twice)
  
  k <- ncol(q)
  kord <- integer(0)
  
  if(is.null(refinds)){
    refpops <- refpops[1:k]
    
    for(p in refpops){
      
      kord <- c(kord, which.max(apply(q[pop==p,],2,mean)))}
  } else {
    
    refinds <- refinds[1:k]
    
    for(i in refinds){
      
      kord <- c(kord, which.max(q[i, ]))}}
  
  # if(any(rowSums(q[,kord]!=1))) warning("reordered admixture proportions don't sum to 1, make sure every refind or refpop defines a unique cluster.")
  
  return(kord)}


# Defines colour palette ~
color_palette <- c("#001260", "#EAEDE9", "#601200")


# Make sure col palette is centered on 0 ~
Min <- -0.1
Max <- 0.1
Thresh <- 0
nHalf <- 10

# Defines the breaks ~
rb1 <- seq(Min, Thresh, length.out = nHalf + 1)
rb2 <- seq(Thresh, Max, length.out = nHalf + 1)[-1]
breaks <- c(rb1, rb2)


# Creates final data frame ~
final_data <- list(cor_mat = cor_mat, mean_cors = mean_cors)


# Merge the melted data frames
cor_mat_df <- melt(final_data$cor_mat)
mean_cors_df <- melt(final_data$mean_cors)


cor_mat_df$value[cor_mat_df$value == 10] <- 0
mean_cors_df$value[mean_cors_df$value == "NaN"] <- 0
mean_cors_df$value[mean_cors_df$value == 10] <- 0


cor_mat_df$CorsType <- "Individual"
mean_cors_df$CorsType <- "Population"


cor_mat_df <- cor_mat_df %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         sorted_pair = ifelse(Var1 < Var2, paste0(Var1, Var2), paste0(Var2, Var1))) %>%
  distinct(sorted_pair, .keep_all = TRUE) %>%
  select(-sorted_pair) %>%
  filter(Var1 != Var2)

mean_cors_df <- mean_cors_df %>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         sorted_pair = ifelse(Var1 < Var2, paste0(Var1, Var2), paste0(Var2, Var1))) %>%
  distinct(sorted_pair, .keep_all = TRUE) %>%
  select(-sorted_pair)

fulldf <- rbind(cor_mat_df, mean_cors_df)


# Creates plot (Heatmap) ~
CareBears <-
  ggplot() +
  geom_tile(data = cor_mat_df, aes(x = Var1, y = Var2, fill = value)) +
  #geom_tile(data = mean_cors_df, aes(x = Var1, y = Var2, fill = value)) +
  scale_x_discrete(limits = rev,
                   expand = c(0, 0)) +
  #scale_y_discrete(limits = rev,
  #                 expand = c(0, 0)) +
  scale_fill_gradientn(name = "Correlation",
                       colours = rampcols,
                       breaks = breaks) + 
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"),
        legend.position = "right",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.margin = margin(t = 0, b = 0, r = 15, l = 15),
        legend.box = "vertical",
        legend.box.margin = margin(t = 20, b = 30, r = 0, l = 0),
        axis.title = element_blank(),
        axis.text.x = element_text(color = "#000000", size = 8.25, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(color = "#000000", size = 8.25, face = "bold"),
        axis.ticks = element_line(color = "#000000", linewidth = .3),
        #strip.text = element_text(colour = "#000000", size = 14, face = "bold", family = "Optima"),
        #strip.background = element_rect(colour = "#000000", fill = "#d6d6d6", linewidth = .3),
        axis.line = element_line(colour = "#000000", linewidth = .3)) +
  guides(fill = guide_legend(title = "", title.theme = element_text(size = 16, face = "bold"),
                             label.theme = element_text(size = 15), reverse = TRUE))


# Saves plot (Boxplot) ~
ggsave(CareBears, file = "CareBears.pdf",
       device = cairo_pdf, limitsize = FALSE, scale = 1, width = 16, height = 16, dpi = 600)


plotCorRes(cor_mat = Corres_2, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 2", max_z = .1, min_z = -.1)


plotCorRes(cor_mat = Corres_3, pop = as.vector(ids[, 1]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 3", max_z = .1, min_z = -.1)
plotCorRes(cor_mat = Corres_4, pop = as.vector(ids[, 1]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 4", max_z = .1, min_z = -.1)
plotCorRes(cor_mat = Corres_5, pop = as.vector(ids[, 1]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 5", max_z = .1, min_z = -.1)
plotCorRes(cor_mat = Corres_6, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 6", max_z = .1, min_z = -.1)
plotCorRes(cor_mat = Corres_7, pop = as.vector(ids[, 2]), ord = NULL,
           title = "Evaluation of 1000G admixture proportions with K = 7", max_z = .1, min_z = -.1)


# Create Ind1 based on the first occurrences of elements in Population
ids$Ind <- ave(seq_along(ids$Population), ids$Population,
               FUN = function(x) paste0(ids$Population[x], "_", sprintf("%02d", seq_along(x))))
levels(ids$Ind <- sub("Y150239_01", "Y150239", ids$Ind))


#
##
### The END ~~~~~