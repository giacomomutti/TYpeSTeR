#!/usr/bin/env Rscript

packages <- c('optparse','ggplot2','ComplexHeatmap', 'magrittr', 'dplyr',
              'patchwork', 'vcfR', 'poppr', 'tidyr', 'viridis')
new_packages <- packages[!packages %in% installed.packages()[,'Package']]

if(length(new_packages)>0) {
  print(paste("please install", new_packages))
  quit(status = 1)
}

library(optparse)

option_list = list(
  make_option(c("-r", "--region"), type="character", default=NA,
              help="STR regions", dest = "regions"),
  make_option(c("-t", "--trf_output"), type="character", default=NA,
              help="STR regions", dest = "trf"),
  make_option(c("-i", "--input"), type = "character", default=NA,
              help="vcf obtained by hipstr", dest = "input"),
  make_option(c("-d","--division"), type = "character", default=NA,
              help="taxonomy df: sample\ttaxonomy_clade", dest = "taxonomy"),
  make_option(c("-o", "--outname"), type="character", default=NA,
              help="names of the output files", metavar="character", dest = "outname"),
  make_option(c("-q","--minquality"), type="double", default=0.9,
              help="minimum quality genotypes", dest = "minqual"),
  make_option(c("-f","--minindels"), type="double", default=0.35,
              help="minimum prop of reads with indels in the flanking regions", dest = "minflank"),
  make_option(c("-s","--minstut"), type="double", default=0.35,
              help="minimum prop of reads with stutter artefacts in the flanking regions", dest = "minstut"),
  make_option(c("-m","--missing"), type="double", default=0.5,
              help="minimum prop of missing samples", dest = "missing"),
  make_option(c("--numclusters"), type="integer", default=2,
              help="number of clusters to split samples", dest = "nclust")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(magrittr)
  library(dplyr)
  library(patchwork)
  library(poppr)
  library(vcfR)
})

source("./scripts/functions.R")
theme_set(theme_bw())

# opt <- NULL
# opt$outname <- "Papio"
# opt$taxonomy <- "taxonomy.csv"
# opt$regions <- "NC_044997.1_filtered_STR.bed"
# opt$trf <- "Panu_Y.fa.2.5.7.80.10.80.6.dat"
# opt$input <- "vcf/papio_highcov_STRs_nn.vcf.gz"
# opt$minqual <- 0.9
# opt$minflank <- 0.35
# opt$minstut <- .35
# opt$missing <- .5
# opt$nclust <- 0

if (is.na(opt$outname)){
  outnm <- "Y_STRs"
} else {
  outnm <- opt$outname
}

if (!is.na(opt$regions) && !is.na(opt$trf)){
  regions <- read.table(opt$regions)
  trf <- read_trf(opt$trf, regions)
  plot_regions <- plot_strs_chrom(trf)
  dir.create("plots", showWarnings = FALSE)
  ggsave(paste0("plots/", outnm, "_regions.pdf"), plot_regions, 
         device = "pdf", width = 8, height = 3)
} else {
  cat("regions file not included: skipping regions plots\n")
}

if (!is.na(opt$taxonomy)){
  taxonomy <- read.delim(opt$taxonomy, col.names = c("sample", "taxon"), header = F)
  n_sps <- length(unique(taxonomy$taxon))
  # find palette with n distinct species
  # eventually by checking if user has given two or three columns you could use the third column as color
  if (n_sps <= 13) {
    col_vector = sample(col_vector, n_sps)
  } else {
    col_vector = c(col_vector, sample(colors(), n_sps-13))
  }
  names(col_vector) <- as.character(unique(taxonomy$taxon))
} else {
  cat("the taxonomy file -d is required to plot vcf EDA files\n")
  quit(status = 1)
}

if (!is.na(opt$input)){
# plot eda vcf
  gt_all <- get_gt(opt$input, taxonomy,
                   q_threshold = opt$minqual, min_dst = opt$minstut, min_dflnk = opt$minflank)
  gt_filtered <- filter_gt(opt$input, taxonomy, trf, missing_data = opt$missing, 
                           q_threshold = opt$minqual, min_dst = opt$minstut, 
                           min_dflnk = opt$minflank)
  n_sam <- length(unique(gt_filtered$Indiv))
  n_site <- length(unique(gt_filtered$Key))
  dimplot <- ifelse(n_sam/5.5<5, 5, n_sam/5.5)
  eda_ps <- plot_eda_vcf(gt_all, col_sp = col_vector)
  ggsave(paste0("plots/", outnm, "_eda.pdf"), eda_ps, 
         device = "pdf", width = dimplot, height = dimplot)
# plot heatmap matrix
  str_cluster <- str_dendro(gt_filtered)
  
  hm <- plot_heatmap(gt = gt_filtered, name_mat = "Y-STRs", dendro = str_cluster$dendro, 
               col_sp = col_vector, row_split = opt$nclust)


  pdf(paste0("plots/", outnm, "_heatmap.pdf"), width = n_site/4, height = n_sam/3.5)
  draw(hm, merge_legend = TRUE,  heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()
} else {
  cat("vcf file not included: skipping vcf plots\n")
}

