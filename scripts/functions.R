# from ggthemes:tableau_color_pal classic cyclic 10
col_vector <- c("#1f83b4", "#12a2a8", "#2ca030", "#78a641",
                "#bcbd22", "#ffbf50", "#ffaa0e", "#ff7f0e",
                "#d63a3a", "#c7519c", "#ba43b4", "#8a60b0", "#6f63bb")


# read TRF data file
read_trf <- function(trf, regions) {
  trf_df <- read.delim(trf, skip = 15, sep = "",
                       col.names = c("START", "END", "size", "ref_gt",
                                     "size_cons", "matches_adj", "indel_adj", "score",
                                     "A", "C", "G", "T", "entropy", "motif", "consensus")) %>%
    dplyr::filter(size>2) %>%
    dplyr::mutate(used=ifelse(END %in% regions$V3, TRUE, FALSE))
  trf_df <- trf_df[!duplicated(trf_df$START),]
  return(trf_df)
}


# plot STRs distribution plot
plot_strs_chrom <- function(trf) {
  n_sz <- length(unique(trf$size))
  col_sizes <- viridis::viridis(n_sz)
  # col_sizes = c( "#8cc8bc", "#7da7ea", "#5773c0", "#1d4497")
  names(col_sizes) <- sort(unique(trf$size))

  bar_period <- trf %>%
    filter(size>2) %>%
    mutate(period_size=as.factor(size)) %>%
    ggplot(aes(period_size, fill=period_size, alpha=used)) +
    geom_bar(position = "dodge") +
    scale_y_continuous(expand = expansion(c(0,0.05))) +
    labs(x="Motif Size") +
    guides(fill="none") +
    scale_fill_manual(values = col_sizes) +
    scale_alpha_manual(values = c(0.3,1)) +
    theme(legend.position = "none")

  position_scatter <- trf %>%
    filter(size>2) %>%
    mutate(period_size=as.factor(size)) %>%
    ggplot(aes((START+END)/2, END-START, color=period_size, alpha=used)) +
    geom_point() +
    # annotate("rect", xmin=0, xmax=par, ymin=0, ymax=150, alpha=0.2, fill="grey") +
    scale_color_manual(values = col_sizes) +
    coord_cartesian(expand = F, ylim = c(0, quantile(trf$END-trf$START, .9))) +
    scale_alpha_manual(values = c(0.2, 1)) +
    labs(x="Position", y="length") +
    scale_x_continuous(breaks = seq(0, 100e6, 1e6),# minor_breaks = seq(0,9e6,2.5e5),
                       labels = scales::unit_format(unit = "Mb", scale = 1e-6)) +
    guides(color = guide_legend(title = "Motif Size")) +
    theme(legend.position = "bottom")

  str_eda <- (bar_period | position_scatter) + plot_layout(widths = c(1,3.5)) + plot_annotation(tag_levels = 'A')

  return(str_eda)
}

# get vcf as dataframe
get_gt <- function(file_nm, meta, q_threshold=.9, min_dst=.35, min_dflnk=.35) {
  vcf <- read.vcfR(file_nm, verbose=FALSE)
  gt <- extract_gt_tidy(vcf, verbose=FALSE)
  gt <- gt %>%
      left_join(meta, by = c("Indiv"="sample")) %>%
      mutate(gt_GB=ifelse(gt_Q<=q_threshold | gt_DSTUTTER/gt_DP > min_dst | gt_DFLANKINDEL/gt_DP > min_dflnk, NA, gt_GB),
             gt_GT_alleles = ifelse(is.na(gt_GB), NA, gt_GT_alleles))

# results <- NULL
  # results$gt <- gt
  # results$info <- info
  return(gt)
}


# obtain gt dataframe with all info
filter_gt <- function(file_nm, meta, trf,
                      missing_data=.5, q_threshold=.9, min_dst=.35, min_dflnk=.35,
                      exclude_samples = c(), exclude=c()) {
  vcf <- read.vcfR(file_nm, verbose=FALSE)
  gt <- extract_gt_tidy(vcf, verbose=FALSE)
  info <- extract_info_tidy(vcf) %>%
    distinct(START, .keep_all = T) %>%
    # left_join(final_names %>% filter(reference==ref)) %>%
    left_join(trf)
  info$size <- as.factor(info$PERIOD)
  min_samples <- missing_data*(dim(vcf)[3]-1)

  gt <- gt %>%
    filter(!Indiv %in% exclude_samples) %>%
    left_join(meta, by = c("Indiv"="sample")) %>%
    mutate(gt_GB=ifelse(gt_Q<=q_threshold | gt_DSTUTTER/gt_DP > min_dst | gt_DFLANKINDEL/gt_DP > min_dflnk, NA, gt_GB),
           gt_GT_alleles = ifelse(is.na(gt_GB), NA, gt_GT_alleles)) %>%
    mutate(length_allele=nchar(gt_GT_alleles)) %>%
    group_by(Key) %>%
    mutate(count_GT=sum(is.na(gt_GB), na.rm = T)) %>%
    filter(count_GT <= min_samples) %>%
    ungroup() %>%
    arrange(Key) %>%
    left_join(info) %>%
    # mutate(length_allele=str_count(gt_GT_alleles, motif)) %>%
    filter(!START %in% exclude)

  # results <- NULL
  # results$gt <- gt
  # results$info <- info
  return(gt)
}


get_matrix <- function(gt, row_names=FALSE) {
  out_matrix <- gt %>%
    mutate(final_name=paste0("Site_", Key)) %>%
    select(final_name, Indiv, taxon, length_allele) %>%
    group_by(final_name) %>%
    tidyr::pivot_wider(names_from = final_name, values_from = length_allele) %>%
    arrange(desc(taxon))

  if (row_names==TRUE){
    out_matrix <- out_matrix %>% tibble::column_to_rownames("Indiv")
  }
  return(out_matrix)
}

get_table_alleles <- function(gt) {
  table <- gt %>%
    select(final_name, Sample, Species, gt_GT_alleles) %>%
    group_by(final_name) %>%
    pivot_wider(names_from = final_name, values_from = gt_GT_alleles) %>%
    arrange(desc(Species))

  vec <- c("Reference", "P. anubis", gt[match(colnames(table)[3:length(colnames(table))], gt$final_name),]$consensus)

  table <- rbind(table, vec)
  return(table)
}


str_dendro <- function(gt){
  out_matrix <- get_matrix(gt, row_names = T)
  # print(out_matrix)
  genind <- out_matrix %>%
    select(-taxon) %>%
    adegenet::df2genind(ploidy = 1L,
                        pop = out_matrix$taxon,
                        strata = select(out_matrix, taxon)) %>%
    poppr::as.genclone()

  poppr_res <- poppr::poppr(genind, strata = ~taxon)
  gt$final_name <- paste0("Site_",gt$Key)
  replens <- pull(gt[match(levels(genind@loc.fac), gt$final_name), "size"])

  dist_matrix <- as.matrix(poppr::bruvo.dist(genind, replen = replens))
  if (length(dist_matrix[is.na(dist_matrix)])>0){
    dist_matrix_miss <- as.matrix(ape::ultrametric(dist_matrix))
    rownames(dist_matrix_miss) <- rownames(dist_matrix)
    colnames(dist_matrix_miss) <- colnames(dist_matrix)
    dist_matrix_miss <- dist_matrix
    dendro_hm <- as.hclust(poppr::upgma(dist_matrix))
  } else{
    dendro_hm <- as.hclust(poppr::upgma(dist_matrix))
  }

  results <- NULL
  results$dist_matrix <- dist_matrix
  results$dendro <- dendro_hm
  results$stats <- poppr_res
  return(results)
}

get_str_dendro <- function(gt){
  out_matrix <- get_matrix(gt, row_names = T)
  genind <- out_matrix %>%
    select(-Species) %>%
    adegenet::df2genind(ploidy = 1L,
                        pop = out_matrix$Species,
                        strata = select(out_matrix, Species)) %>%
    poppr::as.genclone()

  replens <- pull(gt[match(levels(genind@loc.fac), gt$final_name), "size"])

  dist_matrix <- as.matrix(poppr::bruvo.dist(genind, replen = replens))
  if (length(dist_matrix[is.na(dist_matrix)])>0){
    dist_matrix_miss <- as.matrix(ape::ultrametric(dist_matrix))
    rownames(dist_matrix_miss) <- rownames(dist_matrix)
    colnames(dist_matrix_miss) <- colnames(dist_matrix)
    dendro_hm <- poppr::upgma(dist_matrix_miss)
  } else{
    dendro_hm <- poppr::upgma(dist_matrix)
  }
  return(dendro_hm)
}



plot_heatmap <- function(gt, name_mat, dendro, col_sp=col_species, row_split=6) {

  heatmap_matrix <- get_matrix(gt, row_names = T)
  mat <- heatmap_matrix %>% select(-taxon) %>% as.matrix()

  # project of each sample -> bars
  # project_bar <- pull(meta[match(rownames(mat), meta$Sample), "Project"])
  # col_project <- as.character(MetBrewer::met.brewer("Cross", length(unique(project_bar))))
  # names(col_project) <- unique(project_bar)

  gt_red <- gt %>%
    mutate(stutter = gt_DSTUTTER/gt_DP, flanking = gt_DFLANKINDEL/gt_DP) %>%
    mutate(final_name=paste0("Site_", Key)) %>%
    group_by(final_name)

  q_mat <- gt_red %>%
    select(final_name, Indiv, gt_Q) %>%
    tidyr::pivot_wider(names_from = final_name, values_from = gt_Q) %>%
    select(-Indiv) %>%
    as.matrix()

  stut_mat <- gt_red %>%
    select(final_name, Indiv, stutter) %>%
    tidyr::pivot_wider(names_from = final_name, values_from = stutter) %>%
    select(-Indiv) %>%
    as.matrix()

  flank_mat <- gt_red %>%
    select(final_name, Indiv, flanking) %>%
    tidyr::pivot_wider(names_from = final_name, values_from = flanking) %>%
    select(-Indiv) %>%
    as.matrix()

  gt_size <- gt_red %>% group_by(final_name) %>% summarise(size=mean(PERIOD)) %>% tibble::column_to_rownames("final_name")
  sizes <- factor(gt_size[unique(gt_red$final_name),], levels=c(3,4,5,6))
  names(sizes) <- colnames(mat)

  n_sz <- length(unique(gt_size$size))
  col_sizes <- viridis::viridis(n_sz)
  # col_sizes = c( "#8cc8bc", "#7da7ea", "#5773c0", "#1d4497")
  names(col_sizes) <- sort(unique(gt_size$size))

  ht <- Heatmap(mat, name=name_mat, na_col = "grey90",
                heatmap_legend_param = list(direction = "horizontal"),
                cluster_columns = F,
                cluster_rows = dendro,
                row_split = row_split,
                row_title = NULL,
                col = MetBrewer::met.brewer("Tam"),
                # column_labels = rep("", length(colnames(mat))),
                column_names_rot = 60,
                row_gap = unit(3, "mm"),
                # bottom_annotation = HeatmapAnnotation(ggplot = anno_empty(height = unit(height_ggplot, "cm"))),
                top_annotation = HeatmapAnnotation(posterior = anno_boxplot(q_mat, ylim = c(0,1), size = unit(.5,"mm"),
                                                                            height = unit(1.5,"cm"),
                                                                            axis_param = list(at = seq(0,1,.5))),
                                                   stutter = anno_boxplot(stut_mat, ylim = c(0,1), size = unit(.5,"mm"),
                                                                          height = unit(1.5,"cm"),
                                                                          axis_param = list(at = seq(0,1,.5))),
                                                   flanking = anno_boxplot(flank_mat, size = unit(.5,"mm"), ylim = c(0,1),
                                                                           height = unit(1.5,"cm"),
                                                                           axis_param = list(at = seq(0,1,.5))),
                                                   missing_perc = anno_barplot(colSums(is.na(mat)/nrow(mat)), gp = gpar(fill = "grey")),
                                                   sizes = sizes,
                                                   col = list(sizes=col_sizes),
                                                   annotation_legend_param = list(sizes = list(nrow=1)),
                                                   gap = unit(2, "mm")),
                right_annotation = rowAnnotation(species = heatmap_matrix$taxon,
                                                 # project = project_bar,
                                                 col = list(species=col_sp),
                                                 annotation_name_rot = 90, gap = unit(2, "mm"),
                                                 annotation_legend_param = list(species = list(nrow=3, direction = "horizontal"))),
                                                                                # project = list(nrow=3, direction = "horizontal"))),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if (!is.na(mat[i, j])){
                    grid.text(sprintf("%1.0f", mat[i, j]), x, y, gp = gpar(fontsize = 8,
                                                                           col=ifelse(mat[i, j] < quantile(c(min(mat, na.rm = T), max(mat, na.rm = T)), .65),
                                                                                      "black", "white")))
                  } else {
                    NULL
                  }
                }
  )

  # if (coverage_panel==TRUE){
  #   # coverage across reference for each sample -> horizion plot
  #   df_cov <- read_delim("STR/depth_stats.txt", col_names = c("pos", "cov", "Indiv"))
  #
  #   df_cov <- df_cov %>%
  #     left_join(meta %>% select(Sample, Y_id, Species), by = c("Indiv"="Y_id")) %>%
  #     filter(Indiv %in% gt$Indiv) %>%
  #     group_by(Sample) %>%
  #     mutate(scaled_cov = cov/median(cov), mc=median(cov)) %>%
  #     mutate(scaled_cov = ifelse(scaled_cov>3, 3, scaled_cov)) %>%
  #     ungroup()
  #
  #   lt = lapply(rownames(mat), function(x) pull(df_cov[df_cov$Sample==x, "scaled_cov"]))
  #   len_Ychrom <- 8309886
  #
  #   cov_anno <- rowAnnotation(scaled_coverage = anno_horizon(lt, n_slice = 1, gp = gpar(pos_fill="grey"),
  #                                                       axis_param = list(at = (sort(unique(gt$START))/len_Ychrom)*400,
  #                                                                         # facing = "inside",
  #                                                                         labels = c(""), gp = gpar(lwd=.5))))
  #   ht <- ht + cov_anno + rowAnnotation(rn = anno_text(rownames(mat),
  #                                                       location = unit(0, "npc"), just = "left"))
  # }
  return(ht)
}

# plot_alleles <- function(gt) {
#   df_sp_counts <- get_matrix(gt) %>% count(Species)
#
#   geno_sites <- gt %>%
#     mutate(Key2 = fct_reorder(final_name, Key)) %>%
#     left_join(df_sp_counts) %>%
#     select(n, Species, Key2, length_allele) %>%
#     group_by(Key2, Species, length_allele) %>%
#     summarise(n_samp = n(), n_species = min(n)) %>%
#     mutate(n_samp = n_samp/n_species) %>%
#     ungroup(Species) %>%
#     arrange(desc(n_samp), .by_group = T)
#
#   geno_plot <- ggplot(geno_sites, aes(Key2, length_allele, fill=Species, size=n_samp)) +
#       ggbeeswarm::geom_beeswarm(pch=21, alpha=0.7, cex=.25, priority = "none") +
#       scale_size(range = c(.7, 3)) +
#       scale_fill_manual(values = col_species) +
#       theme(legend.position = "none",
#             axis.text.x = element_blank(),
#             axis.title.x = element_blank(),
#             axis.ticks.x = element_blank(),
#             axis.title.y = element_text(size=8))
#   geno_plot
# }


plot_eda_vcf <- function(gt, col_sp=col_species, q_threshold=.9, min_dst=.35, min_dflnk=.35) {

  df_limits <- tibble(y=c(q_threshold, min_dst, min_dflnk), name=c("gt_Q", "flank", "stut"))

  pivoted_gt <- gt %>%
    mutate(flank=gt_DFLANKINDEL/gt_DP, stut=gt_DSTUTTER/gt_DP) %>%
    select(Indiv, gt_Q, flank, stut, gt_FILTER, taxon) %>%
    mutate(filtered=ifelse(gt_Q<=q_threshold | stut > min_dst | flank > min_dflnk, TRUE, FALSE)) %>%
    tidyr::pivot_longer(!c(Indiv, gt_FILTER, filtered, taxon))

  post_filtering <- pivoted_gt %>%
      ggplot(aes(forcats::fct_reorder(Indiv, as.numeric(factor(taxon, levels = names(col_sp)))),
                 value, color=taxon)) +
      geom_boxplot()+
      scale_color_manual(values=col_sp) +
      scale_y_continuous(expand = c(0,0.05)) +
      geom_hline(aes(yintercept=y), data = df_limits) +
      facet_grid(rows = vars(factor(name, levels=c("gt_Q", "flank", "stut")))) +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.ticks.x = element_blank(), legend.position="right")

  filtered_hipstr <- gt %>%
    ggplot(aes(forcats::fct_reorder(Indiv, as.numeric(factor(taxon, levels = names(col_sp)))), fill=gt_FILTER)) +
    geom_bar() +
    scale_fill_manual(values = c("#a00e00","#d04e00","#f6c200","#0086a8","green", "coral", "brown2")) +
    scale_y_continuous(expand = c(0,0.05)) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), legend.position="right")


  final_sites <- gt %>%
    mutate(flank=gt_DFLANKINDEL/gt_DP, stut=gt_DSTUTTER/gt_DP) %>%
    mutate(filtered=factor(ifelse(is.na(gt_GB), TRUE, FALSE), levels=c(TRUE, FALSE))) %>%
    ggplot(aes(forcats::fct_reorder(Indiv, as.numeric(factor(taxon, levels = names(col_sp)))), alpha=filtered)) +
    geom_bar() +
    scale_y_continuous(expand = c(0,0.05)) +
    labs(x="Sample") +
    theme(axis.text.x = element_text(angle=90, hjust=0.95, vjust = 0.5), legend.position="right")

  missing_samples <- filtered_hipstr / post_filtering /final_sites + plot_layout(heights = c(1,4,1)) + plot_annotation(tag_levels = 'A')
  missing_samples
}
