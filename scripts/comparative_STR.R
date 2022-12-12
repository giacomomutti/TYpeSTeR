packages <- c('optparse','ggplot2', 'magrittr', 'dplyr', 'tidyr')
new_packages <- packages[!packages %in% installed.packages()[,'Package']]

if(length(new_packages)>0) {
  print(paste("please install", new_packages))
  quit(status = 1)
}

library(optparse)

option_list = list(
  make_option(c("-b", "--blast"), type="character", default="blast_all_v_all_Y_str.csv",
              help="blast output from find_homology.sh. Deefault blast_all_v_all_Y_str.csv", dest = "blast"),
  make_option(c("-r", "--regions"), type="character", default="regions_metadata.txt",
              help="regions file output from find_homology.sh. Default regions_metadata.txt", dest = "regions"),
  make_option(c("-f","--flank"), type="integer", default=200,
              help="flanking regions used in find_homology.sh Default 200", dest = "flank"),
  make_option(c("-c","--cov"), type="double", default=.65,
              help="minimum query coverage. Default .65", dest = "minqcov"),
  make_option(c("-i","--ident"), type="double", default=80,
              help="min identity percentage. Deafult 80", dest = "minident")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

suppressPackageStartupMessages({
library(ggplot2)
library(magrittr)
library(dplyr)
  })

theme_set(theme_bw())
# file obtained with find_homology.sh
regions_meta <- read.delim(opt$regions, header = F, col.names = c("reference", "start", "end", "length", "name"))
# blast of panu regions +- 100
flank <- opt$flank

blast_col <- c("qseqid", "sseqid", "pident", "length_hit", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
# remove hits between the same sequence
blast_str <- read.delim(opt$blast, header = F, col.names = blast_col) %>% 
  filter(qseqid!=sseqid)

size <- length(unique(regions_meta$reference)) 

# min identity and query coverage
min_ident <- opt$minident
min_qcov <- opt$minqcov

# compute if regions are multi copy, within species and homologous
blast_str <- blast_str %>% 
  left_join(regions_meta, by = c("qseqid"="name")) %>% 
  left_join(regions_meta, by = c("sseqid"="name")) %>% 
  mutate(query_cov = length_hit/(pmax(end.x-start.x, end.y-start.y)+2*flank+1), 
         # query_cov=ifelse(query_cov>1, 1, query_cov),
         within = ifelse(reference.x==reference.y, TRUE, FALSE),
         comparison = paste(pmax(reference.x, reference.y), pmin(reference.x, reference.y), sep = "_vs_"),
         Multicopy = ifelse(reference.x==reference.y & query_cov>min_qcov, TRUE, FALSE),
         Homologous = ifelse(query_cov>=min_qcov & pident>= min_ident, TRUE, FALSE))

# multicopy regions
multi_df <- blast_str %>%
  group_by(qseqid) %>%
  summarise(Multicopy=ifelse(sum(Multicopy)>0, TRUE, FALSE)) %>%
  rename("qid"="qseqid") %>%
  ungroup()

# these are the multicopy str, we recover them in human except YCAIIa and b!!

dir.create("plots", showWarnings = FALSE)

plot_blast <- blast_str %>% 
  ggplot(aes(pident, query_cov, color=Multicopy, alpha=Homologous)) +
  geom_point() +
  geom_hline(yintercept = min_qcov) +
  geom_vline(xintercept = min_ident) +
  scale_alpha_manual(values = c(0.3,1)) +
  facet_grid(reference.x~reference.y) + 
  labs(x="% identity", y="Query coverage") +
  scale_color_manual(values = c("#3B9AB2", "#F21A00")) + 
  theme(legend.position = "bottom",
        text = element_text(family="Helvetica"),
        legend.box="vertical")

ggsave("plots/blast_Y_STRs.pdf", plot_blast, device="pdf", width=size*2, height=size*2)
