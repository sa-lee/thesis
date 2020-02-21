## ----plyranges-setup, include = FALSE------------------------
library(BiocFileCache)
library(ggbio)
library(Biostrings)
library(plyranges)
library(AnnotationHub)
# bam and snp array resources are stored in BiocFileCache
bfc <- BiocFileCache(cache = here::here("data"))
h1_bam_sorted <- bfcrpath(bfc, "GSM433167_BI.H3K27me3.sorted.bam")
h1_bam_sorted_index <- bfcrpath(bfc, "GSM433167_BI.H3K27me3.sorted.bam.bai")
h1_snp_array_path <- bfcrpath(bfc, "h1_snp_array") # this is an rds file
# AnnotationHub resource for T-cell BigWig file
bw_file <- AnnotationHub()[["AH33458"]]

## ----load-bw---------------------------------------
library(plyranges)
chr10_ranges <- bw_file %>%
  get_genome_info() %>%
  filter(seqnames == "chr10")


## ----chr10-scores-read-----------------------------
chr10_scores <- bw_file %>%
  read_bigwig(overlap_ranges = chr10_ranges) %>%
  set_genome_info(genome = "hg19")
chr10_scores


## ----chr10-scores-reduce---------------------------
all_peaks <- chr10_scores %>%
  filter(score > 8) %>%
  reduce_ranges(score = max(score))
all_peaks


## ----max-score-region------------------------------
chr10_max_score_region <- chr10_scores %>%
  filter(score == max(score)) %>%
  anchor_center() %>%
  mutate(width = 5000)


## ----peak_region-----------------------------------
peak_region <- chr10_scores %>%
  join_overlap_inner(chr10_max_score_region)
peak_region


## ----peak-viz--------------------------------------
ggbio::autoplot(peak_region,
                aes(y = score.x),
                geom = "blank") +
  geom_area() +
  ylab("score")


## ----bins------------------------------------------
bam <- read_bam(h1_bam_sorted, index = h1_bam_sorted_index)
locations <- bam %>%
  get_genome_info()


## ----read-bam, cache = TRUE---------------------------------
alignments <- bam %>%
  filter_by_overlaps(locations) %>%
  select(seq) %>%
  mutate(
    score = as.numeric(letterFrequency(seq, "GC", as.prob = TRUE))
  ) %>%
  select(score)
alignments


## ----read-bam-summary, cache = TRUE---------------------------------
bins <- locations %>%
  tile_ranges(width = 10000L)

alignments_summary <- bins %>%
  join_overlap_inner(alignments) %>%
  disjoin_ranges(n = n(), avg_gc = mean(score))
alignments_summary


##----load-array------------------------------------
h1_snp_array <- readRDS(h1_snp_array_path)


##----process-array------------------------------------
h1_snp_array <- h1_snp_array %>%
  filter(!(ref %in% c("I", "D")), seqnames != "M") %>%
  mutate(transition = (ref %in% c("A", "G") & alt %in% c("G","A"))|
                      (ref %in% c("C","T") & alt %in% c("T", "C")))


##----titv------------------------------------
ti_tv_results <- h1_snp_array %>%
  group_by(seqnames) %>%
  summarize(n_snps = n(),
            ti_tv = sum(transition) / sum(!transition))
ti_tv_results


## ----titv-viz, echo = FALSE----
ti_tv_results %>%
  as.data.frame() %>%
  ggplot(aes(forcats::fct_reorder(seqnames, ti_tv), ti_tv)) +
  geom_point() +
  geom_hline(yintercept = 3, colour = "white", size = 4) +
  labs(x = "Chromosome", y = "Transition Transversion Ratio") +
  coord_flip()

