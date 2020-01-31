## ----setup, include = FALSE------------------------
library(tibble)
library(knitr)
opts_chunk$set(message = FALSE,
               warning = FALSE,
               comment = "#>",
               fig.width = 5,
               fig.height = 3,
               fig.align = "center",
               fig.path = "./diagrams/")
# for peak figures
library(ggbio)

# for bam example
library(Biostrings)
library(plyranges)

# -- Data from the Human Epigenomics RoadMap consortium
# retrieving a BigWigFile from AnnotationHub

if (!dir.exists("./data")) dir.create("./data")

library(AnnotationHub)
q <- AnnotationHub()
# T-cell BW file
bw_file <- q[["AH33458"]]

# a BamFile from the H1 cell line
h1_bam <- "./data/GSM433167_BI.H3K27me3.bam"
if (!file.exists(h1_bam))
  download.file(
    url = "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM433nnn/GSM433167/suppl/GSM433167%5FBI%2EH3K27me3%2Ebam",
    destfile = h1_bam
  )

h1_bam_sorted <- "./data/GSM433167_BI.H3K27me3.sorted.bam"

# sort and index bam
if (!file.exists(h1_bam_sorted)) {
  Rsamtools::sortBam(h1_bam, "./data/GSM433167_BI.H3K27me3.sorted")
  Rsamtools::indexBam(h1_bam_sorted)
}

# array data for h1 cell line
# illumina annotations
array_info <- "data/GPL18952_HumanOmni25M-8v1-1_B.annotated.txt.gz"
if (!file.exists(array_info))
  download.file(
    url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL18952&format=file&file=GPL18952%5FHumanOmni25M%2D8v1%2D1%5FB%2Eannotated%2Etxt%2Egz",
    destfile = array_info
  )

snp_info <- "data/GSM1463263_JS-ESCH1.txt.gz"
if (!file.exists(snp_info))
  download.file(
    url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1463263&format=file&file=GSM1463263%5FJS%2DESCH1%2Etxt%2Egz",
    destfile = snp_info
  )



## ----load-bw, cache = TRUE-------------------------
library(plyranges)
chr10_ranges <- bw_file %>%
  get_genome_info() %>%
  filter(seqnames == "chr10")


## --------------------------------------------------
chr10_scores <- bw_file %>%
  read_bigwig(overlap_ranges = chr10_ranges) %>%
  set_genome_info(genome = "hg19")
chr10_scores


## --------------------------------------------------
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


## ----peak-viz, echo = FALSE, fig.cap = "The final result of the \\texttt{plyranges} operations to find a 5000nt region surrounding the peak of normalised coverage scores over chromosome 10, displayed as a density plot.", cache = TRUE----
ggbio::autoplot(peak_region,
                aes(y = score.x),
                geom = "blank") +
  geom_area() +
  ylab("score")


## ----bins------------------------------------------
locations <- h1_bam_sorted %>%
  read_bam() %>%
  get_genome_info()


## ---- cache = TRUE---------------------------------
alignments <- h1_bam_sorted %>%
  read_bam() %>%
  filter_by_overlaps(locations) %>%
  select(seq) %>%
  mutate(
    score = as.numeric(letterFrequency(seq, "GC", as.prob = TRUE))
  ) %>%
  select(score)
alignments


## ---- cache = TRUE---------------------------------
bins <- locations %>%
  tile_ranges(width = 10000L)

alignments_summary <- bins %>%
  join_overlap_inner(alignments) %>%
  disjoin_ranges(n = n(), avg_gc = mean(score))
alignments_summary


## ---- echo = FALSE, warning = FALSE, cache = TRUE----
# prepare a GRanges object from array data
library(readr)
# annotation information available here
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL18952
# chromosome zero are 'problem' probes and should be filtered
anno <- read_tsv(array_info,
                 trim_ws = TRUE,
                 col_types = c("cccc----")) %>%
  dplyr::filter(Chr != "0")

# genotypes for H1 cell line, data starts on line 10
# data obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1463263
geno <- read_tsv(snp_info, skip = 10)

# rename columns
geno <- geno %>%
  dplyr::select(name = `SNP Name`,
         score = `GC Score`,
         baf = `B Allele Freq`,
         logR = `Log R Ratio`)

# merge everything
complete_vranges_df <- dplyr::inner_join(geno, anno, by = c("name" = "Name")) %>%
  dplyr::rename(alleles = Alleles, seqnames = Chr, start = MapInfo)

h1_snp_array <- as_granges(complete_vranges_df,
                           start = as.integer(start),
                           end = as.integer(start))
h1_snp_array <- h1_snp_array %>%
  mutate(ref = stringr::str_replace(alleles, "\\[(.)/.*", "\\1"),
         alt = stringr::str_replace(alleles, ".*/(.)\\]", "\\1"))


## ---- cache = TRUE---------------------------------
h1_snp_array <- h1_snp_array %>%
  filter(!(ref %in% c("I", "D")), seqnames != "M") %>%
  mutate(transition = (ref %in% c("A", "G") & alt %in% c("G","A"))|
                      (ref %in% c("C","T") & alt %in% c("T", "C")))


## ---- cache = TRUE---------------------------------
ti_tv_results <- h1_snp_array %>%
  group_by(seqnames) %>%
  summarize(n_snps = n(),
            ti_tv = sum(transition) / sum(!transition))
ti_tv_results


## ----titv-viz, echo = FALSE, fig.cap = "The final result of computing quality control metrics over the SNP array data with \\texttt{plyranges}, displayed as a dot plot. Chromosomes are ordered by their estimated transition-transversion ratio. A white reference line is drawn at the expected ratio for a human exome."----
ti_tv_results %>%
  as.data.frame() %>%
  ggplot(aes(forcats::fct_reorder(seqnames, ti_tv), ti_tv)) +
  geom_point() +
  geom_hline(yintercept = 3, colour = "white", size = 4) +
  labs(x = "Chromosome", y = "Transition Transversion Ratio") +
  coord_flip()

