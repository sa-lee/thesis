## ----setup, include = FALSE------------------------
# due to bookdown loading everything at once, we
# first detach and reattach to avoid conflicts
detach("package:plyranges", unload = TRUE)
library(fluentGenomics)
library(SummarizedExperiment)
library(tximeta)
library(DESeq2)
library(readr)
library(dplyr)

path <- system.file("extdata", package="macrophage")

makeLinkedTxome(
  indexDir=file.path(path, "gencode.v29_salmon_0.12.0"),
  source="Gencode",
  organism="Homo sapiens",
  release="29",
  genome="GRCh38",
  fasta="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.transcripts.fa.gz",
  gtf=file.path(path, "gencode.v29.annotation.gtf.gz"), # local version
  write=FALSE
)
# # alias to plyranges
# n_distinct <- plyranges::n_distinct
# n <- plyranges::n

## ----setdir----------------------------------------
path <- system.file("extdata", package="macrophage")


## ----coldata-rna-----------------------------------
library(readr)
library(dplyr)
colfile <- file.path(path, "coldata.csv")
coldata <- read_csv(colfile) %>%
  select(
    names,
    id = sample_id,
    line = line_id,
    condition = condition_name
  ) %>%
  mutate(
    files = file.path(path, "quants", names, "quant.sf.gz"),
    line = factor(line),
    condition = relevel(factor(condition), "naive")
  )
glimpse(coldata)


## ----tximeta-run-----------------------------------
library(SummarizedExperiment)
library(tximeta)
se <- tximeta(coldata, dropInfReps=TRUE)
se


## ----gse-------------------------------------------
gse <- summarizeToGene(se)


## ----setup-deseq-----------------------------------
library(DESeq2)
dds <- DESeqDataSet(gse, ~line + condition)
# filter out lowly expressed genes
# at least 10 counts in at least 6 samples
keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep,]


## ----fit-model-------------------------------------
dds <- DESeq(dds)


## ----results-DFrame--------------------------------
res <- results(dds,
               contrast=c("condition","IFNg","naive"),
               lfcThreshold=1, alpha=0.01)


## ----ma-plot---------------------------------------
summary(res)
DESeq2::plotMA(res, ylim=c(-10,10))


## ----results-GRanges-------------------------------
suppressPackageStartupMessages(library(plyranges))
de_genes <- results(dds,
                    contrast=c("condition","IFNg","naive"),
                    lfcThreshold=1,
                    format="GRanges") %>%
  names_to_column("gene_id")
de_genes


## ----de-genes--------------------------------------
de_genes <- de_genes %>%
  filter(padj < 0.01) %>%
  select(
    gene_id,
    de_log2FC = log2FoldChange,
    de_padj = padj
  )


## ----not-de-genes----------------------------------
other_genes <- results(dds,
                       contrast=c("condition","IFNg","naive"),
                       lfcThreshold=1,
                       altHypothesis="lessAbs",
                       format="GRanges") %>%
  filter(padj < 0.01) %>%
  names_to_column("gene_id") %>%
  select(
    gene_id,
    de_log2FC = log2FoldChange,
    de_padj = padj
  )


## ----load-peaks------------------------------------
library(fluentGenomics)
peaks


## ----filter-peaks----------------------------------
da_peaks <- peaks %>%
  filter(da_padj < 0.01)


## ----slice-example---------------------------------
size <- length(de_genes)
slice(other_genes, sample.int(n(), size))


## ----boot-set-01-----------------------------------
# set a seed for the results
set.seed(2019-08-02)
subsamp_genes <- replicate(10,
                           slice(other_genes, sample.int(n(), size)),
                           simplify = FALSE)


## ----boot-set-02-----------------------------------
subsamp_genes <- bind_ranges(subsamp_genes, .id = "resample")


## ----combine-results-------------------------------
all_genes <- bind_ranges(
  de=de_genes,
  not_de = subsamp_genes,
  .id="origin"
) %>%
  mutate(
    origin = factor(origin, c("not_de", "de")),
    resample = ifelse(is.na(resample), 0L, as.integer(resample))
  )
all_genes


## ----resize-01-------------------------------------
all_genes <- all_genes %>%
  anchor_5p() %>%
  mutate(width = 1)


## ----resize-02-------------------------------------
all_genes <- all_genes %>%
  anchor_center() %>%
  mutate(width=2*1e4)


## ----olap-join-------------------------------------
genes_olap_peaks <- all_genes %>%
  join_overlap_left(da_peaks)
genes_olap_peaks


## ----reduce-ex01-----------------------------------
gene_peak_max_lfc <- genes_olap_peaks %>%
  group_by(gene_id, origin)  %>%
  reduce_ranges_directed(
    peak_count = sum(!is.na(da_padj)) / n_distinct(resample),
    peak_max_lfc = max(abs(da_log2FC))
  )


## ----boxplot, fig.cap = "(ref:boxplot)"------------
library(ggplot2)
gene_peak_max_lfc %>%
  filter(peak_count > 0) %>%
  as.data.frame() %>%
  ggplot(aes(origin, peak_max_lfc)) +
  geom_boxplot()


## ----summarize-ex01--------------------------------
origin_peak_lfc <- genes_olap_peaks %>%
  group_by(origin) %>%
  summarize(
    peak_count = sum(!is.na(da_padj)) / n_distinct(resample),
    lfc1_peak_count =
      sum(abs(da_log2FC) > 1, na.rm=TRUE) / n_distinct(resample),
    lfc2_peak_count =
      sum(abs(da_log2FC) > 2, na.rm=TRUE) / n_distinct(resample)
  )
origin_peak_lfc


## ----pivot-enrich----------------------------------
library(tidyr)
origin_peak_lfc %>%
  as.data.frame() %>%
  pivot_longer(cols = -origin) %>%
  pivot_wider(names_from = origin, values_from = value) %>%
  mutate(enrichment = de / not_de)


## ----reduce-summarize------------------------------
genes_olap_peaks %>%
  group_by(gene_id, origin, resample) %>%
  reduce_ranges_directed(
    lfc1 = sum(abs(da_log2FC) > 1, na.rm=TRUE),
    lfc2 = sum(abs(da_log2FC) > 2, na.rm=TRUE)
  ) %>%
  group_by(origin) %>%
  summarize(
    lfc1_gene_count = sum(lfc1 > 0) / n_distinct(resample),
    lfc1_peak_count = sum(lfc1) / n_distinct(resample),
    lfc2_gene_count = sum(lfc2 > 0) / n_distinct(resample),
    lfc2_peak_count = sum(lfc2) / n_distinct(resample)
  )


## ----count-fn--------------------------------------
count_if_above_threshold <- function(var, thresholds) {
  lapply(thresholds, function(.) sum(abs(var) > ., na.rm = TRUE))
}


## ----thresholds------------------------------------
thresholds <- da_peaks %>%
  mutate(abs_lfc = abs(da_log2FC)) %>%
  with(
    seq(min(abs_lfc), max(abs_lfc), length.out = 100)
  )


## ----reduce-ex02-----------------------------------
genes_peak_all_thresholds <- genes_olap_peaks %>%
  group_by(gene_id, origin, resample) %>%
  reduce_ranges_directed(
    value = count_if_above_threshold(da_log2FC, thresholds),
    threshold = list(thresholds)
  )
genes_peak_all_thresholds


## ----expand-summarize------------------------------
origin_peak_all_thresholds <- genes_peak_all_thresholds %>%
  expand_ranges() %>%
  group_by(origin, threshold) %>%
  summarize(
    gene_count = sum(value > 0) / n_distinct(resample),
    peak_count = sum(value) / n_distinct(resample)
  )
origin_peak_all_thresholds


## ----line-chart-data-------
origin_threshold_counts <- origin_peak_all_thresholds %>%
  as.data.frame() %>%
  pivot_longer(cols = -c(origin, threshold),
               names_to = c("type", "var"),
               names_sep = "_",
               values_to = "count") %>%
  select(-var)

## ----line-chart, fig.cap = "(ref:linechart)"-------
origin_threshold_counts %>%
  filter(type == "peak") %>%
  pivot_wider(names_from = origin, values_from = count) %>%
  mutate(enrichment =  de / not_de) %>%
  ggplot(aes(x = threshold, y = enrichment)) +
  geom_line() +
  labs(x = "logFC threshold", y = "Relative Enrichment")


## ----line-chart2, fig.cap = "(ref:linechart2)"-----
origin_threshold_counts %>%
  ggplot(aes(x = threshold,
             y = count + 1,
             color = origin,
             linetype = type)) +
  geom_line() +
  scale_y_log10()
