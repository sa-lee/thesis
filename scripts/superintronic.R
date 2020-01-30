## ---- superintronic-setup
library(plyranges)
library(BiocParallel)
library(superintronic)
library(ggplot2)
library(here)
theme_set(theme_bw())

# Plotting helper functions
set_plot_rng <- function(cvg, target, design, alpha) {
  if (is(cvg, "GRangesList")) {
    stopifnot(is(design, "DataFrame"))
    cvg_rng <- as(lapply(cvg, function(x) superintronic::join_parts(x, target)),
                  "GRangesList")
    md <- S4Vectors::DataFrame(
      Sample = S4Vectors::Rle(names(cvg_rng), lengths(cvg_rng)),
      CellLine = S4Vectors::Rle(design$CellLine, lengths(cvg_rng)),
      Kit = S4Vectors::Rle(as.factor(design$Kit), lengths(cvg_rng))
    )

    cvg_rng <- unlist(cvg_rng, use.names = FALSE)
    S4Vectors::mcols(cvg_rng) <- cbind(S4Vectors::mcols(cvg_rng), md)

    cvg_rng <- plyranges::mutate(cvg_rng,
                                 log_score = log2(score + 0.5),
                                 strand = BiocGenerics::strand(target),
                                 var = paste0("Kit: ", Kit, " Cellline: ", CellLine))
  }

  if (is(cvg, "GRanges")) {
    cvg_rng <- plyranges::filter_by_overlaps(cvg, target)
    cvg_rng <- plyranges::mutate(cvg_rng,
                                 strand = BiocGenerics::strand(target),
                                 var = paste0("Kit: ", Kit, " Cellline: ", CellLine))
  }

  if (alpha) {
    cvg_rng <- plyranges::disjoin_ranges_directed(group_by(cvg_rng, var),
                                                  log_score = BiocGenerics::mean(log_score),
                                                  feature_type = unlist(BiocGenerics::unique(feature_type)))

    cvg_rng <- plyranges::mutate(cvg_rng,
                                 alpha = dplyr::case_when(
                                   score > 3 & feature == "intron" ~ 1,
                                   score <= 3 & feature == "intron"~ 0.5,
                                   TRUE ~ 1
                                 ))
  }

  return(cvg_rng)
}


# pretty track plot
pretty_cov_plot <- function(cvg, parts, target, design = NULL, alpha = FALSE, ...) {

  if(is(target, "data.frame")) {
    target <- plyranges::filter(parts, gene_id == !!target$gene_id)
  }

  cvg_rng <- set_plot_rng(cvg, target, design, alpha)

  if (alpha) {
    p <- ggplot(as.data.frame(cvg_rng),
                aes(x = start, xend = end, y = 0, yend = log_score)) +
      geom_segment(aes(colour = feature_type, alpha = alpha)) +
      superintronic:::rescale_by_width(cvg_rng) +
      guides(alpha = FALSE)
  } else {
    p <- superintronic::view_coverage(cvg_rng,
                                      score = log_score,
                                      colour = feature_type,
                                      facets = dplyr::vars(var))
  }

  p <- p +
    ylim(0, NA) +
    guides(colour = FALSE) +
    scale_color_brewer(palette = "Dark2") +
    labs(subtitle = paste("Coverage over", target$gene_name),
         y = "Log coverage"
    ) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(hjust = 0, size = 8),
          strip.background = element_blank()) +
    expand_limits(y = 0)

  segments <- superintronic::unnest_parts(target)
  track <- view_segments(segments, colour = feature_type) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8))

  patchwork::wrap_plots(p,
                        track,
                        ncol = 1,
                        ...)

}

## ---- load-cvg
cvg <- readRDS(here("data", "superintronic", "parts-coverage.rds")) %>%
  mutate(log_score = log2(score + 0.5),
         CellLine = ifelse(CellLine == "HCC287", "HCC827", CellLine))

## ---- load-annotation
parts <- readRDS(here("data", "superintronic", "filtered-annotation.rds"))


## ---- coverage-histogram
# filter low coverage
cvg2 <- cvg %>%
  filter(CellLine == "HCC827") %>%
  group_by(Kit, gene_id) %>%
  filter(mean(log_score) >= log2(3.5)) %>%
  ungroup()

# compute max per sample
max_cov <- cvg2 %>%
  group_by(Sample, gene_id) %>%
  summarise(max_log_score = max(log_score))

# join it back on and add strand information
cvg2 <- cvg2 %>%
  mutate(strand = feature_strand,
         max_log_score = max_cov[match(
           paste0(Sample, gene_id),
           paste0(max_cov$Sample, max_cov$gene_id)
         ), "max_log_score"],
         normed_score = log_score / max_log_score
  )

# order by coordinates
cvg2 <- sort(cvg2)

# reshape via group split method in plyranges
by_gene <- cvg2 %>%
  dplyr::group_split(gene_id)

# for each gene section it into twenty bins
# reverse order for negative stranded genes
bins <- bplapply(by_gene,
               function(.x) {
                 unlist(
                   GenomicRanges::tile(
                     reduce_ranges_directed(.x),
                     n = 20
                   )
                 ) %>%
                   mutate(bin = ifelse(strand == "-", 20:1, 1:20))
               }) %>%
  as("GRangesList")

# overlap it with genes
olaps <- bplapply(
  seq_along(bins),
  function(i) join_overlap_intersect_directed(bins[[i]], by_gene[[i]])
  ) %>%
  as("GRangesList")


olaps <- unlist(olaps) %>%
  mutate(bin = Rle(bin))

# summarise over bins for each gene and feature
bin_means <- olaps %>%
  group_by(Sample, gene_id, feature_type, bin) %>%
  summarise(mn = mean(normed_score)) %>%
  dplyr::as_tibble()

## ---- categorise-gene-lengths
gene_lengths <- parts %>%
  select(gene_id, width, .drop_ranges = TRUE) %>%
  dplyr::as_tibble() %>%
  mutate(
    cat = dplyr::case_when(
      width >= quantile(width, 2/3) & width <= max(width) ~ "long",
      width >= quantile(width, 1/3) & width < quantile(width, 2/3) ~ "regular",
      width >= min(width) & width < quantile(width, 1/3) ~ "short"
    ),
    cat = factor(cat, levels = c("short", "regular", "long"))
  )

## ---- smoothed-histogram
bin_means %>%
  dplyr::left_join(gene_lengths) %>%
  group_by(Sample, cat, feature_type, bin) %>%
  summarise(lumpy = mean(mn), bumpy = var(mn)) %>%
  ungroup() %>%
  mutate(Kit = gsub("R[1-3]_HCC287_", "", Sample),
         feature_type = factor(feature_type, levels = c("intron", "exon"))) %>%
  ggplot(tbl, aes(x = bin, y = lumpy, colour = Kit, group = Sample)) +
  geom_line() +
  facet_grid(feature_type ~ cat, scales = "free_y") +
  scale_color_manual(values = c("#404040", "#bababa")) +
  theme(aspect.ratio = 1,
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        ) +
  labs(y = "Relative Log-coverage", x = "5'---Genebody---3'")


## ---- pre-mRNA examples
short <- c("ENSG00000136997.17", "ENSG00000117525.13")
long <- c("ENSG00000196937.10", "ENSG00000196428.12")
inx <- seq_len(length(short) + length(long))
names(inx) <- c(short, long)

## ---- short-genes
tracks <- vector("list", length(inx))
for (i in inx) {
  gene <- names(inx)[[i]]
  target <- filter(parts, gene_id == gene)
  tracks[[i]] <- pretty_cov_plot(cvg,
                                 parts,
                                 target,
                                 alpha = TRUE,
                                 heights = c(2, 0.25))
}

patchwork::wrap_plots(tracks[names(inx) %in% short],
                      ncol = 2,
                      guides = "keep")

## ---- long-genes
patchwork::wrap_plots(tracks[names(inx) %in% long],
                      ncol = 2,
                      guides = "keep")
