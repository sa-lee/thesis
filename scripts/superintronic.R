## ----setup, include=FALSE-------------------
# pre load libraries
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(superintronic))
library(ggplot2)

# use hexbinning for pairs plot
hex <- function(data, mapping, ...) {
  ggplot2::ggplot(data = data, mapping = mapping) +
    ggplot2::geom_hex(...) +
    ggplot2::scale_fill_viridis_c(direction = -1)
}

# plot function for pretty coverage traces
pretty_cov_plot <- function(cvg, parts, target, design = NULL, alpha = FALSE, ...) {

  if(is(target, "data.frame")) {
    target <- plyranges::filter(parts, gene_id == !!unique(target$gene_id))
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
    labs(subtitle = paste("Coverage over", target$gene_id),
         y = "Log coverage"
    ) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(hjust = 0, size = 8),
          strip.background = element_blank()) +
    expand_limits(y = 0)

  segments <- superintronic::unnest_parts(target)
  track <- view_segments(segments, colour = feature_type) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank())

  patchwork::wrap_plots(p,
                        track,
                        ncol = 1,
                        ...)

}


set_plot_rng <- function(cvg, target, design, alpha) {
  if (is(cvg, "GRangesList")) {
    stopifnot(is(design, "DataFrame"))
    cvg_rng <- as(lapply(cvg, function(x) superintronic::join_parts(x, target)),
                  "GRangesList")
    md <- S4Vectors::DataFrame(
      Sample = S4Vectors::Rle(names(cvg_rng), lengths(cvg_rng)),
      Genotype = S4Vectors::Rle(design$Genotype, lengths(cvg_rng)),
      Line = S4Vectors::Rle(as.factor(design$Line), lengths(cvg_rng))
    )

    cvg_rng <- unlist(cvg_rng, use.names = FALSE)
    S4Vectors::mcols(cvg_rng) <- cbind(S4Vectors::mcols(cvg_rng), md)

    cvg_rng <- plyranges::mutate(cvg_rng,
                                 log_score = log2(score + 0.5),
                                 strand = BiocGenerics::strand(target),
                                 var = paste0("Line: ", Line, " Genotype: ", Genotype))
  }

  if (is(cvg, "GRanges")) {
    cvg_rng <- plyranges::filter_by_overlaps(cvg, target)
    cvg_rng <- plyranges::mutate(cvg_rng,
                                 strand = BiocGenerics::strand(target),
                                 var = paste0("Line: ", Line, " Genotype: ", Genotype))
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
}

## ----prepare-design, echo= FALSE------------
design <- readRDS(here::here("data", "design.rds"))
design


## ----gff------------------------------------
library(superintronic)
parts <- readRDS(here::here("data", "parts.rds"))
parts


## ----default-filter, echo = FALSE-----------
parts_sub <- parts %>%
  mutate(exonic_parts = unique(exonic_parts),
         intronic_parts = unique(intronic_parts)) %>%
  filter(lengths(exonic_parts) > 1, n_olaps == 1, lengths(intronic_parts) >= 1)
parts_sub


## ---- load-features, echo = FALSE-----------
cvg <- readRDS(here::here("data", "complete-coverage.rds"))
cvg_over_features <- superintronic::join_parts(cvg, parts_sub)


## ----rangenostics-01------------------------
# compute intron/exon features
cvg_over_features <- cvg_over_features %>%
  mutate(log_score = log2(score + 0.5))

cvg_by_features <- cvg_over_features %>%
  group_by(Genotype, Line, gene_id, feature_type)

E_vals <- cvg_by_features %>%
  summarise(mn = Hmisc::wtd.mean(log_score, width),
            sd = sqrt(Hmisc::wtd.var(log_score, width)),
            raw_mn = Hmisc::wtd.mean(score, width),
            score = score,
            n_bases = width,
            start = min(start),
            end = max(end),
            seqnames = unlist(unique(seqnames)),
            strand = unlist(unique(feature_strand))) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)



## ----rangenostics-02------------------------
lumpy <- E_vals %>%
  filter(sd > 0, feature_type == "exon") %>%
  group_by(Genotype, Line) %>%
  summarise(
    E_mu = mean(mn),
    E_sd = mean(sd),
    E_raw = mean(raw_mn)
  )

lumpy

mcols(E_vals) <- cbind(mcols(E_vals),
                       lumpy[match(paste0(E_vals$Genotype, E_vals$Line),
                                   paste0(lumpy$Genotype, lumpy$Line)), -c(1,2)])

E_vals <- E_vals %>%
  mutate(
    n_bases_E_raw = sum(n_bases[score  > E_raw]),
  )


## ----rangenostics-03------------------------
rango <- E_vals %>%
  plyranges::select(-n_bases, -score, -seqnames, -start, -end, -width, -strand,
                    .drop_ranges = TRUE) %>%
  dplyr::as_tibble() %>%
  group_by(Genotype, Line, feature_type) %>%
  tidyr::nest() %>%
  mutate(
    smooth = lapply(data, function(x) {
      mgcv::gam(sd ~ s(mn), data = x)
    }),
    augment = lapply(smooth, broom::augment),
  ) %>%
  tidyr::unnest(data, augment) %>%
  dplyr::select(-dplyr::ends_with("1"), -.sigma)


## ----voom-like-plot-------------------------
# rangenostics feature_average vs feature_sd coloured by n of bases above
# average raw exon coverage
library(ggplot2)

voom_plot <- function(.x, .y) {
  p <- ggplot(data = .x, aes(x = mn, y = sd)) +
    geom_point() +
    geom_line(aes(y = .fitted), colour = "blue") +
    geom_vline(aes(xintercept = E_mu),
               data  = dplyr::distinct(.x, feature_type, E_mu)) +
    geom_hline(aes(yintercept = E_sd),
               data = dplyr::distinct(.x, feature_type, E_sd)) +
    facet_wrap(~feature_type) +
    labs(
      subtitle = paste(.y$Line, .y$Genotype),
      x = "mean log-coverage",
      y = "sd log-covreage"
    )
  fname <- here::here("img",
                      paste0(.y$Line, "-", .y$Genotype,
                             "voom-like-plot.png")
  )
  ggsave(fname, p)

}

rango %>%
  group_by(Genotype, Line) %>%
  dplyr::group_walk(voom_plot)



## ----rangenostics-pairs---------------------
all_features <- tidyr::gather(rango, "key", "value",
                              -Line, -Genotype, -gene_id, -feature_type,
                              -E_mu, -E_sd, -E_raw) %>%
  mutate(var = sub("\\.", "", paste0(feature_type, "_", key))) %>%
  plyranges::select(-feature_type, -key) %>%
  arrange(gene_id) %>%
  tidyr::spread(var, value)

pairs_plot <- function(.x, .y) {
  sub <- plyranges::select(.x,
                           exon_mn,
                           intron_mn,
                           intron_sd,
                           bases_above = intron_n_bases_E_raw)
  fname <- here::here("img",
                      paste0(.y$Line, "-", .y$Genotype,
                             "-pairs-plot.png")
  )

  p <- GGally::ggpairs(sub, columns = 1:4, lower = list(continuous = hex))
  ggsave(fname, p)
}

all_features %>%
  group_by(Genotype, Line) %>%
  dplyr::group_walk(pairs_plot)




## ----rangenositcs-hits----------------------
hits <- all_features %>%
  group_by(Genotype, Line) %>%
  filter(exon_mn > E_mu)

dplyr::count(hits)


## ---- echo  = FALSE-------------------------
sets <- hits %>%
  dplyr::distinct(gene_id) %>%
  arrange(gene_id) %>%
  mutate(value = 1L, var = paste(Genotype, Line)) %>%
  ungroup() %>%
  plyranges::select(gene_id, var, value) %>%
  tidyr::spread(var, value) %>%
  dplyr::mutate_at(2:5, .funs = ~ dplyr::if_else(is.na(.), 0L, .))

sets <- as.data.frame(sets)
library(UpSetR)
pdf(here::here("img", "olaps-superintronic-hits.pdf"))
upset(sets)
dev.off()


## ----all-features-out, include = FALSE------
saveRDS(all_features, here::here("data", "superintronic-features.rds"))
saveRDS(hits, here::here("data", "superintronic-hits.rds"))
readr::write_csv(hits, here::here("data", "zfish-hits.csv"))


## ----superintronic-coverage-plots, include = FALSE----
# look at combinations and save those coverage plots
# this is really clunky code but gets the job done
sets <- sets %>%
  mutate(combo =
           `Rnpc3_KO Cal` + `Rnpc3_KO ZM_Cal` +  `WT Cal` + `WT ZM_Cal` )

all_sets <- dplyr::filter(sets, combo == 4)

for (gene in seq_len(nrow(all_sets))) {
  track_plot <- pretty_cov_plot(cvg_over_features, parts_sub, all_sets[gene,], heights = c(0, 4, 0.25))
  ggsave(here::here("img",
                    paste0("all-sets-hit-", all_sets[gene, "gene_id"], ".pdf")),
         plot = track_plot
  )
}

three_sets <- dplyr::filter(sets, combo == 3)

for (gene in seq_len(nrow(three_sets))) {
  input <- three_sets[gene,]
  cols <- names(
    Filter(function(x) length(x) ==  1,
           lapply(input[, 2:5], function(x) which(x == 1))
    )
  )
  cols <- paste(gsub(" ", "-", cols), collapse = "_")
  track_plot <- pretty_cov_plot(cvg_over_features, parts_sub, input, heights = c(0, 4, 0.25))
  ggsave(here::here("img",
                    paste0("three-sets-hit-", cols, "-", input$gene_id, ".pdf")),
         plot = track_plot
  )
}

two_sets <- dplyr::filter(sets, combo == 2)

for (gene in seq_len(nrow(two_sets))) {
  input <- two_sets[gene,]
  cols <- names(
    Filter(function(x) length(x) ==  1,
           lapply(input[, 2:5], function(x) which(x == 1))
    )
  )
  cols <- paste(gsub(" ", "-", cols), collapse = "_")
  track_plot <- pretty_cov_plot(cvg_over_features, parts_sub, input, heights = c(0, 4, 0.25))
  ggsave(here::here("img",
                    paste0("two-sets-hit-", cols, "-", input$gene_id, ".pdf")),
         plot = track_plot
  )
}

one_sets <- dplyr::filter(sets, combo == 1)

for (gene in seq_len(nrow(one_sets))) {
  input <- one_sets[gene,]
  cols <- names(
    Filter(function(x) length(x) ==  1,
           lapply(input[, 2:5], function(x) which(x == 1))
    )
  )
  cols <- paste(gsub(" ", "-", cols), collapse = "_")
  track_plot <- pretty_cov_plot(cvg_over_features, parts_sub, input, heights = c(0, 4, 0.25))
  ggsave(here::here("img",
                    paste0("one-sets-hit-", cols, "-", input$gene_id, ".pdf")),
         plot = track_plot
  )
}
