# Prepare GRanges object for SNP array resource
library(readr)
library(BiocFileCache)
library(plyranges)

bfc <- BiocFileCache(cache = here::here("data"))

rname <- "h1_snp_array"
res <- bfcquery(bfc, rname)

# if the resource doesn't exist start downloading files
if (bfccount(res) == 0) {
  array_info <- here::here("data-raw", "GPL18952_HumanOmni25M-8v1-1_B.annotated.txt.gz")
  if (!file.exists(array_info))
    download.file(
      url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL18952&format=file&file=GPL18952%5FHumanOmni25M%2D8v1%2D1%5FB%2Eannotated%2Etxt%2Egz",
      destfile = array_info
    )
  snp_info <- here::here("data-raw", "GSM1463263_JS-ESCH1.txt.gz")
  if (!file.exists(snp_info))
    download.file(
      url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1463263&format=file&file=GSM1463263%5FJS%2DESCH1%2Etxt%2Egz",
      destfile = snp_info
    )
  # prepare a GRanges object from array data
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
  saveRDS(h1_snp_array, bfcnew(bfc, rname))
} else {
  bfcrpath(bfc, rname)
}

