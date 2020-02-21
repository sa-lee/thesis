# download and prepare BAM file example

library(BiocFileCache)
bfc <- BiocFileCache(cache = here::here("data"))

# ftp link on GEO
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM433nnn/GSM433167/suppl/GSM433167%5FBI%2EH3K27me3%2Ebam"
fpath <- "GSM433167_BI.H3K27me3.bam"

res <- bfcquery(bfc, fpath)

if (bfccount(res) == 0) {
  download.file(url = url, destfile = fpath)
  ans <- bfcadd(bfc, rname = fpath, action = "move" )
} else {
  ans <- bfcrpath(bfc, fpath)
}


sortedBam <- "GSM433167_BI.H3K27me3.sorted.bam"
res <- bfcquery(bfc, sortedBam)

if (bfccount(res) == 0) {
  Rsamtools::sortBam(ans, sub("\\.bam", "", sortedBam))
  Rsamtools::indexBam(sortedBam)
  indexBai <- paste0(sortedBam, ".bai")
  bfcadd(bfc, rname = indexBai, action = "move", rtype = "local")
  ans <- bfcadd(bfc, rname = sortedBam, action = "move", rtype = "local")
} else {
  ans <- bfcrpath(bfc, sortedBam)
}
