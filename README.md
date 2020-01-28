# PhD thesis

This repo contains source files for my PhD thesis at Monash University.

The R packages used in this thesis can be installed via the BiocManager package:

```r
BiocManager::install("sa-lee/thesis")
```

## Caching local resources 

The `data-raw/` directory is used to download and cache large 
data files contained in this thesis (using a local cache created with
`BiocFileCache`). This requires a local internet connection in order
to perform the downloads. 


## Directories

* `scripts/`: R code to reproduce tables, 
    figures and analysis for each chapter.
* `Rmd/`: R Markdown source documents for thesis document.
* `data/`: Cleaned data used for thesis document.
* `data-raw/`: R code to generate data in `data/`.
* `img/`: Images made with other tools to illustrate ideas. 
* `bib/`: Bibliography files.
* `template/`: Monash thesis template from [robjhydman/MonashThesis](https://github.com/robjhyndman/MonashThesis).


## Acknowledgements

Shout out to [Earo Wang](https://earo.me), who initially created this amazing
reproducible thesis template.
