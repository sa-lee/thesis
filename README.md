# PhD thesis

This repo contains source files for my PhD thesis at Monash University.

The R packages used in this thesis can be installed via

```r
BiocManager::install("sa-lee/thesis")
```

## Clone with `git-lfs`

To clone this repo, you need to first download and install a git plugin called [`git-lfs`](https://git-lfs.github.com) for versioning large files, and set up Git LFS using command `git lfs install` in console.

## Directories

* `R/`: R code to reproduce tables, figures and analysis.
* `Rmd/`: R Markdown source documents for thesis document.
* `data/`: Cleaned data used for thesis document.
* `data-raw/`: R code to generate data in `data/`.
* `img/`: Images made with other tools to illustrate ideas. 
* `bib/`: Bibliography files.
* `template/`: Monash thesis template from [robjhydman/MonashThesis](https://github.com/robjhyndman/MonashThesis).


## Acknowledgements

Shout out to [Earo Wang](https://earo.me), who initially created this amazing
reproducible thesis template.
