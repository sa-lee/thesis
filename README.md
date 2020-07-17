# PhD thesis 

This repo contains source files for my PhD thesis, entitled 

> Fluent statistical computing interfaces for biological data analysis    

The R packages used in this thesis can be installed 
via the BiocManager package:

```r
BiocManager::install("sa-lee/thesis")
```

## Caching local resources 

The `data-raw/` directory provides scripts used to download and 
cache large/raw files contained in this thesis. Thesea are then managed
and transfered into R friendly data formats and managed with
`BiocFileCache`. 

Processed data files are stored as `.rds` and are located in the `data`
directory. Since some files are greater than 50Mb, in order to clone this repository you will need to install `git lfs` and run `git lfs clone https://github.com/sa-lee/thesis`


## Project Structure

* `Rmd/`: R Markdown source documents for thesis document.
* `scripts/`: R code to reproduce tables, figures and analyses for each chapter.
* `data/`: Cleaned data used for thesis document.
* `data-raw/`: R code to generate data in `data/`..
* `img/`: Images made with other tools to illustrate ideas. 
* `bib/`: Bibliography files.
* `template/`: Monash thesis template from[robjhydman/MonashThesis](https://github.com/robjhyndman/MonashThesis).
* `docs/`: the compiled thesis and pdf

Videos are rendered with `media9` for PDF and or hosted on vimeo for html.


## Acknowledgements

Shout out to [Earo Wang](https://earo.me), who initially created this amazing
reproducible thesis template that I forked and modified.
