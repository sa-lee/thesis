# PhD thesis 

This repo contains source files for my PhD thesis, entitled 

> Fluent statistical computing interfaces for biological data analysis    


## Summary

This work contributes new computing interfaces for the analysis of data 
collected from technologies that measure changes in biology, like gene 
expression. The first interface is a new tool that unites data measured along
the genome under a common language of transformation. This interface enables
analysts to interpret and combine multiple results from different technologies
in a clear and simple way. The second interface provides a framework for 
identifying interesting regions along the genome. The third interface empowers
analysts to interrogate complex models using visual diagnostics. All 
are provided as open source software packages. 

## Caching local resources 

The `data-raw/` directory provides scripts used to download and 
cache large/raw files contained in this thesis. These are then managed
and transferred into R friendly data formats and managed with
`BiocFileCache`. 

Processed data files are stored as `.rds` and are located in the `data`
directory. Since some files are greater than 50Mb, in order to clone this 
repository you will need to install `git lfs` and run `git lfs clone https://github.com/sa-lee/thesis`

## Reproducibility

The environment and R packages used to construct this thesis
can be recovered using the `renv` package. Run the following
code to install the packages used in this thesis:

```r
# install.packages("renv")
renv::restore()
```

## Project Structure

* `Rmd/`: R Markdown source documents for thesis document.
* `scripts/`: R code to reproduce tables, figures and analyses for each chapter.
* `data/`: Cleaned data used for thesis document.
* `data-raw/`: R code to generate data in `data/`..
* `img/`: Images made with other tools to illustrate ideas. 
* `bib/`: Bibliography files.
* `template/`: Monash thesis template from [robjhydman/MonashThesis](https://github.com/robjhyndman/MonashThesis).
* `video/`: Videos made with other tools to illustrate ideas. These are in mp4 format, screenshots are provided in the pdf version.
* `docs/`: the compiled thesis as pdf and website.

## Acknowledgements

Thank you [Earo Wang](https://earo.me), who initially created this  reproducible thesis template that I have forked and modified.
