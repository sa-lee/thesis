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

## Reproducibility

The environment and R packages used to construct this thesis
can be recovered using the **renv** package. Run the following
code to install the packages used in this thesis:

```r
# install.packages("renv")
renv::restore()
```

The thesis can be compiled into both html and pdf format as follows:

```zsh
# uncomment if cache isn't setup
# make hub
# make bfcbam
# make bfcsnp
make both
```

## Caching local resources 

The `data-raw/` directory provides scripts used to download and 
cache large/raw files contained in this thesis. These are then managed
and transferred into R friendly data formats and managed with
**BiocFileCache** and **AnnotationHub**. 

Processed data files are stored as `.rds` and are located in the `data` directory. 

These files can be created via **BiocFileCache**
by running the following at the command line:

```zsh
make hub
make bfcbam
make bfcsnp
```


## Project Structure

* `Rmd/`: R Markdown source documents for thesis document.
* `scripts/`: R code to reproduce tables, figures and analyses for each chapter.
* `renv/`: Package environment captured by **renv**
* `data/`: Cleaned data used for thesis document.
* `data-raw/`: R code to generate data in `data/`..
* `img/`: Images made with other tools to illustrate ideas. 
* `video/`: Videos made with other tools to illustrate ideas. These are in mp4 format, screenshots are provided in the pdf version.
* `bib/`: Bibliography files.
* `template/`: Monash thesis template from [robjhydman/MonashThesis](https://github.com/robjhyndman/MonashThesis).
* `docs/`: the compiled thesis as pdf and website.

## Acknowledgements

Thank you [Earo Wang](https://earo.me), who initially created this  reproducible thesis template that I have forked and modified.
