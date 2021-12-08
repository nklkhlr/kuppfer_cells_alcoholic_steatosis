## Loss of Kupffer cells in alcoholic steatosis is associated with increased hepatic ceramide and free cholesterol levels

Code used in the paper by [K&#246;hler,  H&#246;ring, Czepukojc, Rose, Buechler et al.]() to generate ensemble biclustering results and figures.

#### Files
* DataPrep.R: Data pre-processing
* MoSBi.R: (Ensemble) bicluster computation
* PaperFigures.R: Visualisation of ensemble biclustering results

To install the [MoSBi package](https://www.biorxiv.org/content/10.1101/2021.09.30.462567v1) run

```r
if (!"devtools" %in% installed.packages()) install.packages("devtools")
devtools::install_github("tdrose/mosbi")
```

Please check the data analysis section in the paper for R version used and other package details.

