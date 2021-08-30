Bayesian Phylogeographic and Ecological Clustering
================
Achyut K Banerjee
August 23, 2021

**Premise:** To investigate the environmental influence of the sampling
locations on population structure of a mangrove species Aegiceras
corniculatum in the Indo-West Pacific (populations examined = 44).

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/mol.bio/bpec\]

*Install and load package*

``` r
library(BPEC)
```

*Read and format data*

``` r
#Note: Here the covariates considered are: lat,long,bio1,2,7,10,12,13 and DEM. Other models can be fit with more covariates like soil parameters, vegetation parameters etc.#
coordsLocs <- bpec.loadCoords("C:/Users/Achyut/Desktop/web/mol.bio/bpec/coordsLocsFile.txt", header = TRUE) #coordinates and covariates of the populations#
coordsLocs<-as.data.frame(coordsLocs)
seq<-bpec.loadSeq("C:/Users/Achyut/Desktop/web/mol.bio/bpec/haplotypes.nex") #haplotypes present in each population#
```

*Assess clustering pattern*

``` r
bpecout <- bpec.mcmc(seq, coordsLocs, maxMig = 6, iter = 10000000,
                     ds = 0, postSamples = 1000, dims = 9)
bpec.Tree <- bpec.treePlot(bpecout)
par(mar=c(0,0,0,0),pty="m",mfrow=c(2,4))
bpec.covariatesPlot(bpecout)
par(mar=c(0,0,0,0),pty="m",mfrow=c(1,1))
bpec.contourPlot(bpecout,GoogleEarth=0) 
output.clust(bpecout)
```

**END**

**References**

-   \[BPEC package\]
    (<https://cran.r-project.org/web/packages/BPEC/BPEC.pdf>)
