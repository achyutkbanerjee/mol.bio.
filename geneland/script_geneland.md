Genetic structure using Geneland
================
Achyut K Banerjee
August 23, 2021

**Premise:** To investigate the genetic structure (**by integrating the
spatial information in the clustering algorithm**) of a mangrove
associate species Pluchea indica in the Indo-West Pacific (populations
examined = 31, individuals = 348; SSR data).

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/mol.bio/geneland\]

*Install and load package*

``` r
library(Geneland)
```

*Read and format data*

``` r
geno <- read.table("C:/Users/Achyut/Desktop/web/mol.bio/geneland/data.txt",na.string="0")
coord <- read.table("C:/Users/Achyut/Desktop/web/mol.bio/geneland/loc.txt")
dim(coord)
nrow(geno)
ncol(geno)
plot(coord,xlab="Eastings",ylab="Northings",asp=1)
```

*Assess clustering pattern*

``` r
MCMC(coordinates=coord,
     geno.dip.codom=geno,
     varnpop=TRUE,
     npopmax=10,
     spatial=TRUE,
     freq.model="Uncorrelated",
     nit=1000000,
     thinning=1000,
     path.mcmc="./")
```

*Post-processing and plotting*

``` r
PostProcessChain(coordinates=coord,
                 geno.dip.codom=geno,
                 path.mcmc="./",
                 nxdom=100,
                 nydom=100,
                 burnin=200)
PostProcessChain(coordinates = coord,path.mcmc = "./",nxdom = 50,nydom = 50,burnin = 200)
PosteriorMode(coordinates=coord,
              path.mcmc="./",
              file="map.pdf")
Plotnpop(path.mcmc="./",
         burnin=500)
```

**END**

**References**

-   \[Geneland tutorial\]
    (<http://www.peterbeerli.com/classes/index.php?title=GenelandTutorial>)
