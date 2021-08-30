Relationship of genetic differentiation with IBD, IBE and IBR (Linear
regression)
================
Achyut K Banerjee
August 23, 2021

**Premise:** To assess the relationship of genetic differentiation with
isolation by distance, environment and resistance surfaces for Mikania
micrantha (populations examined = 46, SSR data).

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/mol.bio/mrm\]

*Install and load package*

``` r
library(ecodist)
```

*Read and format genetic distance data*

``` r
gen_dist<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/mrm/gen_dist.csv")#(FST/1-FST) genetic distance without ENA#
gen_dist
gen_dist<-as.dist(gen_dist)
gen_dist
```

*Read and format environmental distance data*

``` r
data=read.csv("predict_occ3_PC.csv")#Read environmental distance data#
#The environmental data has been generated following the same procedure as gdm. Refer to the process#

#Computing euclidean distances between points for 3 PCs, together first and then individually#
euc_envdist<-distance(data[,2:4], "euclidean")
data1=read.csv("predict_occ3_PC1.csv")
euc_envdist1<-distance(data1[,2], "euclidean")
data2=read.csv("predict_occ3_PC2.csv")
euc_envdist2<-distance(data2[,2], "euclidean")
data3=read.csv("predict_occ3_PC3.csv")
euc_envdist3<-distance(data3[,2], "euclidean")
euc_envdist
```

*Preparing geographic distance data*

``` r
geo_dist<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/mrm/geo_dist.csv")#distance as calculated in GenAlex in Log(1+GGD)#
geo_dist
geo_dist<-as.dist(geo_dist)
geo_dist
```

*Preparing altitude distance data*

``` r
alt_dist<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/mrm/alt_dist.csv")#distance generated from the circuitscape#
alt_dist<-as.dist(alt_dist)
```

*Preparing anthropogenic distance data*

``` r
anthro_dist<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/mrm/anthro_dist.csv")#distance generated from the circuitscape#
anthro_dist<-as.dist(anthro_dist)
```

*Fitting models*

``` r
#Full model#
MRM(gen_dist ~ geo_dist + euc_envdist1+euc_envdist2+euc_envdist3+
      alt_dist+anthro_dist, nperm = 10000, mrank = TRUE)
MRM(gen_dist ~ geo_dist + euc_envdist+alt_dist+anthro_dist, nperm = 10000, mrank = TRUE)

#Geo only#
MRM(gen_dist ~ geo_dist, nperm = 10000, mrank = TRUE)
#Env only#
MRM(gen_dist ~ euc_envdist, nperm = 10000, mrank = TRUE)
MRM(gen_dist ~ euc_envdist1+euc_envdist2+euc_envdist3, nperm = 10000, mrank = TRUE)
#Alt+Anthro#
MRM(gen_dist ~ alt_dist+anthro_dist, nperm = 10000, mrank = TRUE)
#Geo+Env#
MRM(gen_dist ~ geo_dist + euc_envdist1+euc_envdist2+euc_envdist3, nperm = 10000, mrank = TRUE)
#Geo+Alt#
MRM(gen_dist ~ geo_dist + alt_dist, nperm = 10000, mrank = TRUE)
#Geo+Anthro#
MRM(gen_dist ~ geo_dist + anthro_dist, nperm = 10000, mrank = TRUE)
#Alt+Env#
MRM(gen_dist ~ euc_envdist1+euc_envdist2+euc_envdist3+ alt_dist, nperm = 10000, mrank = TRUE)
#Env+Anthro#
MRM(gen_dist ~ euc_envdist1+euc_envdist2+euc_envdist3+ anthro_dist, nperm = 10000, mrank = TRUE)
#Alt only#
MRM(gen_dist ~ alt_dist, nperm = 10000, mrank = TRUE)
#Anthro only#
MRM(gen_dist ~ anthro_dist, nperm = 10000, mrank = TRUE)
#Geo+Alt+Anthro#
MRM(gen_dist ~ geo_dist + alt_dist+anthro_dist, nperm = 10000, mrank = TRUE)
#Env+Alt+Anthro#
MRM(gen_dist ~ euc_envdist1+euc_envdist2+euc_envdist3+alt_dist+anthro_dist, nperm = 10000, mrank = TRUE)
```

**END**

**References** \* \[Paper - Annals in Botany\]
(<https://doi.org/10.1093/aob/mcaa044>)
