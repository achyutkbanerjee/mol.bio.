Relationship of genetic differentiation with IBD, IBE and IBR (Nonlinear
regression)
================
Achyut K Banerjee
August 23, 2021

**Premise:** To assess the relationship of genetic differentiation with
isolation by distance, environment and resistance surfaces for Mikania
micrantha (populations examined = 46, SSR data).

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/mol.bio/gdm\]

*Install and load package*

``` r
library(gdm)
```

*Read and format genetic distance data*

``` r
gen_dist<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/gdm/gen_dist.csv") ####(FST/1-FST) genetic distance without ENA#####
gen_dist_mat<-as.matrix(gen_dist)
biodata<-gen_dist_mat
```

*Preparing environmental distance data*

``` r
#Install and load packages#
library(biomod2)
library(ggplot2)
library(gridExtra)
library(rgdal)
library(raster)
library(ade4)

#Making stack of the bioclimatic variables#
bioclim_world<-stack(list.files("C:/Users/Achyut/Desktop/web/mol.bio/gdm/climate/current/bio",pattern="bio_",full.names = T),RAT=FALSE)

#Loading occurrence data#
data<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/gdm/occ_points.csv")
occ<-data[data$Site,]

#Creating masks#
mask_back<-shapefile("C:/Users/Achyut/Desktop/web/mol.bio/gdm/shp/invasive_kop_dissolve.shp")
bioclim_back<-mask(bioclim_world,mask_back)

#To obtain the identifiers where the species occurs#
points_occ<-data.frame(occ[1:46,c("LONG","LAT")])
occ_cell_id<-cellFromXY(subset(bioclim_world,1),points_occ)
bioclim_occ<-extract(bioclim_world,points_occ)
bioclim_occ_df<-na.omit(as.data.frame(bioclim_occ))
bioclim_occ_df<-as.data.frame(bioclim_occ)

####PCA####
library(bindrcpp)
library(factoextra)
library(rlang)
library(FactoMineR)
library(corrplot)
#Converting raster to data frame#
bioclim_back_df<-na.omit(as.data.frame(bioclim_back))
head(bioclim_back_df)

res.pca<- PCA(bioclim_back_df, scale.unit = TRUE, ncp = 3, graph = TRUE)

eig.val <- get_eigenvalue(res.pca)
eig.val
var <- get_pca_var(res.pca)
var
print(var)
head(var$coord)
head(var$cos2)
head(var$contrib,19)
head(var$coord, 19)
contrib<-data.frame(var$contrib)
coord<-data.frame(var$coord)
fviz_pca_var(res.pca, col.var = "black")

#Predict values of the locations across three PC axes#
predict_occ<-predict(res.pca,bioclim_occ_df)
predict_occ_df<-as.data.frame(predict_occ)
write.csv(predict_occ_df,"predict_occ3.csv")

#Make a file 'site.csv' with the coord.Dim.1-coord.Dim.3 values#
#Read the environmental data#
expdata<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/gdm/site.csv")
spdata<-expdata[,c(1,2,6,7)]
envdata<-expdata[,c(2:7)]
site<-unique(spdata$site)
dissim<-cbind(site,biodata)
#Format the environmental data#
exformat<-formatsitepair(dissim,3,XColumn = "long",YColumn = "lat",
                         predData = envdata,siteColumn = "site")
```

*Preparing geographic distance data*

``` r
envdata1<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/gdm/geo_dist.csv") #distance as calculated in GenAlex as Log(1+GGD)#
geo_dist_mat<-as.matrix(envdata1)
envdata1<-geo_dist_mat
envsim1<-cbind(site,envdata1)
```

*Preparing altitude distance data*

``` r
envdata3<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/gdm/alt_dist.csv") #distance generated from the circuitscape#
alt_dist_mat<-as.matrix(envdata3)
envdata3<-alt_dist_mat
envsim3<-cbind(site,envdata3)
```

*Preparing anthropogenic distance data*

``` r
envdata5<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/gdm/anthro_dist.csv") #distance generated from the circuitscape#
anthro_dist_mat<-as.matrix(envdata5)
envdata5<-anthro_dist_mat
envsim5<-cbind(site,envdata5)
```

*Data formatting for the models*

``` r
#Full model#
exformat2a<-formatsitepair(exformat,4,
                          predData = envdata,siteColumn = "site",
                          distPreds = list(as.matrix(envsim1),
                                           as.matrix(envsim3),as.matrix(envsim5)))
#Env+Geo#
exformat2b<-formatsitepair(exformat,4,
                           predData = envdata,siteColumn = "site",
                           distPreds = list(as.matrix(envsim1)))

#Env+Alt#
exformat2c2<-formatsitepair(exformat,4,
                            predData = envdata,siteColumn = "site",
                            distPreds = list(as.matrix(envsim3)))

#Env+Anthro#
exformat2c3<-formatsitepair(exformat,4,
                            predData = envdata,siteColumn = "site",
                            distPreds = list(as.matrix(envsim5)))

#Env only#
exformat2d<-formatsitepair(exformat,4,
                           predData = envdata,siteColumn = "site")


#Env+Resistance(Alt+Anthro)#
exformat2e<-formatsitepair(exformat,4,
                            predData = envdata,siteColumn = "site",
                            distPreds = list(as.matrix(envsim3),as.matrix(envsim5)))
```

*Calculate variable importance*

``` r
varimp<-gdm.varImp(exformat2a, geo=FALSE, splines = NULL, knots = NULL, 
                   fullModelOnly = FALSE,
                   nPerm = 500, parallel = FALSE, cores = 2, 
                   sampleSites = 1, sampleSitePairs = 1,
                   outFile = NULL)
##NOTE: Repeat for "exformat2b,c,d,e"

barplot(sort(varimp[[2]][,1], decreasing=T))
varimp[1] #Full model evaluation#
varimp[2] #Individual variable importance#
varimp[3] #Individual variable's p-value#
```

*Data formatting for the geographic models*

``` r
expdata1<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/gdm/site_1.csv")#all environmental PCs with same weight#
spdata1<-expdata1[,c(1,2,6,7)]
envdata6<-expdata1[,c(2:7)]
site<-unique(spdata1$site)
dissim<-cbind(site,biodata)
exformat2e<-formatsitepair(dissim,3,XColumn = "long",YColumn = "lat",
                           predData = envdata6,siteColumn = "site")
exformat2f<-formatsitepair(exformat2e,4,
                           predData = envdata6,siteColumn = "site",
                           distPreds = list(as.matrix(envsim3),as.matrix(envsim5)))
exformat2g<-exformat2f[,c(1:6,10:11,15:16)] #Resistance only#
exformat2j<-formatsitepair(exformat2e,4,
                           predData = envdata6,siteColumn = "site",
                           distPreds = list(as.matrix(envsim1),as.matrix(envsim3),as.matrix(envsim5)))
exformat2k<-exformat2j[,c(1:6,10:12,16:18)]#Geo+Resistance#

exformat2h<-formatsitepair(exformat2e,4,
                           predData = envdata5,siteColumn = "site",
                           distPreds = list(as.matrix(envsim1)))
exformat2i<-exformat2h[,c(1:6,10,14)] #Geo only#
```

*Calculate variable importance*

``` r
varimp<-gdm.varImp(exformat2j, geo=FALSE, splines = NULL, knots = NULL, fullModelOnly = FALSE,
                   nPerm = 500, parallel = FALSE, cores = 2, sampleSites = 1, sampleSitePairs = 1,
                   outFile = NULL)
##NOTE: Repeat for "exformat2k"

barplot(sort(varimp[[2]][,1], decreasing=T))
varimp[1] #Full model evaluation#
varimp[2] #Individual variable importance#
varimp[3] #Individual variable's p-value#


gdm<-gdm(exformat2i, geo=FALSE, splines=NULL, knots=NULL) #variable importance for 2i model#
summary.gdm(gdm)
gdm
```

**END**

**References** \* \[Paper - Annals in Botany\]
(<https://doi.org/10.1093/aob/mcaa044>) \* \[gdm package manual\]
(<https://cran.microsoft.com/snapshot/2015-07-22/web/packages/gdm/gdm.pdf>)
