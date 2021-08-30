Standardized measure of genetic differentiation
================
Achyut K Banerjee
August 23, 2021

**Premise:** To estimate population differentiation based on cpDNA
sequence data, Hedrick’s G/ST has been estimated for two mangrove
species of genus Lumnitzera - L. littorea (LL) and L. racemosa (LR)
(populations examined = 27 for LL, 32 for LR).

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/mol.bio/gen.diversity\]

*Install and load package*

``` r
library("apex")
library("adegenet")
library("pegas")
library("mmod")
library("poppr")
```

*Read data and estimate parameters for species\_1*

``` r
data<-read.multiFASTA("C:/Users/Achyut/Desktop/web/mol.bio/gen.diversity/data_LL.fas")
plot(data, cex = 0.2)
getLocusNames(data)
(setLocusNames(data) <- gsub(".fas", "", getLocusNames(data)))
data.gid <- multidna2genind(data, mlst = TRUE)
data.gid
my_strata <- data.frame(populations = rep(c("BAS","MDI","RNT","KPT","PNT","LKW","TJM","SJM","SCM","KNT","CPT","TYT","CTT","KKC","KCM","KPM","SDM","TWM","CRP","PAP","BPP","TLC","IBP","BUI","SRI","DRA","EVI"),
                                       c(5,8,16,18,7,10,8,12,13,12,6,8,13,8,23,18,15,9,13,10,16,13,6,16,18,12,16)))
strata(data.gid) <- my_strata                        
setPop(data.gid) <- ~populations
data.gid
diff_stats(data.gid)
Phi_st_Meirmans(data.gid)
pairwise_Gst_Nei(data.gid, linearized = FALSE)
Hedrick_Gst<-pairwise_Gst_Hedrick(data.gid, linearized = FALSE)
Hedrick_Gst_matrix<-as.matrix(dist(Hedrick_Gst))
write.csv(Hedrick_Gst_matrix,"Hedrick_Gst_LL.csv")
```

*Read data and estimate parameters for species\_2*

``` r
data<-read.multiFASTA("C:/Users/Achyut/Desktop/web/mol.bio/gen.diversity/data_LR.fas")
plot(data, cex = 0.2)
getLocusNames(data)
(setLocusNames(data) <- gsub(".fas", "", getLocusNames(data)))
data.gid <- multidna2genind(data, mlst = TRUE)
data.gid
my_strata <- data.frame(populations = rep(c("MIK","RSL","PTS","BSL","RNT","LKW","KSM","TJM","SGP","KNT","TNT","STT","CPT","KKT","BKT","CTT","KKC","KEC","KPM","LKM","SDM","KBI","SMC","SCC","FGC","TMC","YHC","IBP","KLP","SRI","LCA","CQA"),
                                          c(9,12,12,14,3,18,8,14,9,8,13,15,7,8,12,14,18,11,16,9,12,14,13,12,14,13,13,8,6,10,27,13)))
strata(data.gid) <- my_strata                        
setPop(data.gid) <- ~populations
data.gid
diff_stats(data.gid)
Phi_st_Meirmans(data.gid)
pairwise_Gst_Nei(data.gid, linearized = FALSE)
Hedrick_Gst<-pairwise_Gst_Hedrick(data.gid, linearized = FALSE)
Hedrick_Gst_matrix<-as.matrix(dist(Hedrick_Gst))
write.csv(Hedrick_Gst_matrix,"Hedrick_Gst_LR.csv")
```

**END**

**References**

-   \[Paper - Frontiers in Plant Science\]
    (<https://www.frontiersin.org/articles/10.3389/fpls.2021.637009/full>)
