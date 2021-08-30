Discriminant Analysis of Principal Components (DAPC)
================
Achyut K Banerjee
August 23, 2021

**Premise:** To identify genetic clusters and describe the relationships
between them in a mangrove associate species Scaevola taccada in the
Indo-West Pacific (populations examined = 42, number of individuals =
744, SSR data).

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/mol.bio/dapc\]

*Install and load package*

``` r
library(adegenet)
library("poppr")
```

*Read and format data*

``` r
data_com <- read.fstat("C:/Users/Achyut/Desktop/web/mol.bio/dapc/data.dat", quiet = TRUE)
data_com
```

*Find clusters*

``` r
grp <- find.clusters(data_com, max.n.clust=42) 
#maximum number of clusters should be equal to number of populations# 
#Choose the number of PCs = maximum (say, 200)#
#number of clusters from where the BIC value keeps increasing, here 10#

find.clusters
names(grp)
head(grp$Kstat, 8)
head(grp$grp, 10)
grp$size
table(pop(data_com), grp$grp)
table.value(table(pop(data_com), grp$grp), col.lab=paste("inf", 1:3),
            row.lab=paste("Pop", 1:42))
table.value(table(pop(data_com), grp$grp))
```

*Selection of number of PCA axes*

``` r
genalex<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/dapc/data_com_genalex.csv")
data<-read.genalex("C:/Users/Achyut/Desktop/web/mol.bio/dapc/data_com_genalex.csv", 
                   ploidy = 2, geo = FALSE, region = FALSE, 
                   genclone = TRUE, sep = ",", recode = FALSE)
data
set.seed(999)
datax <- xvalDapc(tab(data, NA.method = "mean"), pop(data))
set.seed(999)
system.time(datax <- xvalDapc(tab(data, NA.method = "mean"), pop(data),
                              n.pca = 40:70, n.rep = 1000,
                              parallel = "multicore", ncpus = 4L))
names(datax)
datax[-1]
scatter(datax$DAPC)
class(datax$DAPC$posterior)
dim(datax$DAPC$posterior)
round(head(datax$DAPC$posterior),3)
summary(datax$DAPC)
```

*Compute DAPC*

``` r
dapc1 <- dapc(data_com, grp$grp,var.contrib = TRUE, scale = FALSE, n.pca = 64)
dapc1
scatter(dapc1)
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)

myCol <- c("darkblue","purple","green","orange","red","blue","chocolate","firebrick","forestgreen","red")
scatter(dapc1, posi.da="bottomleft", bg="white",
        pch=15:17, cstar=0, col=myCol, scree.pca=TRUE,
        posi.pca="bottomright")


scatter(dapc1,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)


scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:10))

class(dapc1$posterior)
dim(dapc1$posterior)
round(head(dapc1$posterior),3)
summary(dapc1)
assignplot(dapc1, subset=1:50)
compoplot(dapc1, posi="bottomright",
          txt.leg=paste("Cluster", 1:10), lab="",
          ncol=1, xlab="individuals", col=myCol)

temp <- which(apply(dapc1$posterior,1, function(e) all(e<0.9)))
temp
compoplot(dapc1, subset=temp, posi="bottomright",
          txt.leg=paste("Cluster", 1:10),
          ncol=10, col=myCol)

pop(data_com) <- data_com$pop
dapc2 <- dapc(data_com, var.contrib = TRUE, scale = FALSE, n.pca = 24, n.da = nPop(data_com) - 1)
scatter(dapc2, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2)
```

**END**

**References**

-   \[DAPC tutorial\]
    (<http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf>)
