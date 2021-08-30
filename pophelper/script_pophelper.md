Plotting genetic structure
================
Achyut K Banerjee
August 23, 2021

**Premise:** To plot the output of BAPS admixture model. The data is
generated from a phylogeographic study of a mangrove associate species
Derris trifoliata in the Indo-West Pacific, having 356 individuals
distributed among 30 populations.

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/mol.bio/pophelper\]

*Install and load package*

``` r
library(pophelper)
library(gridExtra)
library(grid)
```

*Read data*

``` r
sfiles <- "C:/Users/Achyut/Desktop/web/mol.bio/pophelper/data_baps.mat.txt"
readQ(files=sfiles,filetype="baps")
slist <- readQ(files=sfiles,indlabfromfile=T)

##Load individual population name##
inds <- read.delim("C:/Users/Achyut/Desktop/web/mol.bio/pophelper/label.txt",header=FALSE,stringsAsFactors = FALSE)
rownames(slist[[1]]) <- inds$V1 #attach to the data#

##Load predefined group names##
group<-read.delim("C:/Users/Achyut/Desktop/web/mol.bio/pophelper/group.txt",header=TRUE,stringsAsFactors = FALSE)
#Making subsets#
onelabset <- group[,2,drop=F]
twolabset<-group[,2:3]
twolabset1<-group[,1:2]
```

*Add color palette (optional)*

``` r
clist <- list(
  "shiny"=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
  "strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"),
  "oceanfive"=c("#00A0B0", "#6A4A3C", "#CC333F", "#EB6841", "#EDC951"),
  "keeled"=c("#48B098", "#91CB62", "#FFEE3B", "#FB9013", "#FF3C28"),
  "vintage"=c("#400F13", "#027368", "#A3BF3F", "#F2B950", "#D93A2B"),
  "muted"=c("#46BDDD","#82DDCE","#F5F06A","#F5CC6A","#F57E6A"),
  "teal"=c("#CFF09E","#A8DBA8","#79BD9A","#3B8686","#0B486B"),
  "merry"=c("#5BC0EB","#FDE74C","#9BC53D","#E55934","#FA7921"),
  "funky"=c("#A6CEE3", "#3F8EAA", "#79C360", "#E52829", "#FDB762","#ED8F47","#9471B4"),
  "retro"=c("#01948E","#A9C4E2","#E23560","#01A7B3","#FDA963","#323665","#EC687D"),
  "cb_paired"=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
  "cb_set3"=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
  "morris"=c("#4D94CC","#34648A","#8B658A","#9ACD32","#CC95CC","#9ACD32","#8B3A39","#CD6601","#CC5C5B","#8A4500"),
  "wong"=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#006699","#D55E00","#CC79A7"),
  "krzywinski"=c("#006E82","#8214A0","#005AC8","#00A0FA","#FA78FA","#14D2DC","#AA0A3C","#FA7850","#0AB45A","#F0F032","#A0FA82","#FAE6BE"))

# add length of palettes
lengths <- sapply(clist,length)
names(clist) <- paste0(names(clist),"_",lengths)

par(mar=c(0.2,6,0.2,0))
par(mfrow=c(length(clist),1))

for(i in 1:length(clist))
{
  {barplot(rep(1,max(lengths)),col=c(clist[[i]],rep("white",max(lengths)-length(clist[[i]]))),axes=F,border=F)
    text(x=-0.1,y=0.5,adj=1,label=names(clist)[i],xpd=T,cex=1.2)}
}
```

*Plot*

``` r
#export plot = TRUE#
plotQMultiline(slist[1],barbordercolour="black",barsize = 0.8,lpp = 6,
                     spl=90,exportplot=T,returnplot=F,grplab=twolabset,ordergrp=TRUE,selgrp="region",
                     useindlab=T,indlabcol="black",showlegend=T,clustercol=clist$shiny,grplabcol = "black",
                     outputfilename="plotq3",imgtype="png",exportpath=getwd())
#export plot = FALSE#
plotQMultiline(slist[1],barbordercolour="black",barsize = 0.8,lpp = 6,
                     spl=90,exportplot=F,returnplot=F,grplab=twolabset,ordergrp=TRUE,selgrp="region",
                     useindlab=T,indlabcol="black",showlegend=T,clustercol=clist$shiny,grplabcol = "black")#export plot TRUE#
```

**END**

**References**

-   \[POPHELPER tutorial\]
    (<http://www.royfrancis.com/pophelper/index.html>)
