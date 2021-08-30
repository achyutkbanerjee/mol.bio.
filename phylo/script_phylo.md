Phylogenetic tree
================
Achyut K Banerjee
August 22, 2021

**Premise:** To generate phylogenetic tree (species and genus level).
This data is a subset of the ILORA database version 1.0.

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/mol.bio/phylo\]

*Install and load package*

``` r
##Install V.PHYLOMAKER from GITHUB##
install.packages("remotes")
remotes::install_github("jinyizju/V.PhyloMaker")
install_github("jinyizju/V.PhyloMaker")

##Install V.PHYLOMAKER locally##
unzip("V.PhyloMaker-master.zip")
file.rename("V.PhyloMaker-master", "V.PhyloMaker")
shell("R CMD build V.PhyloMaker")
install.packages("V.PhyloMaker_0.1.0.tar.gz", repos = NULL)

##Load library##
library("V.PhyloMaker")
library("phytools")
```

*Read data and generate the tree*

``` r
data<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/phylo/species_list.csv")
tree.a<-phylo.maker(sp.list = data,tree = GBOTB.extended,nodes = nodes.info.1,scenarios = "S3")
##Write the tree##
write.tree(tree.a$scenario.3,"tree1.tre")
```

*Create genus level tree*

``` r
tree<-read.tree("tree1.tre")
tips<-tree$tip.label
genera<-unique(sapply(strsplit(tips,"_"),function(x) x[1]))
genera
ii<-sapply(genera,function(x,y) grep(x,y)[1],y=tips)
tree<-drop.tip(tree,setdiff(tree$tip.label,tips[ii]))
plotTree(tree,ftype="i")
tree$tip.label<-sapply(strsplit(tree$tip.label,"_"),function(x) x[1])
plotTree(tree,ftype="i")
write.tree(tree,"tree1.genus.tre")
```

**END**

**Premise\_2:** To check phylogenetic signals for the traits by
estimating Blomberg’s K and Pagel’s λ. For the categorical growth form
variable, we used the phylo.signal.disc() function to get an analogous
estimate of Blomberg’s K. This data is a subset of the ILORA database
version 1.0.

**Data availability:** The associated data (and the function) is
available in \[C:/Users/Achyut/Desktop/web/mol.bio/phylo\]

*For discrete trait* *Install and load package*

``` r
library(phytools)
```

*Assess phylogenetic signal*

``` r
tree<-read.tree("tree1.tre")
trait<-read.csv("C:/Users/Achyut/Desktop/web/mol.bio/phylo/species_data.csv")
trait<-as.data.frame(trait)
k.trait<-phylosig(tree,trait$Seed,test = TRUE)#Replace 'Seed' with other traits#
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree,trait$Seed,method = "lambda",test = TRUE)#Replace 'Seed' with other traits#
print(lambda.trait)
plot(lambda.trait)
```

*For categorical trait (growth form)* *Install and load packages*

``` r
library(ape)
library(geiger)
library(phangorn)
library(phylobase)
```

*Assess phylogenetic signal*

``` r
tree <-read.tree("tree1.tre")
trait<-read.csv("species_data_gf.csv")
trait<-as.data.frame(trait)
tmpTr <- drop.tip(tree, tree$tip.label[! tree$tip.label %in% trait[, 1]])
tr4 <- phylo4d(tmpTr, trait, label.type="column")
pruned <- as(extractTree(tr4), "phylo")

ldata <- tdata(tr4, "tip")
ldata <- ldata[pruned$tip.label, 1]
names(ldata) <- pruned$tip.label
phylo.signal.disc(ldata, pruned) ->result
write(result$.Randomization.Results, file="tree.gf.out")
```

**END**

**References**

-   \[Paper - Journal of Environmental Management\]
    (<https://doi.org/10.1016/j.jenvman.2021.113054>)
