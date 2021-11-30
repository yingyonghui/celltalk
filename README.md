# celltalk
## Installation
celltalk R package can be easily installed from Github using devtools:
```
devtools::install_github("yingyonghui/celltalk")
```
### Dependencies 
- [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
- [GSVA](https://cran.r-project.org/web/packages/reshape2/index.html)

## Tutorials
Dependencies:

```
library(circlize)
library(reshape2)
library(dplyr)
library(ggplot2)
library(GSVA)
```
celltalk内置数据集，用于展示celltalk计算过程：

```
load('cellTalk.sample.RData')
```
sample.expr: expression matrix gene * cell, used as an example
sample.lable: vector  of identity classes of cells in the expression matrix
sample.marker: data frame of marker genes for each identity class, usually calculated by FindAllMarkers from Seurat
```
# find significant LR pairs
species = 'mmusculus'
logFC.thre = 0.25
p.thre = 0.01
Interact <- findLRpairs(sample.marker,species=species,logFC.thre=logFC.thre,p.thre=p.thre)

# type ?findLRpairs to get more information about each parameter
?findLRpairs
```
函数返回Interact为list对象，交互数量信息存储在 Interact[['InteractNumer']],具体的交互基因存储在Interact[['InteractGeneUnfold']]

根据交互数量做circos plot
```
circosPlot(Interact=Interact)
```
LR作用的点图
```
# present a dot plot of LR pairs for specific clusters
receptor.ident=6
dotPlot(all.marker.dat=sample.marker, Interact=Interact, receptor.ident=receptor.ident)
```

```
# find pathways in which genesets show overlap with the ligands and receptors in the exsample dataset
Interact <- findLRpath(Interact=Interact)
```
Now genesets show overlap with the ligands and receptors in the exsample dataset are saved in Interact[['pathwayLR']]

GSVA analysis
```
# computed gsva score by gsva function from the GSVA package
# to save time, we have precomputed gsva score and saved it in the varible *gsva.mat*
# gsva.mat <- gsva(sample.expr, Interact[['pathwayLR']], min.sz=10, verbose=FALSE, parallel.sz=10)
```

Pathway differential enrichment analysis 
```
# to find the different enriched pathways for cells in the selected identity class and the receptor and ligand in the pathway
ident.lable = sample.lable
select.ident.1 = 6
test.res.dat <- diffLRpath(Interact=Interact, gsva.mat=gsva.mat, ident.lable=ident.lable, select.ident.1=6)
head(test.res.dat)
```
Columns ***mean.diff***, ***mean.1***, ***mean.2***, ***t***, ***df***, ***p.val***, ***p.val.adj*** show the statistic result; *description* shows the name of pathway; 
Columns ***cell.up*** and ***ligand.up*** show the upstream identity classes which would release specific ligands to interact with the receptors from the current identity class; 
Column ***receptor.in.path*** shows the marker receptors expressed by the current identity class and these receptors are included in the current pathway;
Column ***ligand.in.path*** shows the marker ligands released by the current identity class and these ligands are also included in the current pathway.