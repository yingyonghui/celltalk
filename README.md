# celltalk
celltalk is an R package for inference and analysis of ligand-receptor interactions from single cell RNA sequencing data
## Installation
celltalk R package can be easily installed from Github using devtools:
```
devtools::install_github("yingyonghui/celltalk")
library(celltalk)
```
### Dependencies
- [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
- [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
- [GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html)

## Tutorials
#### Dependencies:

```
library(circlize)
library(reshape2)
library(dplyr)
library(ggplot2)
library(GSVA)
```
#### celltalk内置数据集，用于展示celltalk计算过程：

```
data("cellTalkSample",package='celltalk')
```
***sample.expr*** : expression matrix gene * cell, used as an example

***sample.lable*** : vector  of identity classes of cells in the expression matrix

***sample.marker*** : data frame of marker genes for each identity class, usually calculated by FindAllMarkers from Seurat

***gsva.mat*** : precomputed gsva scores for the example dataset

#### LR关系鉴定
```
# find significant LR pairs
marker.dat = sample.marker
species = 'mmusculus'
logFC.thre = 0.25
p.thre = 0.01
Interact <- findLRpairs(marker.dat=sample.marker, 
    species=species, 
    logFC.thre=logFC.thre, 
    p.thre=p.thre)

# type ?findLRpairs to get more information about each parameter
?findLRpairs
```
函数返回***Interact***为list对象，细胞之间交互数量信息存储在Interact[['InteractNumer']]，具体的交互基因存储在Interact[['InteractGeneUnfold']]

根据交互数量做circos plot：
```
circosPlot(Interact=Interact)
```
LR相互作用dotplot：
```
# present a dot plot of LR pairs for specific clusters
receptor.ident=6
dotPlot(sample.marker, Interact=Interact, receptor.ident=receptor.ident)
```

For specific upstream cells and their ligands, find the downstream cells and their receptors：
```
select.ident = 6
select.ligand = 'Dkk3'
ident.down.dat <- findReceptor(Interact=Interact, 
    select.ident=select.ident, 
    select.ligand=select.ligand)

head(ident.down.dat)
```
#### 通路分析
```
# find pathways in which genesets show overlap 
# with the ligands and receptors in the example dataset
Interact <- findLRpath(Interact=Interact, category='all')
```
Now genesets show overlap with the ligands and receptors in the exsample dataset are saved in Interact[['pathwayLR']]

gsva analysis：
```
# to compute gsva score by gsva function from the GSVA package
# to save time, we have precomputed gsva score and saved it in the varible *gsva.mat*
# gsva.mat <- gsva(sample.expr, Interact[['pathwayLR']], min.sz=10, parallel.sz=10)
```
gsva pathway heatmap, to display the highly variable pathways among all cells (待实现)


Pathway differential enrichment analysis：
```
# to find the different enriched pathways for cells in the selected identity class 
# and the receptor and ligand in the pathway
ident.lable = sample.lable
select.ident.1 = 6
test.res.dat <- diffLRpath(Interact=Interact, 
    gsva.mat=gsva.mat, 
    ident.lable=ident.lable, 
    select.ident.1=6,
    method='t.test')

head(test.res.dat)
```
Columns ***mean.diff***, ***mean.1***, ***mean.2***, ***t***, ***df***, ***p.val***, ***p.val.adj*** show the statistic result; *description* shows the name of pathway; 

Columns ***cell.up*** and ***ligand.up*** show the upstream identity classes which would release specific ligands to interact with the receptors from the current identity class; 

Column ***receptor.in.path*** shows the marker receptors expressed by the current identity class and these receptors are included in the current pathway;

Column ***ligand.in.path*** shows the marker ligands released by the current identity class and these ligands are also included in the current pathway.

Then we can find the second interactions mediated by specific pathways by identifying the downstream receptor of ligands in the ***ligand.in.path*** columns via findReceptor:
```
select.ident = 6
select.ligand = 'Dkk3'
ident.down.dat <- findReceptor(Interact=Interact, 
    select.ident=select.ident, 
    select.ligand=select.ligand)

head(ident.down.dat)
```