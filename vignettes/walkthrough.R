## ---- message=FALSE, warning=FALSE--------------------------------------------
library(conos)
library(dplyr)

## -----------------------------------------------------------------------------
panel <- readRDS(file.path(find.package('conos'),'extdata','panel.rds'))

## -----------------------------------------------------------------------------
str(panel, 1)

## -----------------------------------------------------------------------------
head(colnames(panel[[1]]))

## -----------------------------------------------------------------------------
any(duplicated(unlist(lapply(panel,colnames))))

## -----------------------------------------------------------------------------
library(pagoda2)
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=4, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)

## -----------------------------------------------------------------------------
str(panel.preprocessed, 1)

## ---- eval=FALSE--------------------------------------------------------------
#  library(Seurat)
#  panel.preprocessed <- lapply(panel, basicSeuratProc)

## -----------------------------------------------------------------------------
con <- Conos$new(panel.preprocessed, n.cores=4)

## -----------------------------------------------------------------------------
str(con$samples,1)

## ---- fig.height=8, fig.width=8-----------------------------------------------
con$plotPanel(clustering="multilevel", use.local.clusters=TRUE, title.size=6)

## -----------------------------------------------------------------------------
con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)

## ---- fig.height=6, fig.width=8-----------------------------------------------
plotComponentVariance(con, space='PCA')

## -----------------------------------------------------------------------------
con$findCommunities(method=leiden.community, resolution=1)

## ---- fig.height=8, fig.width=8-----------------------------------------------
con$plotPanel(font.size=4)

## ---- fig.height=8, fig.width=8-----------------------------------------------
plotClusterBarplots(con, legend.height = 0.1)

## ---- fig.height=8, fig.width=8-----------------------------------------------
con$plotPanel(gene = 'GZMK')

## ---- fig.height=6, fig.width=6, message=FALSE, warning=FALSE-----------------
con$plotGraph(alpha=0.1)

## ---- fig.height=6, fig.width=8-----------------------------------------------
con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)

## ---- fig.height=6, fig.width=6-----------------------------------------------
con$plotGraph(gene='GZMK', title='GZMK expression')

## -----------------------------------------------------------------------------
con$findCommunities(method = igraph::walktrap.community, steps=7)

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotPanel(clustering='walktrap', font.size=4)

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(clustering='walktrap')

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$embedGraph(alpha=0.001, sgd_batched=1e8)  

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(clustering='walktrap', size=0.1)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=4, min.prob.lower=1e-3)

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(clustering='walktrap', size=0.1)

## ---- fig.width=8, fig.height=8, message=FALSE, warning=FALSE-----------------
con$plotPanel(clustering='walktrap', size=0.1, use.common.embedding=TRUE)

## -----------------------------------------------------------------------------
fc <- greedyModularityCut(con$clusters$walktrap$result, 40)

## ---- fig.height=6, fig.width=6-----------------------------------------------
con$plotGraph(groups=fc$groups, size=0.1)

## ---- fig.width=8, fig.height=6-----------------------------------------------
# fc$hc is an hclust structure ... here we will convert it to a dendrogram
dend <- as.dendrogram(fc$hc)
plot(dend)

## -----------------------------------------------------------------------------
samf <- con$getDatasetPerCell()
str(samf)

## ---- fig.width=8, fig.height=6-----------------------------------------------
dend <- dendSetWidthByBreadth(dend, samf, fc$leafContent, min.width=1, max.width=4)
plot(dend)

## -----------------------------------------------------------------------------
tissue.factor <- as.factor(setNames(ifelse(grepl('BM',names(samf)),'BM','CB'), names(samf)))
str(tissue.factor)

## ---- fig.width=8, fig.height=6-----------------------------------------------
dend <- dendSetColorByMixture(dend, tissue.factor, fc$leafContent)
plot(dend)

## ---- eval=FALSE--------------------------------------------------------------
#  library(conosViz)
#  conosViz::conosShinyApp(con, N=30)

## -----------------------------------------------------------------------------
cellannot <- read.table(file.path(find.package('conos'),'extdata','cellannot.txt'),header=FALSE,sep='\t')
cellannot <- setNames(cellannot[,2], cellannot[,1])

## ---- fig.width=8, fig.height=8-----------------------------------------------
con$plotPanel(groups = cellannot)

## -----------------------------------------------------------------------------
new.label.info <- con$propagateLabels(labels = cellannot, verbose=TRUE)

## ---- fig.width=8, fig.height=8-----------------------------------------------
con$plotPanel(colors=new.label.info$uncertainty, show.legend=TRUE, legend.title="Uncertainty", legend.pos=c(1, 0))
con$plotPanel(groups=new.label.info$labels, show.legend=FALSE)

## -----------------------------------------------------------------------------
head(new.label.info$label.distribution)

## -----------------------------------------------------------------------------
new.annot <- new.label.info$labels
de.info <- con$getDifferentialGenes(groups=new.annot, n.cores=4, append.auc=TRUE)
head(de.info$`B cells`)

## ---- fig.width=8, fig.height=6-----------------------------------------------
cowplot::plot_grid(con$plotGraph(groups=new.annot), con$plotGraph(gene="CD74"))

## -----------------------------------------------------------------------------
de.info$monocytes %>% filter(AUC > 0.75) %>% arrange(-Precision) %>% head()

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(gene="CD14")

## ----fig.width=8, fig.height=8------------------------------------------------
plotDEheatmap(con,as.factor(new.annot),de.info, n.genes.per.cluster = 5, column.metadata=list(samples=con$getDatasetPerCell()), row.label.font.size = 7)

## ----fig.width=8, fig.height=8------------------------------------------------
gns <- c("GZMB","IL32","CD3E","LYZ","HLA-DRA","IGHD","GNLY","IGHM","GZMK")
plotDEheatmap(con,new.annot,de.info[-c(3,10)], n.genes.per.cluster = 30, column.metadata=list(samples=con$getDatasetPerCell()), row.label.font.size = 7, labeled.gene.subset = gns)

## -----------------------------------------------------------------------------
str(con$getClusterCountMatrices(), 1)

## -----------------------------------------------------------------------------
samplegroups <- list(
  bm = c("MantonBM1_HiSeq_1","MantonBM2_HiSeq_1"),
  cb = c("MantonCB1_HiSeq_1","MantonCB2_HiSeq_1")
)

## -----------------------------------------------------------------------------
de.info <- getPerCellTypeDE(con, groups=as.factor(new.annot), sample.groups = samplegroups, ref.level='bm', n.cores=4)

## -----------------------------------------------------------------------------
str(de.info[1:3], 2)

