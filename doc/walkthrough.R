## ---- message=FALSE, warning=FALSE--------------------------------------------
library(conos)
library(dplyr)

## -----------------------------------------------------------------------------
install.packages('conosPanel', repos='https://kharchenkolab.github.io/drat/', type='source')

## -----------------------------------------------------------------------------
panel <- conosPanel::panel

## -----------------------------------------------------------------------------
str(panel, 1)

## -----------------------------------------------------------------------------
head(colnames(panel[[1]]))

## -----------------------------------------------------------------------------
any(duplicated(unlist(lapply(panel,colnames))))

## -----------------------------------------------------------------------------
install.packages('p2data', repos='https://kharchenkolab.github.io/drat/', type='source')

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(pagoda2)
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=1, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)

## -----------------------------------------------------------------------------
typeof(panel.preprocessed)
names(panel.preprocessed)

## ---- eval=FALSE--------------------------------------------------------------
#  library(Seurat)
#  panel.preprocessed.seurat <- lapply(panel, basicSeuratProc)

## -----------------------------------------------------------------------------
con <- Conos$new(panel.preprocessed, n.cores=1)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)

## ---- fig.height=4, fig.width=6-----------------------------------------------
plotComponentVariance(con, space='PCA')

## ---- fig.height=6, fig.width=6-----------------------------------------------
con$plotPanel(clustering="multilevel", use.local.clusters=TRUE, title.size=6)

## -----------------------------------------------------------------------------
con$findCommunities(method=leiden.community, resolution=1)

## ---- fig.height=6, fig.width=6-----------------------------------------------
con$plotPanel(font.size=4)

## ---- fig.height=6, fig.width=6-----------------------------------------------
plotClusterBarplots(con, legend.height = 0.1)

## ---- fig.height=6, fig.width=6-----------------------------------------------
con$plotPanel(gene = 'GZMK')

## ---- fig.height=6, fig.width=6, message=FALSE, warning=FALSE-----------------
## create new variable for panel.preprocessed
preprocessed_panel_example <- panel.preprocessed
## set NAs within all cell sin MantonCB2_HiSeq_1
preprocessed_panel_example$MantonCB2_HiSeq_1$clusters$PCA$multilevel <- NA
## create new Conos object
con_example <- Conos$new(panel.preprocessed, n.cores=1)
## construct joint graph
con_example$buildGraph()
## now plot the panel with NA values as "X"
con_example$plotPanel(clustering="multilevel", use.local.clusters=TRUE, title.size=6)

## -----------------------------------------------------------------------------
con$embedGraph(method='largeVis')

## ---- fig.height=6, fig.width=6, message=FALSE, warning=FALSE-----------------
con$plotGraph(alpha=0.1)

## ---- fig.height=4, fig.width=6-----------------------------------------------
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
con$embedGraph(alpha=0.001, embedding.name="example_embedding", sgd_batched=1e8)  

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(clustering='walktrap', size=0.1)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$embedGraph(method="UMAP", min.dist=0.01, spread=15, min.prob.lower=1e-3)

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(clustering='walktrap', size=0.1)

## ---- fig.width=6, fig.height=6, message=FALSE, warning=FALSE-----------------
con$plotPanel(clustering='walktrap', size=0.1, use.common.embedding=TRUE)

## -----------------------------------------------------------------------------
fc <- greedyModularityCut(con$clusters$walktrap$result, 40)

## ---- fig.height=6, fig.width=6-----------------------------------------------
con$plotGraph(groups=fc$groups, size=0.1)

## ---- fig.width=6, fig.height=4-----------------------------------------------
# fc$hc is an hclust structure ... here we will convert it to a dendrogram
dend <- as.dendrogram(fc$hc)
plot(dend)

## -----------------------------------------------------------------------------
cellannot <- read.table(file.path(find.package('conos'), 'extdata', 'cellannot.txt'), header=FALSE, sep='\t')
cellannot <- setNames(cellannot[,2], cellannot[,1])

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotPanel(groups = cellannot)

## -----------------------------------------------------------------------------
new.label.info <- con$propagateLabels(labels = cellannot, verbose=TRUE)

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotPanel(colors=new.label.info$uncertainty, show.legend=TRUE, legend.title="Uncertainty", legend.pos=c(1, 0))
con$plotPanel(groups=new.label.info$labels, show.legend=FALSE)

## -----------------------------------------------------------------------------
head(new.label.info$label.distribution)

## ---- message=FALSE, warning=FALSE--------------------------------------------
new.annot <- new.label.info$labels
de.info <- con$getDifferentialGenes(groups=new.annot, append.auc=TRUE)

## -----------------------------------------------------------------------------
head(de.info$`B cells`)

## ---- fig.width=6, fig.height=4-----------------------------------------------
cowplot::plot_grid(con$plotGraph(groups=new.annot), con$plotGraph(gene="CD74"))

## -----------------------------------------------------------------------------
de.info$monocytes %>% filter(AUC > 0.75) %>% arrange(-Precision) %>% head()

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(gene="CD14")

## ----fig.width=6, fig.height=6------------------------------------------------
plotDEheatmap(con,as.factor(new.annot),de.info, n.genes.per.cluster = 5, column.metadata=list(samples=con$getDatasetPerCell()), row.label.font.size = 7)

## ----fig.width=6, fig.height=6------------------------------------------------
gns <- c("GZMB","IL32","CD3E","LYZ","HLA-DRA","IGHD","GNLY","IGHM","GZMK")
plotDEheatmap(con,new.annot,de.info[-c(3,10)], n.genes.per.cluster = 30, column.metadata=list(samples=con$getDatasetPerCell()), row.label.font.size = 7, labeled.gene.subset = gns)

## -----------------------------------------------------------------------------
str(con$getClusterCountMatrices(), 1)

## -----------------------------------------------------------------------------
samplegroups <- list(
  bm = c("MantonBM1_HiSeq_1","MantonBM2_HiSeq_1"),
  cb = c("MantonCB1_HiSeq_1","MantonCB2_HiSeq_1")
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
de.info <- getPerCellTypeDE(con, groups=as.factor(new.annot), sample.groups = samplegroups, ref.level='bm', n.cores=1)

## -----------------------------------------------------------------------------
str(de.info[1:3], 2)

## -----------------------------------------------------------------------------
res <- de.info[['B cells']]$res
head(res[order(res$padj,decreasing = FALSE),])

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$buildGraph(k=15, k.self=5, alignment.strength=0.3, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)

## ---- fig.height=4, fig.width=6, message=FALSE, warning=FALSE-----------------
con$embedGraph(embedding.name="new_embedding")
con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)

## ---- fig.height=6, fig.width=6-----------------------------------------------
con$findCommunities()
plotClusterBarplots(con, legend.height = 0.1)

## -----------------------------------------------------------------------------
# library(pagoda2)
# p2app = p2app4conos(conos=con, file="conosApp1.bin", save=TRUE)
# show.app(app=p2app, name='conos_app')

