## ---- message=FALSE, warning=FALSE--------------------------------------------
install.packages('p2data', repos='https://kharchenkolab.github.io/drat/', type='source')
install.packages('conosPanel', repos='https://kharchenkolab.github.io/drat/', type='source')

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(conos)
library(pagoda2)
library(parallel)
library(ggplot2)

## ---- message=FALSE, warning=FALSE--------------------------------------------
load(url("http://pklab.med.harvard.edu/peterk/conos/atac_rna/data.RData"))

## -----------------------------------------------------------------------------
str(data,1)

## ---- message=FALSE, warning=FALSE--------------------------------------------
p2l <- lapply(data,basicP2proc,n.odgenes=3e3,min.cells.per.gene=-1,nPcs=30,make.geneknn=FALSE,n.cores=1)

## ---- message=FALSE, warning=FALSE--------------------------------------------
## instantiate Conos object
con <- Conos$new(p2l, n.cores=1)

## build joint graph
con$buildGraph(k=15, k.self=5, k.self.weigh=0.01, ncomps=30, n.odgenes=5e3, space='PCA') 

## find communities
con$findCommunities(resolution=1.5)

## generate embedding
con$embedGraph(alpha=1/2)

## ---- fig.height=8, fig.width=8, message=FALSE, warning=FALSE-----------------
p1 <- con$plotGraph(font.size=c(3,5),title='conos clusters',alpha=0.2) #+ annotate("text", x=-Inf, y = Inf, label = "clusters", vjust=1, hjust=0)

p2 <- con$plotGraph(groups=rna.annotation,mark.groups=TRUE,alpha=0.2,plot.na=FALSE,title='annotation: RNA',font.size=c(3,5))+xlim(range(con$embedding[,1]))+ylim(range(con$embedding[,2]))

p2c <- con$plotGraph(groups=atac.annotation,mark.groups=T,alpha=0.2,plot.na=FALSE,title='annotation: ATAC',font.size=c(3,5))+xlim(range(con$embedding[,1]))+ylim(range(con$embedding[,2]))

p3 <- con$plotGraph(color.by='sample',mark.groups=F,alpha=0.1,show.legend=TRUE,title='platform',raster=TRUE)+theme(legend.position=c(1,1),legend.justification = c(1,1))+guides(color=guide_legend(ncol=2,override.aes = list(size=3,alpha=0.8)))

pp <- cowplot::plot_grid(plotlist=list(p2,p2c,p1,p3), ncol=2) 
print(pp)

