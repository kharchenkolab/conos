# Default workflow

The following code presents a default workflow for working with multiple pagoda2 objects calling joint clusters
and generating differential expression based on them.

```r
## Make ensemble p2 object
p2ens <- clusterMatch:::Pagoda2ensemble$new()

## Set the apps
p2ens$setObjects(list.of.p2.objects) 

## Run 3 different joint clustering on the objects (these take some time to run)
jcl1 <- cpcaJCp3(p2ens$p2objs,k=20,return.details=T,n.odgenes = 2e3,ncomps=30,k.self=5,k.self.weight=0.1,community.detection.method = multilevel.community,n.cores=34,neighborhood.average=TRUE,reduction.method='CPCA',var.scale=T,neighborhood.average.k=10,xl=NULL)

jcl2 <- cpcaJCp3(p2ens$p2objs,k=20,return.details=T,n.odgenes = 2e3,ncomps=30,k.self=5,k.self.weight=0.1,community.detection.method = cluster_louvain,n.cores=34,neighborhood.average=TRUE,reduction.method='CPCA',var.scale=T,neighborhood.average.k=10,xl=jcl1$xl)

## When using walktrap we can prorocess
jcl3 <- cpcaJCp3(p2ens$p2objs,k=20,return.details=T,n.odgenes = 2e3,ncomps=30,k.self=5,k.self.weight=0.1,community.detection.method = walktrap.community,n.cores=34,neighborhood.average=TRUE,reduction.method='CPCA',var.scale=T,neighborhood.average.k=10,xl=jcl1$xl,steps=8)

## postprocess walktrap to select granularity
jcl3.cleanup <- postProcessWalktrapClusters(jcl3, p2ens, no.cl=200, size.cutoff=150,n.cores=32) 
levels(jcl3.cleanup) <- seq_along(levels(jcl3.cleanup))


## Plot all object in the set with specified clustering
p2ens$plotWithGroups(jcl1$groups) 
p2ens$plotWithGroups(jcl2$groups) 
p2ens$plotWithGroups(jcl3$groups) 
p2ens$plotWithGroups(jcl3.cleanup) 

## Plot of cell counts in each app
v <- melt(p2ens$getClusterCountsPerApp(groups=jcl3.cleanup.rename))
v$app <- factor(as.character(v$app),levels=levels(v$app)[order(strpart(levels(v$app),'-',2))])
ggplot(v, aes(x=app,y=cluster,fill=log10(value+1),label=value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradient(low='grey90',high='red') + geom_text()  
ggsave('cellcounts.png',width=16,height=9)

## Plot of cell percentages in each app
vc <- acast(v, cluster ~ app, value.var='value')                                                                          
v.prop <- melt(sweep(vc ,2,apply(vc,2,sum),FUN='/'))                                                    
colnames(v.prop) <- c('cluster','app','value')                                                         
ggplot(v.prop, aes(x=app,y=cluster,fill=value,label=scales::percent(value))) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_gradient(low='grey90',high='red') + geom_text(size=3) + ggtitle('Cell type abundance as fraction of sample') 
ggsave('cellpercentages.png',width=16,height=9)   

                                                 
## Prepare p2 ensembl for differential expression by summarising into clusters:sample combinations
p2ens$extractCountMatrices()  
p2ens$makeFullMatrix()
                                                                                    
## Setup metadata (to be automated)
cid <- rownames(p2ens$full.matrix)
sample <- strpart(cid, '_',1)
individual <- strpart(sample,'-',1) 
sample.type <- strpart(sample,'-',2)

row.meta <- data.frame(
    row.names=cid,
    sample=sample,
    individual=individual,
    sample.type=sample.type)                                                                                                                                                                                                                                                   
row.meta$cluster <- jcl3.cleanup.rename[cid]                                                                                                                        
## Set metadata in the object
p2ens$setCellMetadata(row.meta)  

## Aggregate cells by desired metadata columns
## results get saved in p2ens under aggregateMatrices
p2ens$aggregateCells(c('cluster','sample'))

ccm.aggr2 <- p2ens$aggregateMatrices[['cluster:sample']]                                                                                                             

## Metadata for aggregated samples (automate)
ccm.aggr2.meta <- data.frame(
    celltype = strpart(rownames(ccm.aggr2),':',1),
    sample = strpart(rownames(ccm.aggr2),':',2)
)
ccm.aggr2.meta$sample.type <- nbHelpers::strpart(ccm.aggr2.meta$sample,'-',2)
ccm.aggr2.meta$sample.name <- rownames(ccm.aggr2)                                                                                                                                                                                                                              
## save the metadata in the matrix                                                
p2ens$setAggregateMatrixMeta('cluster:sample',ccm.aggr2.meta)                                                                                                                                                 
## character vector of extracellular genes(make into gene metadata)
hs.membrane.hgnc <- readRDS('hs.membrane.hgnc.rds')                                                                                                                                                                   
## setup the contrasts to run, values must be in sample.type                                                                                      
comparisons <- list(       
    A.vs.B = c('A','B'), 
    B.vs.C = c('B','C')
)                                                           

## Run differential expression for each cell type in jcl3         
de.global <- lapply(comparisons, function(c1) {                             
    try({
        res <- getCorrectedDE.allTypes(ens.p2=p2ens, ## pagoda ensemble object to use
                                       cellfactor=jcl3, ## factor specifing cell types
                                       sample.type.comparison=c1, ## compare what with what
                                       membrane.gene.names=hs.membrane.hgnc, ## which genes are membrane gens
                                       n.cores=32, ## number of cores
                                       ## perform global('global') or cell type specific correction 'exclcurrent'?
                                       correction.method='global',  
                                       ## method for differential expression correction
                                       ## current deseq2 or t.test
                                       de.method = 'deseq2',
                                       ## fc.method: method for calculation of correction vector
                                       ## deseq2 -- use deseq2 moderated fold changes
                                       ## simple -- just get fold changes from cpms
                                       ## dummy -- don't perform any correction
                                       fc.method = 'deseq2',
                                       ## factor by which to weight the corrections globally
                                       ## 0 is same as setting dummy to fc.method
                                       correction.global.weight,
                                       ## cell types not to consider when generatin FC correction
                                       cell.types.fc.exclude=c('Tcells'), 
                                       verbose=T) ## verbosity
        res
    })
})   

## save results in the p2 ensembl object
p2ens$de.results[['global']] <- de.global

## export results as json (requires viewer front end)
p2ens$saveDEasJSON('global',saveprefix='web.de.global/') 

## export an index for the results
p2ens$saveDEindex('global',filename='de.global.index.html',defile.prefix='web.de.global/')    


