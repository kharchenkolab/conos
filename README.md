# conos
Clustering on Network of Samples

# Usage example
```r
# construct conos object, where pl is a list of pagoda2 objects 
con <- Conos$new(pl)

# build graph
con$buildGraph(verbose=T)

# find communities
con$findCommunities()

# plot panel with joint clustering results
con$plotPanel()
  
# show a gene
con$plotPanel(gene='GZMK')

# embed jonit graph
con$embedGraph()

# plot joint clusters on the joint graph embedding
con$plotGraph()
# plot sample colors on the joint graph embedding
con$plotGraph(color.by='sample')
```
