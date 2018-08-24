1. Download directory ShinyApp from git into a working directory. 
2. Prepare a conos object with walktrap community detection results.
3. Prepare leaf.labels factor according to greedy.modularity.cut function requirements. 
4. Open app.R in the saved ShinyApp directory in a running RStudio session.
5. Set parametres for greedy.modularity.cut function (it is not recommended to set "N" more than 50). 
6. Run the app.
7. Current cut membership may be saved into the conos object (con$clusters$walktrap) through pressing "Save membership" button.

Test data:
con <- readRDS("/home/yarloz/RDAs_DFs/0813_conosobj_multilevel_hca_bm_unfiltered.rds")
leaf.labels <- con$clusters$walktrap$result$names
leaf.labels <- sapply(leaf.labels, function(x) strsplit(strsplit(x,"_")[[1]][1],"Manton")[[1]][2])

