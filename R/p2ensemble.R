
#' An ensemble of p2 objects that can be manipulated together
#' @exportClass Pagoda2ensemble
Pagoda2ensemble <- setRefClass(
    "Pagoda2ensemble",
  
    fields=c('p2objs','rawMatrices','n.cores','jointClustering','samplesheet','de.results','de.results.json','full.matrix',
             'cell.metadata', 'aggregateMatrices', 'aggregateMatrixMeta'),
    
    methods=list(
        #' constructor
        #' @name Pagoda2ensemble object constructor
        initialize=function(x, ..., n.cores=parallel::detectCores(logical=F)) {
            p2objs <<- list();
            rawMatrices <<- list();
            jointClustering <<- list();
            de.results <<- list();
            de.results.json <<- list();
            aggregateMatrices <<- list();
            aggregateMatrixMeta <<- list();


            if(!missing(x) && class(x) == 'Pagoda2ensemble') {
                ## copy constructor
                callSuper(x, ..., n.cores=n.cores);
            } else {
                callSuper(..., n.cores = n.cores);
            }
        },
        
        #' @name Set the list of pagoda2 objects to work on
        #' @param p2objects list of pagoda2 objects
        setObjects=function(p2objects) {
            if(!is.list(p2objects))
                stop("'p2objects' is not a list");
            if(length(p2objects) < 2)
                stop("The length of p2objects needs to be at least 1");
            if(!all(unlist(lapply(p2objects, class)) == 'Pagoda2'))
                stop("Some of the p2objects are not of class Pagoda2");
            p2objs <<- p2objects
        },

        #' @name Extract count matrices for differential expression
        #' and proportion comparisons from pagoda2 apps
        #' @return a list of dgCMatrices
        extractCountMatrices=function() {
            ## Get common genes
            gl <- lapply(p2objs, function(r) colnames(r$counts))
            all.genes <- unique(unlist(gl))
            gc <- do.call(rbind, lapply(gl, function(x) all.genes %in% x))
            common.genes <- all.genes[apply(gc,2,sum) > 3]
            ## Extract count matrices for common gene list
            ccm.raw <- parallel::mclapply(p2objs,function(r) {
                om <- as(r$misc$rawCounts,'dgTMatrix')
                om <- om[,colnames(om) %in% common.genes]
                mi <- match(colnames(om), common.genes)
                x <- new("dgTMatrix",i = om@i,j = as.integer(mi[om@j+1]-1),x=om@x,Dim=c(nrow(om),length(common.genes)))
                rownames(x) <- rownames(om);
                colnames(x) <- common.genes
                as(x,'dgCMatrix')
            }, mc.cores =30)
            rawMatrices <<- ccm.raw
            invisible(ccm.raw)
        },

        #' @name Calculate joint clustering between these apps
        calcJointClustering=function() {
            stop('Not implemented')
        },

        #' @name Specify clustering that encompases cells from all the apps
        #' @param name name for this clustering
        #' @param jc named factor specifying the cluster membership
        setJointClustering=function(name=NULL, jc=NULL){
            if(is.null(name))
                stop("'name' must be specified");
            if(is.null(jc))
                stop("'jc' must be specified");
            if(!is.character(name))
                stop("'name must be a character value");
            ## TODO check integrity of jc
            jointClustering[[name]] <<- jc;
        },

        #' @name Get the genes that are in all the p2 objects in this ensemble
        #' @return names of common genes
        getCommonGenes=function() {
            app.genes <- lapply(p2objs, function(o) { colnames(o$counts) })
            Reduce(intersect, app.genes);
        },

        #' @name Helper function for getting panel dimentions
        #' @param n.items items that we want to show
        #' @param square logical: force square shape
        #' @return integer vector with 2 numbers
        getParMfrow = function(n.items, square=FALSE) {
            n <- ceiling(sqrt(n.items))
            if(square) {
                c(n,n);
            } else {
                c(n,ceiling(n.items/n))
            }
        },
        
        #' @name Plot all the apps with the designated groups
        #' @param groups a factor of groups
        #' @param filename if not NULL save to file
        #' @param panel.size panel size for saving to file
        #' @param mark.cluster.cex cex for marking clusters
        #' @param mar mar parameter for par
        #' @param mgp mgp parameter for par
        #' @param cex cex parameter for par
        #' @param type type of embedding to plot
        #' @param embeddingType embeddingType of embedding to plot
        #' @param mark.clusters show names of clusters?
        plotWithGroups = function(groups=NULL,filename=NULL,panel.size=600,mark.cluster.cex=0.8,
                                  mar=c(0.5,0.5,0.5,0.5), mgp = c(2,0.65,0), cex = 0.85, type='PCA',
                                  embeddingType='tSNE',mark.clusters=TRUE) {
            require(Cairo)
            panel.dims <- getParMfrow(length(p2objs))
            if(!is.null(filename))
                CairoPNG(file=filename,height=paneldims[1],width=paneldims[2])
            par(mfrow=c(panel.dims[1],panel.dims[2]), mar = mar, mgp = mgp, cex = cex);
            lapply(names(p2objs),function(dn) {
                d <- p2objs[[dn]];
                g1 <- as.factor(groups)
                colors <- NULL
                ## If no cells present fix
                if (!any(names(g1) %in% rownames(d$counts))) {
                    g1 <- NULL
                    cell.names <- rownames(d$counts)
                    colors <- rep('grey70',length(cell.names))
                    names(colors) <- cell.names
                }
                d$plotEmbedding(type=type,embeddingType=embeddingType,groups=g1,alpha=0.2,
                                min.group.size=0,mark.clusters = mark.clusters,
                                mark.cluster.cex=mark.cluster.cex,do.par=F,colors=colors);
                legend(x='topleft',bty='n',legend=dn)
            })
            if(!is.null(filename))
                dev.off()
            invisible(NULL);
        },
        
        #' @name Plot all apps with the clusters specified by the joint clustering
        #' @param ... parameters to pass to plotWithGroups
        plotWithJointclustering=function(jc.name=NULL, ...) {
            if(is.null(jc.name))
                stop("'jc.name' cannot be NULL")
            if(is.null(jointClustering[[jc.name]]))
                stop("The specified joint clustering does not exist");
            plotWithGroups(groups=jointClustering[[jc.name]]$groups, ...)
        },

        #' @name Helper function for calculation of common zlim
        #' @param vs values
        #' @param gradient.range.quantile parameter for trimming
        #' @return zlim values
        calcZlim = function(vs, gradient.range.quantile = 0.95) {
            zlim <- as.numeric(quantile(vs,p=c(1-gradient.range.quantile,gradient.range.quantile)))
            if(diff(zlim)==0) {
                zlim <- as.numeric(range(vs))
            }
            zlim
        },

        #' @name plot embeddings of all objects colored by depth
        #' @param filename if not NULL save to file
        #' @param panel.size panel size for saving to file
        #' @param mark.cluster.cex cex for marking clusters
        #' @return NULL
        plotWithDepth=function(filename=NULL,panel.size=600,mark.cluster.cex=0.8){
            require(Cairo)
            n <- ceiling(sqrt(length(p2objs)))
            if (!is.null(filename)) {
                CairoPNG(file=filename,height=n*panel.size,width=n*panel.size)
            }
            ## Get all depth values
            depthvalues <- unlist(lapply(p2objs, function(o) {o$depth}))
            zlim <- calcZlim(depthvalues)
            ## Do the plotting
            par(mfrow=c(n,n), mar = c(0.5,0.5,0.5,0.5), mgp = c(2,0.65,0), cex = 0.85);
            invisible(lapply(names(p2objs),function(dn) {
                d <- p2objs[[dn]];
                ## If no cells present fix
                d$plotEmbedding(type='PCA',embeddingType='tSNE',colors=d$depth,alpha=0.2,do.par=F, zlim=zlim);
                legend(x='topleft',bty='n',legend=dn)
            }))
            if(!is.null(filename)) {
                dev.off()
            }
        },
        
        #' @name set the sample sheet
        setSamplesheet=function(sample.sheet=NULL) {
            if(!is.data.frame(sample.sheet))
                stop('sample.sheet is not a dataframe')
            samplesheet <<- sample.sheet;
        },
        
        #' @name compare a set of cells accross the specified matrices with DESeq2
        #' @param rl1 list of test count matrices
        #' @param rl2 list of control count matrices
        #' @param n1 name for condition 1
        #' @param n2 name for condition 2
        #' @param min.cell minumum number of cells per cluster/app to include in the comparison
        #' @param membrane.gene.names names of membrane genes
        t.two.exp.comp.2  = function(rl1,rl2,cells,n1='test',n2='control',min.cells=10,
                                     membrane.gene.names = NULL) {
            require(Matrix)
            ## Subset apps to cells
            rl1 <- lapply(rl1,function(x) x[rownames(x) %in% cells,,drop=F])
            rl2 <- lapply(rl2,function(x) x[rownames(x) %in% cells,,drop=F])
            ## Check for empty samples
            ## Remove apps with < min.cells
            rl1 <- rl1[unlist(lapply(rl1,nrow))>min.cells]
            rl2 <- rl2[unlist(lapply(rl2,nrow))>min.cells]
            ## Don't run if there are no samples in one condition after filtering
            if (length(rl1)  == 0 | length(rl2) == 0) {
                return(NULL)
            }
            cts <- do.call(cbind,lapply(c(rl1,rl2),Matrix::colSums))
            ## Generate comparison factor with n2 as base level
            type.factor <- rep(c(n1,n2),c(length(rl1),length(rl2)));
            type.factor <- factor(type.factor, levels=c(n2,n1));
            coldata <- data.frame(type=type.factor);
            rownames(coldata) <- c(names(rl1),names(rl2))
            ## Make DESeq2 object
            dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,colData = coldata,design = ~ type)
            ## Do diff expression and get results
            dds <- DESeq2::DESeq(dds)
            ## Report the results without filtering
            res <- DESeq2::results(dds,cooksCutoff=FALSE,independentFiltering=FALSE)
            ## Get DESeq2 results
            allgenes <- rownames(res);
            res <- as.data.frame(res);
            res <- cbind(gene= allgenes, res);
            res <- res[!is.na(res$pvalue),]
            ## Mark everything with a q-val less the 0.05 as significant
            res$significant <- res$padj < 0.05
            ## Order by significance
            res <- res[order(res$padj),]
            ## TODO fix this into any kind of metadata
            if(!is.null(membrane.gene.names)) {
                res$membrane <- res$gene %in% membrane.gene.names;
            }
            ## Calculate Z and Z corrected score
            res$Z <- (-1) * qnorm(res$pvalue)*sign(res$log2FoldChange)
            res$Za <- (-1) * qnorm(res$padj)*sign(res$log2FoldChange)
            ## calculate mean expression values for the plots
            if(nrow(res)>0) {
                cts <- t(t(cts[as.character(res$gene),,drop=FALSE])/colSums(cts)*1e6)
                ilev <- tapply(1:ncol(cts),coldata$type,function(ii) {
                    cts[,ii,drop=F]
                })
            } else {
                ilev <- NULL;
            }
            ## return
            return(list(
                res=res,
                genes=allgenes,
                ilev=ilev,
                snames=c(n1,n2)
            ))
        },
        
        #' @name helper function that take a factor and breaks it down
        factorBreakdown = function(f) {
            if(!is.factor(f)) stop('not a factor!')
            lvls <- levels(f);
            names(lvls) <- lvls;
            lapply(lvls, function(l) {
                r <- names(f)[f == l]
                names(r) <- r
                r
            })
        },
        
        #' @name Calculated differential expression between apps of the specified type using deseq2
        #' @param app.types a named factor allocating each app to a type
        #' @param contrasts a named list of vectors that specified comparisons to run in terms of app.types
        #' @param jc.name name of joint clustering to used
        #' @param groups alternative clustering factor, used only when jc.name is NULL
        #' @param membreane.gene.names list of genes to be annotated as membrane genes
        #' @param de.name the name to use to save these differential results
        calcDiffExpression = function(app.types, contrasts, jc.name=NULL,groups=NULL,
                                      membrane.gene.names=NULL, de.name='default'){
            if(is.null(jc.name) & is.null(groups))
                stop('One of jc.names or groups needs to be defined');
            if(!is.null(jc.name)) {
                if(! jc.name %in% names(jointClustering))
                    stop('The joint clustering you specified was not found');
                groups <- jointClustering[[jc.name]]$groups;
            }
            if(is.null(rawMatrices) | length(rawMatrices) == 1)
                stop('"rawMatrices" is empty, run extractCountMatrices()');
            ccm <- rawMatrices
            clusters <- groups
            ## end param check
            clusters <- as.factor(clusters[!is.na(clusters)])
            clusters.split <- factorBreakdown(clusters)
            x <- lapply(contrasts, function(cc) {
                grp1 <- names(app.types)[app.types == cc[1]]
                grp2 <- names(app.types)[app.types == cc[2]]
                parallel::mclapply(clusters.split, function(cells.compare) {
                    res = t.two.exp.comp.2(ccm[grp1], ccm[grp2], cells.compare,
                                           n1=cc[1], n2=cc[2],
                                           membrane.gene.names=membrane.gene.names)
                    res
                }, mc.cores = n.cores)
            })
            de.results[[de.name]] <<- x;
            invisible(x)
        },

        #' @name Return the specified differential expression results
        #' @param de.name name of results
        #' @return differential expression results
        getDEresults = function(de.name='default') {
            if(!de.name %in% names(de.results))
                stop('Specified results not found');
            de.results[[de.name]]
        },
        
        #' @name Get a table that shows the number of cells in each cluster in each app
        #' @param jc.name name of joint clustering to use
        #' @param groups alternative groups to use, jc.name must be NULL
        getClusterCountsPerApp = function(jc.name=NULL,groups=NULL) {
            if(is.null(jc.name)) {
                if(is.null(groups)) {
                    stop("one of 'jc.name' and 'groups' needs to be specified");
                } else {
                    grp <- groups
                }
            } else {
                grp <- jointClustering[[jc.name]]$groups;
            }
            ## tmp data frame with information about cells and samples
            x <- do.call(rbind, mapply(function(o,n) {
                data.frame(
                    cellid=rownames(o$counts),
                    app=c(n),
                    stringsAsFactors=FALSE
                )
            }, p2objs, names(p2objs), SIMPLIFY = FALSE))
            x$cluster <- as.character(grp[x$cellid]);
            table(x[,c('app','cluster')])
        },
        
        #' @name get a table that shows the percent of cells in each cluster in each app
        #' @param ... parameters passed to getClusterCountsPerApp
        getClusterPercentPerApp = function(...) {
            z <- getClusterCountsPerApp(...);
            sweep(z, 1, Matrix::rowSums(z), '/');
        },

        #' @name get JSON file for interactive visualisation of expression
        #' @param jc.name name of joint clustering to use
        #' @param groups alternative factor to use for cell grouping, jc.name must be NULL to use this
        #' @param json.output output file name
        saveCLAVasJSON = function(jc.name = NULL, groups = NULL, json.output =  NULL) {
            if(is.null(jc.name) & is.null(groups))
                stop('jc.name or group must not be NULL');
            if(!is.null(jc.name)) {
                if(!jc.name %in% names(jointClustering))
                    stop('jc.name not found');
                groups <- jointClustering[[jc.name]]$groups;
            }
            if(is.null(json.output))
                stop('json.output is not specified')
            ## per-cluster average matrix
            get.cluster.pooled.averages <- function(r,groups,mult=1e3) {
                groups <- as.factor(groups[rownames(r$misc[['rawCounts']])])
                lvec <- pagoda2:::colSumByFac(r$misc[['rawCounts']],as.integer(groups))[-1,,drop=F];
                ## correct for the rare case when the last levels had no counts
                ## (should be corrected in the C implementation?
                if(nrow(lvec)<length(levels(groups)))
                    lvec <- rbind(lvec,matrix(0,ncol=ncol(lvec),nrow=length(levels(groups))-nrow(lvec)));
                lvec <- t(lvec/pmax(1,Matrix::rowSums(lvec)))*mult;
                rownames(lvec) <- colnames(r$misc[['rawCounts']]); colnames(lvec) <- levels(groups);
                lvec
            }
            get.cluster.means <- function(r,groups) {
                groups <- as.factor(groups[rownames(r$counts)])
                lvec <- pagoda2:::colSumByFac(r$counts,as.integer(groups))[-1,,drop=F];
                ## correct for the rare case when the last levels had no counts
                ## should be corrected in the C implementation?
                if(nrow(lvec)<length(levels(groups)))
                    lvec <- rbind(lvec,matrix(0,ncol=ncol(lvec),nrow=length(levels(groups))-nrow(lvec)));
                lvec <- t(lvec/as.integer(table(groups)))
                lvec[is.nan(lvec)] <- 0;
                rownames(lvec) <- colnames(r$counts); colnames(lvec) <- levels(groups);
                lvec
            }
            sn <- function(x) { names(x) <- x; x }
            ## Prepare the cell count table
            x <- reshape2::melt(.self$getClusterCountsPerApp(groups=groups))
            colnames(x) <- c('sample','type','cells')
            cell.count.table <- as.data.frame(x[,c('type','sample','cells')])
            rm(x)
            rownames(cell.count.table) <- paste(cell.count.table[,2],cell.count.table[,1],sep=":")
            cell.count.table <- cell.count.table[cell.count.table$cells>0,]   
            ## Calculate per cluster per sampels expression fro every gene
            gcl <- lapply(p2objs,get.cluster.means,groups=groups)
            common.genes <- Reduce(intersect, lapply(gcl,rownames))
            gcl <- lapply(gcl, function(o) {o[common.genes,]})

            gcl <- lapply(sn(names(gcl)), function(n) {
                x <- gcl[[n]];
                colnames(x) <- paste(n,colnames(x),sep=':');
                x
            })
            gcl <- do.call(cbind,gcl)
            gcl <- gcl[,rownames(cell.count.table)]    
            ## Get the sample names
            sample.names <- names(p2objs)
            ## TODO: allow custom order
            write(
                jsonlite::toJSON(
                              list(
                                  cellCounts=data.frame(cell.count.table),
                                  geneCounts=data.frame(t(gcl)),
                                  sampleNames=sample.names
                              ),
                              dataframe='columns',
                              digits=5
                          ),
                file=json.output
            )   
        },
        
        #' @name save differential expression results as JSON file for online viewing
        #' @description save the differential expression results as JSON, also generates
        #' the internal structure required to make an index for the files
        #' @param de.name name of differential expression results to save
        #' @param saveprefix prefix of json files
        saveDEasJSON = function(de.name = 'default', saveprefix = '') {
            if (! de.name %in% names(de.results))
                stop('"de.name" not found in de.results');
            comps <- de.results[[de.name]]; 
            ret <- lapply(namedNames(comps), function(ncc) {
                cc <- comps[[ncc]]
                lapply(namedNames(cc), function(nccc) {
                    xe <- cc[[nccc]]
                    if (!is.null(xe$res) && nrow(xe$res) > 0) {
                        ## Make the filename
                        file <- paste0(make.names(ncc), "__", make.names(nccc), ".json")
                        xe$res$rowid <- 1:nrow(xe$res)
                        xe$res <- data.frame(xe$res)
                        xe$ilev <- lapply(xe$ilev, function(x) list(snames=colnames(x),val=x))
                        y <- jsonlite::toJSON(xe)
                        write(y, file=paste0(saveprefix,file))
                        list(file = file, contrast = ncc, cluster = nccc)
                    } else {
                        NULL
                    }
                })
            })
            de.results.json[[de.name]] <<- ret
            invisible(ret)     
        },
        
        #' @name make an index html file for the designated comparisons
        #' @param de.name name of de expression field to use
        #' @param filename name of file to save
        #' @param deviewfile file of the html viewing interface to link to
        #' @param defile.prefix prefix for the links to the json files
        saveDEindex = function(de.name = 'default',filename='index_comparisons.html',
                               deviewfile = 'deview.2.html',defile.prefix='') {
            if(!any(de.name == names(de.results.json)))
                stop('de.name not found in de.results.json');
            g <- de.results.json[[de.name]]
            g0 <- lapply(g, function(g1) {
                g1 <- g1[!unlist(lapply(g1, is.null))]
                lapply(g1, function(x) {
                    paste0( '<a href="',deviewfile,'?d=',defile.prefix,x$file,'">', x$contrast, ': ',x$cluster,'</a>' )
                })
            })
            g0 <- reshape2::melt(g0)
            t0 <- reshape2::dcast(g0, L1 ~ L2, value.var='value')
            t0 <- t(t0)
            print(xtable::xtable(t0), sanitize.text.function=identity, type='html',file=filename)
        },

        ## '@name make a matrix by merging the individual object and store it in full.matrix
        ## '@return invisble full.matrix
        makeFullMatrix = function() {
            full.matrix <<- do.call(rbind, rawMatrices)
            invisible(full.matrix)
        },

        getFullMatrix = function() {
            invisible(full.matrix)
        },

        setCellMetadata = function(newMeta) {
            cell.metadata <<- newMeta
        },

        ## Sum raw cell counts according the cell.metadata fields aggr1
        aggregateCells = function(fields=NULL) {
            aggr1 <- apply(cell.metadata[,fields], 1, function(x) {paste(x,collapse=":")})
            ccm.aggr1 <- rowsum(as.matrix(full.matrix), aggr1)
            aggr.id <- paste(fields,collapse=':')
            aggregateMatrices[[aggr.id]] <<- ccm.aggr1
            invisible(aggregateMatrices[aggr.id])
        },

        # set metadata for aggregate matrix
        setAggregateMatrixMeta = function(aggregateName, meta) {
            aggregateMatrixMeta[[aggregateName]] <<- meta
        }
    )
)



