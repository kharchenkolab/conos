#' @export cpcaJCp4
cpcaJCp4 <- function(r.n, k=30, k.self=0, k.self.weight=1,community.detection.method = multilevel.community, reduction.method='CPCA', var.scale =TRUE, min.group.size = 10,ncomps=100, n.odgenes=1000, n.cores=30, return.details=F,xl=NULL,neighborhood.average=FALSE,neighborhood.average.k=5,groups=NULL, common.centering=TRUE, ...) {
    k1 <- k2 <- k;

    require(igraph)
    require(irlba)
    require(cpca)
    require(Matrix)

    ## neighborhood average
    if(neighborhood.average) {
        cat("neighborhood averaging ")
        r.n <- lapply(r.n,function(r) {
            xk <- clusterMatch:::crossNN(r$reductions$PCA,r$reductions$PCA,neighborhood.average.k,2,2.0,FALSE,n.cores)
            xk@x <- pmax(1-xk@x,0);
            diag(xk) <- 1;
            xk <- t(t(xk)/colSums(xk))
            colnames(xk) <- rownames(xk) <- rownames(r$reductions$PCA)
            r$misc$edgeMat$quickCPCA <- xk;
            cat(".")
            r
        })
        cat(" done\n")
    }

    ## Get combinations of sampes
    cis <- combn(names(r.n),2)

    ## Calculate pairwise alignments or use xl object
    if(is.null(xl)) {
        cat('pairwise',reduction.method,'(',ncol(cis),'pairs) ')
        xl <- pagoda2:::papply(1:ncol(cis), function(i) {
            if(reduction.method=='CPCA') {
                xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average)
            } else if(reduction.method=='JNMF') {
                xcp <- quickJNMF(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average,maxiter=3e3)
            } else if(reduction.method=='PCCPCA') {
                xcp <- quickPCCPCA(r.n[cis[,i]],ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale)
            } else {
                stop(paste("unknown reduction method",reduction.method))
            }
            cat('.')
            xcp
        },n.cores=n.cores,mc.preschedule=T);
        names(xl) <- apply(cis,2,paste,collapse='.vs.');
        cat("done\n")
    } else {
        ## match for all the pairs
        mi <- rep(NA,ncol(cis));
        nm <- match(apply(cis,2,paste,collapse='.vs.'),names(xl));
        mi[which(!is.na(nm))] <- na.omit(nm);
        nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(xl));
        mi[which(!is.na(nm))] <- na.omit(nm);
        cat('matched',sum(!is.na(mi)),'out of',length(mi),'CPCA pairs ... ')
        if(any(is.na(mi))) {
            cat('running',sum(is.na(mi)),'additional CPCAs ')
            xl2 <- pagoda2:::papply(which(is.na(mi)), function(i) {
                if(reduction.method=='CPCA') {
                    xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average)
                } else if(reduction.method=='JNMF') {
                    xcp <- quickJNMF(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average,maxiter=3e3)
                }
                cat('.')
                xcp
            },n.cores=n.cores);
            names(xl2) <- apply(cis[,which(is.na(mi)),drop=F],2,paste,collapse='.vs.');
            xl <- c(xl,xl2);
        }
        
        ## re-do the match and order
        mi <- rep(NA,ncol(cis));
        nm <- match(apply(cis,2,paste,collapse='.vs.'),names(xl));
        mi[which(!is.na(nm))] <- na.omit(nm);
        nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(xl));
        mi[which(!is.na(nm))] <- na.omit(nm);
        if(any(is.na(mi))) { stop("unable to get complete set of CPCA results") }
        xl <- xl[mi]
        cat(" done\n");
    }
    
    ## run mNN separatly as it can't deal with multithreading
    cat('mNN ')
    mnnres <- lapply(1:ncol(cis), function(i) {
        r.ns <- r.n[cis[,i]]
        if(!is.null(xl[[i]]$rot1)) {
            n12 <- clusterMatch:::crossNN(xl[[i]]$rot1,xl[[i]]$rot2,k,2,2.0,FALSE,n.cores)
            n21 <- clusterMatch:::crossNN(xl[[i]]$rot2,xl[[i]]$rot1,k,2,2.0,FALSE,n.cores)
            mnn <- drop0(n21*t(n12))
            mnn <- as(n21*t(n12),'dgTMatrix')
            cat(".")
            return(data.frame('mA.lab'=rownames(xl[[i]]$rot1)[mnn@i+1],'mB.lab'=rownames(xl[[i]]$rot2)[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
        } else if(!is.null(xl[[i]]$rotations)) {
            ## PCCPCA
            n12 <- clusterMatch:::crossNN(xl[[i]]$rotations[[1]],xl[[i]]$rotations[[2]],k,2,2.0,FALSE,n.cores)
            n21 <- clusterMatch:::crossNN(xl[[i]]$rotations[[2]],xl[[i]]$rotations[[1]],k,2,2.0,FALSE,n.cores)
            mnn <- drop0(n21*t(n12))
            mnn <- as(n21*t(n12),'dgTMatrix')
            cat(".")
            return(data.frame('mA.lab'=rownames(xl[[i]]$rotations[[1]])[mnn@i+1],'mB.lab'=rownames(xl[[i]]$rotations[[2]])[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
        } else {
            if(!is.null(xl[[i]]$CPC)) {
                ## CPCA
                rot <- xl[[i]]$CPC;
                odgenes <- rownames(xl[[i]]$CPC)
            } else if(!is.null(xl[[i]]$o$Q)) {
                ## GSVD
                rot <- xl[[i]]$o$Q;
                odgenes <- rownames(rot) <- colnames(xl[[i]]$o$A);
            } else {
                stop("unknown reduction provided")
            }
            if (var.scale) {
                cgsf <- do.call(cbind,lapply(r.ns,function(x) x$misc$varinfo[odgenes,]$gsf))
                cgsf <- exp(rowMeans(log(cgsf)))
            }
            ## create matrices, adjust variance
            cproj <- lapply(r.ns,function(r) {
                x <- r$counts[,odgenes];
                if(var.scale) {
                    x@x <- x@x*rep(cgsf,diff(x@p))
                }
                if(neighborhood.average) {
                    xk <- r$misc$edgeMat$quickCPCA;
                    x <- t(xk) %*% x
                }
                x
            })
            if(common.centering) {
                ncells <- unlist(lapply(cproj,nrow));
                centering <- colSums(do.call(rbind,lapply(cproj,colMeans))*ncells)/sum(ncells)
            } else {
                centering <- NULL;
            }
            cpproj <- lapply(cproj,function(x) {
                if(is.null(centering)) { centering <- colMeans(x) }
                x <- t(as.matrix(t(x))-centering)
                x %*% rot;
            })
            n1 <- cis[1,i]; n2 <- cis[2,i]
            cat(".")
            
            n12 <- clusterMatch:::crossNN(cpproj[[n1]],cpproj[[n2]],k,2,2.0,FALSE,n.cores)
            n21 <- clusterMatch:::crossNN(cpproj[[n2]],cpproj[[n1]],k,2,2.0,FALSE,n.cores)

            mnn <- drop0(n21*t(n12))
            mnn <- as(n21*t(n12),'dgTMatrix')
            return(data.frame('mA.lab'=rownames(cpproj[[n1]])[mnn@i+1],'mB.lab'=rownames(cpproj[[n2]])[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
            
        }
        mnnres
    })
    cat("done\n")
    
    ## Merge the results into a edge table
    el <- do.call(rbind,mnnres)

    ## Add self edges 
    if(k.self>0) {
        cat('kNN pairs ')
        x <- data.frame(do.call(rbind,lapply(r.n,function(x) {
            xk <- clusterMatch:::crossNN(x$reductions$PCA,x$reductions$PCA,k.self,2,2.0,FALSE,n.cores)
            diag(xk) <- 0;
            xk <- as(xk,'dgTMatrix')
            cat(".")
            return(data.frame('mA.lab'=rownames(x$reductions$PCA)[xk@i+1],'mB.lab'=rownames(x$reductions$PCA)[xk@j+1],'w'=pmax(1-xk@x,0),stringsAsFactors=F))
        })),stringsAsFactors = F)
        x$w <- k.self.weight
        cat(' done\n')
        el <- rbind(el,x)
    }

    ## Make graph and add weights
    g  <- graph_from_edgelist(as.matrix(el[,c(1,2)]), directed =FALSE)
    E(g)$weight <- el[,3]
    
    ## Do community detection on this graph
    cat('detecting clusters ...');
    cls <- community.detection.method(g, ...)
    cat('done\n')

    ## Extract groups from this graph
    cls.mem <- membership(cls)
    cls.groups <- as.character(cls.mem)
    names(cls.groups) <- names(cls.mem)
    
    ## Filter groups
    lvls.keep <- names(which(table(cls.groups)  > min.group.size))
    cls.groups[! as.character(cls.groups) %in% as.character(lvls.keep)] <- NA
    cls.groups <- as.factor(cls.groups)

    ## return
    if(return.details) {
        return(list(groups=cls.groups,xl=xl,cls=cls,mnnres=mnnres,g=g))
    } else {
        cls.groups
    }
}

