
generateSerializedAppsWithJC <- function(ens.p2, jc, prefix="") {
    lapply(names(ens.p2$p2objs), function(n) {
        p2 <- ens.p2$p2objs[[n]]
        jc2 <- jc[rownames(p2$counts)]
        extra.meta = list(
            jointClustering=p2.metadata.from.factor(jc2, displayname='joint clustering')
        );
        p2web <- basicP2web(p2,app.title=n,extra.meta,n.cores=4)
        p2web$serializeToStaticFast(binary.filename=paste0(prefix,n,'.bin'));
    })
    invisible(NULL)
}
