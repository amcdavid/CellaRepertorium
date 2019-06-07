get_cluster_tbls = function(x, index){
    if(missing(index)) x@cluster_tbls else x@cluster_tbls[[index]]
}
setGeneric('cluster_tbls', function(x, ...) standardGeneric('cluster_tbls'))
setMethod('cluster_tbls', signature = c(x = 'ContigCellDB'), get_cluster_tbls)

set_cluster_tbls = function(x, index, value){
    if(missing(index)) x@cluster_tbls = value else x@cluster_tbls[[index]] = value
    x
}

setGeneric('cluster_tbls<-', function(x, ..., value) standardGeneric('cluster_tbls<-'),  signature=c("x", "value"))

setReplaceMethod('cluster_tbls', signature = c(x = 'ContigCellDB'), set_cluster_tbls)
