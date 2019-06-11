setGeneric('cluster_tbls', function(x, ...) standardGeneric('cluster_tbls'))
setGeneric('cluster_tbls<-', function(x, ..., value) standardGeneric('cluster_tbls<-'),  signature=c("x", "value"))


get_cluster_tbls = function(x, index){
    if(missing(index)) x@cluster_tbls else x@cluster_tbls[[index]]
}

#' @export
setMethod('cluster_tbls', signature = c(x = 'ContigCellDB'), get_cluster_tbls)

set_cluster_tbls = function(x, index, value){
    if(missing(index)) x@cluster_tbls = value else x@cluster_tbls[[index]] = value
    x
}


#' @export
setReplaceMethod('cluster_tbls', signature = c(x = 'ContigCellDB'), set_cluster_tbls)
