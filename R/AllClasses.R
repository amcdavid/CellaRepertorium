setClass('FineCluster', contains = 'list', slots = c(cluster = 'ANY', distance = 'matrix', medoid = 'integer', homology = 'numeric'))

valid_KeyedTbl = function(tbl, keys){
    tbl_nm = deparse(substitute(tbl))
    if( length(missing_fields <- setdiff(keys, names(tbl))) > 0){
        stop(sprintf("%s fields were named as primary keys but were missing from %s", paste(missing_fields, collapse = ', '), tbl_nm))
    }

    if(any(dups <- duplicated(tbl[keys]))){
        stop(sprintf("In %s, rows %s... have identical `keys`, which must uniquely identify a row.", tbl_nm,
                       paste(head(which(dups)), collapse = ',')))
    }
    invisible(TRUE)
}

#
# setValidity('KeyedTblS4', valid_KeyedTbl)
#
# KeyedTbl = function(tbl = tibble(), primary_keys = character()){
#     #new('KeyedTbl', tbl, primary_keys = primary_keys)
#     object = new_tibble(tbl, primary_keys = primary_keys, subclass = 'KeyedTbl')
#     if(!isTRUE(validmsg <- valid_KeyedTbl(object))) stop(validmsg)
#     object
# }
#
# primary_keys = function(object, ...){
#     UseMethod("primary_keys")
# }
#
# `primary_keys<-` = function(object, ...){
#     UseMethod("primary_keys<-")
# }
#
# primary_keys.KeyedTbl = function(object){
#     attr(object, 'primary_keys')
# }
#
# `primary_keys<-.KeyedTbl` = function(object, value){
#     attr(object, 'primary_keys') = value
#     if(!isTRUE(validmsg <- valid_KeyedTbl(object))) stop(validmsg)
#     object
# }

# setOldClass("KeyedTbl", S4Class = 'KeyedTblS4')

setClass("ContigCellDB", slots = c(contig_tbl = 'data.frame', contig_pk = 'character', cell_tbl = 'data.frame', cell_pk = 'character', cluster_tbls = 'SimpleList', cluster_pk = 'SimpleList', equalized = 'logical'))

#' Construct a ContigCellDB
#'
#' @param contig_tbl a data frame of contigs, and additional fields describing their properties
#' @param contig_pk character vector naming fields in `contig_tbl` that uniquely identify a row/contig
#' @param cell_tbl a data frame of cell barcodes, and (optional) additional fields describing their properties
#' @param cell_pk character vector naming fields in `cell_tbl` that uniquely identify a cell barcode
#' @param cluster_tbls An optional list of data frames that provide cluster assignments for each contig
#' @param cluster_pk If `cluster_tbls` was provided, a list of character vector naming fields in `cluster_tbls` that uniquely identify a cluster
#'
#' @return \code{ContigCellDB}
#' @export
#' @importFrom S4Vectors List SimpleList
#' @importFrom tibble as_tibble
#' @rdname ContigCellDB-fun
#'
#' @examples
#' data(contigs_qc)
#' ContigCellDB(contigs_qc, contig_pk = c('barcode', 'pop', 'sample', 'contig_id'), cell_pk = c('barcode', 'pop', 'sample'))
ContigCellDB = function(contig_tbl, contig_pk, cell_tbl, cell_pk, cluster_tbls = List(), cluster_pk = List()){
    valid_KeyedTbl(contig_tbl, contig_pk)
    equalized = FALSE
    if(missing(cell_tbl)){
        if(missing(cell_pk) || !is.character(cell_pk)) stop("If `cell_tbl` missing then `cell_pk` must name columns in `contig_tbl` that identify cells")
        cell_tbl = as_tibble(unique(contig_tbl[cell_pk]))
        equalized = TRUE
    } else{
        valid_KeyedTbl(cell_tbl, cell_pk)
    }
    obj = new('ContigCellDB', contig_tbl = contig_tbl, contig_pk = contig_pk, cell_tbl = cell_tbl, cell_pk = cell_pk, cluster_tbls = cluster_tbls, cluster_pk = cluster_pk, equalized = equalized)
    if(!equalized) equalize(obj) else obj
}

valid_ContigCellDB = function(object){
    if( length(missing_fields <- setdiff(object$cell_pk, names(object$contig_tbl))) > 0){
        return(sprintf("%s fields were named as primary keys in `cell_tbl` but were missing from `contig_tbl`", paste(missing_fields, collapse = ', ')))
    } else{
        return(TRUE)
    }
}

setValidity('ContigCellDB', valid_ContigCellDB)

#' @describeIn ContigCellDB-fun provide defaults that correspond to identifiers in 10X VDJ data
#' @export
ContigCellDB_10XVDJ = function(contig_tbl, contig_pk = c('barcode', 'contig_id'), cell_pk = 'barcode'){
    ContigCellDB(contig_tbl = contig_tbl, contig_pk = contig_pk, cell_pk = cell_pk)
}

setMethod("$", signature = c(x = 'ContigCellDB'), function(x, name){
    if(name %in% c('contig_tbl', 'cell_tbl', 'contig_pk', 'cell_pk')){
        slot(x, name)
    } else{
        stop("Cannot access member", name)
    }
})

setReplaceMethod("$", signature = c(x = 'ContigCellDB'), function(x, name, value){
    if(name %in% c('contig_tbl', 'cell_tbl', 'contig_pk', 'cell_pk')){
        slot(x, name) <- value
        x@equalized = FALSE
    } else{
        stop("Cannot access member", name)
    }
    invisible(x)
})

setMethod('show', signature = c(object = 'ContigCellDB'), function(object){
    cat(class(object), "of", nrow(object$contig_tbl), "contigs")
    if((ncells <- nrow(object$cell_tbl)) > 0) cat(";", ncells, "cells;")
    cat(" and", length(cluster_tbls(object)), "cluster tables")
    cat(".")
})

equalize_ContigCellDB = function(x, ...){
    if(length(list(...)) > 0) stop('No additional arguments accepted')
    x$cell_tbl = semi_join(x$cell_tbl, x$contig_tbl, by = x$cell_pk)
    x$contig_tbl = semi_join(x$contig_tbl, x$cell_tbl, by = x$cell_pk)
    x@equalized = TRUE
    x
}

setGeneric('equalize', function(x, ...) standardGeneric('equalize'))

setMethod('equalize', signature = c(x = 'ContigCellDB'), equalize_ContigCellDB)


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

# Methods
# /primary_keys, primary_keys<-
# canonicalize contigs on cells, canonicalize contigs on clusters
# subset
# combine
