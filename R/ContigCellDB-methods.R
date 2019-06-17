# To Document
# primary_keys, primary_keys<- -- use `$` to access/replace. Table not checked for validity, but I think should be?
# canonicalize contigs on cells
# subset (subset a table, then equalize...)


# Todo
# canonicalize contigs on clusters
# combine methods

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
#' @importFrom methods new slot slot<- validObject
#' @rdname ContigCellDB-fun
#'
#' @examples
#' data(contigs_qc)
#' ContigCellDB(contigs_qc, contig_pk = c('barcode', 'pop', 'sample', 'contig_id'), cell_pk = c('barcode', 'pop', 'sample'))
ContigCellDB = function(contig_tbl, contig_pk, cell_tbl, cell_pk, cluster_tbl, cluster_pk = character()){
    valid_KeyedTbl(contig_tbl, contig_pk)
    equalized = FALSE
    if(missing(cell_tbl)){
        if(missing(cell_pk) || !is.character(cell_pk)) stop("If `cell_tbl` missing then `cell_pk` must name columns in `contig_tbl` that identify cells")
        cell_tbl = as_tibble(unique(contig_tbl[cell_pk]))
        equalized = TRUE
    } else{
        valid_KeyedTbl(cell_tbl, cell_pk)
    }
    if(missing(cluster_tbl)){
        cluster_tbl = tibble()
    } else{
        valid_KeyedTbl(cluster_tbl, cluster_pk)
    }
    obj = new('ContigCellDB', contig_tbl = contig_tbl, contig_pk = contig_pk, cell_tbl = cell_tbl, cell_pk = cell_pk, cluster_tbl = cluster_tbl, cluster_pk = cluster_pk, equalized = equalized)
    if(!equalized) equalize_ccdb(obj) else obj
}

#' @describeIn ContigCellDB-fun provide defaults that correspond to identifiers in 10X VDJ data
#' @export
ContigCellDB_10XVDJ = function(contig_tbl, contig_pk = c('barcode', 'contig_id'), cell_pk = 'barcode'){
    ContigCellDB(contig_tbl = contig_tbl, contig_pk = contig_pk, cell_pk = cell_pk)
}

setMethod("$", signature = c(x = 'ContigCellDB'), function(x, name){
    if(name %in% c('contig_tbl', 'cell_tbl', 'contig_pk', 'cell_pk', 'cluster_tbl', 'cluster_pk')){
        slot(x, name)
    } else{
        stop("Cannot access member ", name)
    }
})

setReplaceMethod("$", signature = c(x = 'ContigCellDB'), function(x, name, value){
    if(name %in% c('contig_tbl', 'cell_tbl', 'contig_pk', 'cell_pk', 'cluster_tbl', 'cluster_pk')){
        slot(x, name) <- value
        x@equalized = FALSE
    } else{
        stop("Cannot access member ", name)
    }
    if(name == 'contig_tbl'){
        valid_KeyedTbl(x$contig_tbl, x$contig_pk)
        x = equalize_ccdb(x)
    }
    if(name == 'cell_tbl'){
        valid_KeyedTbl(x$cell_tbl, x$cell_pk)
        x = equalize_ccdb(x)
    }
    if(name == 'cluster_tbl'){
        valid_KeyedTbl(x$cluster_tbl, x$cluster_pk)
        x = equalize_ccdb(x)
    }
    invisible(x)
})

setMethod('show', signature = c(object = 'ContigCellDB'), function(object){
    cat(class(object), "of", nrow(object$contig_tbl), "contigs")
    if((ncells <- nrow(object$cell_tbl)) > 0) cat(";", ncells, "cells;")
    cat(" and", nrow(object$cluster_tbl), "clusters")
    cat(".\n")
    cat('Contigs keyed by ', paste(object@contig_pk, collapse = ', '), '; cells keyed by ', sep = '')
    cat(paste(object@cell_pk, collapse = ', '), '.\n', sep = '')
})

#' @importFrom dplyr semi_join left_join right_join
equalize_ccdb = function(x){
    # Must use @ to avoid infinite loop!
    x@cell_tbl = semi_join(x$cell_tbl, x$contig_tbl, by = x$cell_pk)
    x@contig_tbl = semi_join(x$contig_tbl, x$cell_tbl, by = x$cell_pk)
    if(nrow(x$cluster_tbl) > 0) x@cluster_tbl = semi_join(x$cluster_tbl, x$contig_tbl, by = x$cluster_pk)
    x@equalized = TRUE
    x
}

replace_cluster_tbl = function(ccdb, cluster_tbl, contig_tbl, cluster_pk){
    if(nrow(ccdb$cluster_tbl)>0 && !missing(cluster_pk)){
        warning("Replacing `cluster_tbl` with key ccdb$cluster_pk")
    }
    if(!missing(cluster_pk)) ccdb$cluster_pk = cluster_pk
    ccdb@cluster_tbl = cluster_tbl
    if(!missing(contig_tbl)){
        ccdb@contig_tbl = contig_tbl
        valid_KeyedTbl(ccdb$contig_tbl, ccdb$contig_pk)
    }
    valid_KeyedTbl(ccdb$cluster_tbl, ccdb$cluster_pk)
    ccdb = equalize_ccdb(ccdb)
    validObject(ccdb)
    ccdb
}

