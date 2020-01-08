# Contig -- always equalize to cell
# Cluster -- equalize to contig
# Cell -- equalize after calling

valid_KeyedTbl = function(tbl, keys){
    tbl_nm = deparse(substitute(tbl))
    if( length(missing_fields <- setdiff(keys, names(tbl))) > 0){
        stop(sprintf("%s fields were named as primary keys but were missing from %s",
                     paste(missing_fields, collapse = ', '), tbl_nm))
    }

    if(any(dups <- duplicated(tbl[keys]))){
        stop(sprintf("In %s, rows %s... have identical `keys`,
                     which must uniquely identify a row.", tbl_nm,
                     paste(head(which(dups)), collapse = ',')))
    }
    if(tibble::has_rownames(tbl)){
        warning(sprintf('rownames in %s will be ignored.', tbl_nm))
    }
    invisible(TRUE)
}

#' Construct a ContigCellDB
#'
#' @param contig_tbl a data frame of contigs, and additional fields describing their properties
#' @param contig_pk character vector naming fields in `contig_tbl` that uniquely identify a row/contig
#' @param cell_tbl a data frame of cell barcodes, and (optional) additional fields describing their properties
#' @param cell_pk character vector naming fields in `cell_tbl` that uniquely identify a cell barcode
#' @param cluster_tbl A data frame that provide cluster assignments for each contig
#' @param cluster_pk If `cluster_tbl` was provided, a character vector naming fields in `cluster_tbl` that uniquely identify a cluster
#' @param equalize `logical`. Should the contig, cells and clusters be equalized by taking the intersection of their common keys?
#' @return \code{ContigCellDB}
#'
#' @section Accessors/mutators:
#' See \code{\link[=ContigCellDB-mutate]{$,ContigCellDB-method}} for more on how to access and mutate slots.
#' See [mutate_cdb()] and [filter_cdb()] for endomorphic filtering/mutation methods
#' See [split_cdb()] to split into a list, and [rbind.ContigCellDB()] for the inverse operation.
#' @export
#' @importFrom S4Vectors List SimpleList
#' @importFrom tibble as_tibble
#' @importFrom methods new slot slot<- validObject
#' @rdname ContigCellDB-fun
#' @seealso \code{\link[=ContigCellDB-mutate]{$,ContigCellDB-method}}
#' @examples
#' data(contigs_qc)
#' contigs_qc
#'
#' cdb = ContigCellDB(contigs_qc, contig_pk = c('barcode', 'pop', 'sample', 'contig_id'),
#'  cell_pk = c('barcode', 'pop', 'sample'))
#'  cdb
#'
#'  # everything that was in contigs_qc
#'  cdb$contig_tbl
#'
#'  # Only the cell_pk are included by default (until clustering/canonicalization)
#'  cdb$cell_tbl
#'
#'  # Empty, since no cluster_pk was specified
#'  cdb$cluster_tbl
#'
#'  # Keys
#'  cdb$contig_pk
#'  cdb$cell_pk
#'  cdb$cluster_pk
ContigCellDB = function(contig_tbl, contig_pk, cell_tbl, cell_pk,
                        cluster_tbl, cluster_pk = character(),  equalize = TRUE){
    contig_tbl = as_tibble(contig_tbl)
    valid_KeyedTbl(contig_tbl, contig_pk)
    if(missing(cell_tbl)){
        if(missing(cell_pk) || !is.character(cell_pk))
            stop("If `cell_tbl` missing then `cell_pk` must name columns in `contig_tbl`
                 that identify cells")
        cell_tbl = as_tibble(unique(contig_tbl[cell_pk]))
        equalize = TRUE
    } else{
        cell_tbl = as_tibble(cell_tbl)
        valid_KeyedTbl(cell_tbl, cell_pk)
    }
    if(missing(cluster_tbl)){
        if(length(cluster_pk)>0){
            cluster_tbl = as_tibble(unique(contig_tbl[cluster_pk]))
        } else{
            cluster_tbl = tibble()
        }
    }
    valid_KeyedTbl(cluster_tbl, cluster_pk)
    obj = new('ContigCellDB', contig_tbl = contig_tbl, contig_pk = contig_pk,
              cell_tbl = cell_tbl, cell_pk = cell_pk, cluster_tbl = cluster_tbl,
              cluster_pk = cluster_pk, equalized = equalize)
    if(equalize) equalize_ccdb(obj) else obj
}

#' @describeIn ContigCellDB-fun provide defaults that correspond to identifiers in 10X VDJ data
#' @param ... passed to [ContigCellDB()]
#' @export
ContigCellDB_10XVDJ = function(contig_tbl, contig_pk = c('barcode', 'contig_id'), cell_pk = 'barcode', ...){
    ContigCellDB(contig_tbl = contig_tbl, contig_pk = contig_pk, cell_pk = cell_pk, ...)
}

access_cdb = function(x, name){
    if(name %in% c('contig_tbl', 'cell_tbl', 'contig_pk', 'cell_pk', 'cluster_tbl', 'cluster_pk')){
        slot(x, name)
    } else if(name %in% names(x@cell_tbl)) {
        x@cell_tbl[[name]]
    } else{
        stop("Cannot access member ", name)
    }
}

replace_cdb = function(x, name, value){
    if(name %in% c('contig_tbl', 'cell_tbl', 'contig_pk', 'cell_pk', 'cluster_tbl', 'cluster_pk')){
        slot(x, name) <- value
        x@equalized = FALSE
    } else if(name %in% names(x@cell_tbl)) {
        x@cell_tbl[[name]] <- value
    } else{
        stop("Cannot access member ", name)
    }
    if(name %in% c('contig_tbl','contig_pk')){
        valid_KeyedTbl(x$contig_tbl, x$contig_pk)
        x = equalize_ccdb(x, contig = TRUE, cell = FALSE, cluster = TRUE)
    }
    if(name %in% c('cell_tbl','cell_pk')){
        valid_KeyedTbl(x$cell_tbl, x$cell_pk)
        x = equalize_ccdb(x, contig = TRUE, cell = FALSE, cluster = FALSE)
    }
    if(name %in% c('cluster_tbl','cluster_pk')){
        valid_KeyedTbl(x$cluster_tbl, x$cluster_pk)
        x = equalize_ccdb(x, contig = TRUE, cell = FALSE, cluster = TRUE)
    }
    invisible(x)
}

#' Access public members of ContigCellDB object
#'
#' @param x A ContigCellDB object
#' @param name a slot of a ContigCellDB object (one of  `c('contig_tbl', 'cell_tbl', 'contig_pk', 'cell_pk', 'cluster_tbl', 'cluster_pk')`)
#'
#' @return Update or return a slot of [ContigCellDB()]
#' @export
#' @aliases ContigCellDB-mutate
#' @examples
#' ccdb_ex$contig_tbl
#' ccdb_ex$cell_tbl
#' ccdb_ex$cluster_tbl
setMethod("$", signature = c(x = 'ContigCellDB'), access_cdb)

#' @param value The value assigned to a slot of ContigCellDB object
#' @rdname cash-ContigCellDB-method
#' @export
#'
#' @examples
#' ccdb_ex$contig_pk = c("sample","barcode","contig_id") # 'pop' is technically redundant with 'sample'
#' # Take a subset of ccdb_ex
#' ccdb_ex
#' ccdb_ex$contig_tbl = dplyr::filter(ccdb_ex$contig_tbl, pop == 'b6')
#' ccdb_ex
setReplaceMethod("$", signature = c(x = 'ContigCellDB'),
                 function(x, name, value) replace_cdb(x, name, value))

setMethod('show', signature = c(object = 'ContigCellDB'), function(object){
    cat(class(object), "of", nrow(object$contig_tbl), "contigs")
    if((ncells <- nrow(object$cell_tbl)) > 0) cat(";", ncells, "cells;")
    cat(" and", nrow(object$cluster_tbl), "clusters")
    cat(".\n")
    cat('Contigs keyed by ', paste(object@contig_pk, collapse = ', '),
        '; cells keyed by ', sep = '')
    cat(paste(object@cell_pk, collapse = ', '), '.\n', sep = '')
})


#' `data.frame`-like mutation/accessor generics for `ContigCellDB` objects
#'
#' A `ContigCellDB` pretend to be a `cell_tbl` data.frame in several regards.
#' This is to enable nesting `ContigCellDB` objects in the `colData` of a `SingleCellExperiment`
#' and so that various plotting functionality in `scater` can do something sensible.
#'
#' If `x` a `ContigCellDB`, then `dim(x)` and `dimnames(x)` return `dim(x$cell_tbl)` and `dimnames(x$cell_tbl)`, respectively, and `x[[col]]`  returns `x$cell_tbl[[col]]`.
#' Likewise indexing with `x[i,]` returns cells indexed by `i`.
#' Finally `as.data.frame(x)` returns `x$cell_tbl`.
setMethod('[[', signature = c(x = 'ContigCellDB', i = 'character', j = 'missing'), function(x, i, ...){
    x@cell_tbl[[i]]
})


#' @rdname sub-sub-ContigCellDB-character-missing-method
#' @param j ignored
#' @param drop ignored
#' @param x `ContigCellDB`
#' @param i integer or character index
#' @param ... ignored
#' @return See details.
#' @aliases [,ContigCellDB,ANY,missing-method
#' @examples
#'  data(ccdb_ex)
#'  ccdb_ex[1:10,]
#'  head(ccdb_ex[['barcode']])
#'  dim(ccdb_ex)
#'  dimnames(ccdb_ex)
setMethod('[', signature = c(x = 'ContigCellDB', i = 'ANY', j = 'missing'),
          function(x, i, ...){
            i = S4Vectors::NSBS(i, x)
            y = x
            y$cell_tbl = x$cell_tbl[i@subscript,]
            y
        })


# Should this be c(ncol(x$cell_tbl), nrow(x$cell_tbl))? Seems unlikely..
#' @rdname sub-sub-ContigCellDB-character-missing-method
setMethod('dim', signature = c(x = 'ContigCellDB'), function(x) dim(x$cell_tbl))

#' @rdname sub-sub-ContigCellDB-character-missing-method
setMethod('dimnames', signature = c(x = 'ContigCellDB'), function(x){
   dimnames(x$cell_tbl)
})

#' @rdname sub-sub-ContigCellDB-character-missing-method
setMethod('nrow', signature = c(x = 'ContigCellDB'), function(x){
    nrow(x$cell_tbl)
})

#' @rdname sub-sub-ContigCellDB-character-missing-method
setMethod('ncol', signature = c(x = 'ContigCellDB'), function(x){
    ncol(x$cell_tbl)
})

setMethod('NROW', signature = c(x = 'ContigCellDB'), function(x){
    NROW(x$cell_tbl)
})

setMethod('NCOL', signature = c(x = 'ContigCellDB'), function(x){
    NCOL(x$cell_tbl)
})

setMethod('showAsCell', signature = c(object = 'ContigCellDB'), function(object){
    not_cellkey = setdiff(names(object$cell_tbl), object$cell_pk)
    showAsCell(object$cell_tbl[not_cellkey])
})

setAs('ContigCellDB', 'data.frame', function(from){
    from$cell_tbl
})

as.data.frame.ContigCellDB = function(object) as(object, 'data.frame')


#' Take the intersection of keys in tables in `x`
#'
#' The cells in `cell_tbl`, and clusters in `cluster_tbl` can potentially be a superset of the `contig_tbl`.
#'
#' *  `equalize_ccdb(x, cell = TRUE)` trims cells that aren't in `contig_tbl` or  `cluster_tbl`.
#' *  `equalize_ccdb(x, cluster = TRUE)` trims clusters that aren't in `contig_tbl`.
#' *  `equalize_ccdb(x, contig = TRUE)` trims contigs that aren't `cell_tbl` or `cluster_tbl`.
#' @param x [ContigCellDB()]
#'
#' @param cell `logical` equalize cells
#' @param contig `logical` equalize contigs
#' @param cluster `logical` equalize clusters
#' @return [ContigCellDB()]
#'
#' @importFrom dplyr semi_join left_join right_join
#' @export
equalize_ccdb = function(x, cell = TRUE, contig = TRUE, cluster = TRUE){
    # Must use @ to avoid infinite loop!
    if(contig){
        x@contig_tbl = semi_join(x$contig_tbl, x$cell_tbl, by = x$cell_pk)
        if(nrow(x$cluster_tbl) > 0)
        x@contig_tbl = semi_join(x$contig_tbl, x$cluster_tbl, by = x$cluster_pk)
    }
    if(cell){
        if(nrow(x$cluster_tbl) > 0 && !contig){
            contig_tbl = semi_join(x$contig_tbl, x$cluster_tbl, by = x$cluster_pk)
        } else {
            contig_tbl = x$contig_tbl
        }
        x@cell_tbl = semi_join(x$cell_tbl, x$contig_tbl, by = x$cell_pk)

    }
    if(nrow(x$cluster_tbl) > 0 && cluster){
        x@cluster_tbl = semi_join(x$cluster_tbl, x$contig_tbl, by = x$cluster_pk)
    }
    x@equalized = (cell & contig & cluster) | x@equalized
    x
}

replace_cluster_tbl = function(ccdb, cluster_tbl, contig_tbl, cluster_pk){
    if(nrow(ccdb$cluster_tbl)>0 && !missing(cluster_pk)){
        warning("Replacing `cluster_tbl` with ", paste(ccdb$cluster_pk, sep = ', '), '.')
    }
    if(!missing(cluster_pk)) ccdb@cluster_pk = cluster_pk
    ccdb@cluster_tbl = cluster_tbl
    if(!missing(contig_tbl)){
        if(nrow(contig_tbl) != nrow(ccdb@contig_tbl)) stop("Length mismatch; this is a bug!")
        ccdb@contig_tbl = contig_tbl
        valid_KeyedTbl(ccdb$contig_tbl, ccdb$contig_pk)
    }
    valid_KeyedTbl(ccdb$cluster_tbl, ccdb$cluster_pk)
    validObject(ccdb)
    ccdb
}

#' @describeIn mutate_cdb Filter rows of a table in a `ContigCellDB` object
#' @export
#' @examples
#' data(ccdb_ex)
#' subset_contig = filter_cdb(ccdb_ex,full_length, productive == 'True',
#' high_confidence, chain != 'Multi', nchar(cdr3) > 5)
#' subset_cell = filter_cdb(ccdb_ex, sample == 4, tbl = 'cell_tbl')
filter_cdb <- function(ccdb, ..., tbl='contig_tbl'){
    thetbl <- access_cdb(ccdb,tbl)
    thetbl <- dplyr::filter(.data=thetbl,!!!rlang::quos(...))
    ccdb <- replace_cdb(ccdb,tbl,thetbl)
    return(ccdb)
}

#' Create new or update existing columns of `ContigCellDB` tables
#'
#' @param ccdb [ContigCellDB()]
#' @param ... name and value pair of column that will be updated
#' @param tbl `character.` One of `contig_tbl`, `cell_tbl` or `cluster_tbl`, naming the table to be updated.
#'
#' @return ContigCellDB object with updated table
#' @seealso [dplyr::mutate()]
#' @seealso [dplyr::filter()]
#' @export
#'
#' @examples
#' data(ccdb_ex)
#' new_contig = mutate_cdb(ccdb_ex, new_col = 1)
#' new_cell = mutate_cdb(ccdb_ex, new_col = 1, tbl = 'contig_tbl')
mutate_cdb <- function(ccdb, ..., tbl='contig_tbl'){
    thetbl <- access_cdb(ccdb,tbl)
    thetbl <- thetbl %>% dplyr::mutate(!!!rlang::quos(...))
    ccdb <- replace_cdb(ccdb,tbl,thetbl)
    return(ccdb)
}


#' Split into a list of [ContigCellDB()] by named fields
#'
#' @param ccdb [ContigCellDB()]
#' @param fields `character` naming fields in `tbl`
#' @param tbl one of `contig_tbl`, `cell_tbl` or `cluster_tbl`
#' @inheritParams base::split
#' @inheritParams ContigCellDB
#' @return list of `ContigCellDB`
#' @export
#'
#' @examples
#' data(ccdb_ex)
#' splat = split_cdb(ccdb_ex, 'chain', 'contig_tbl')
#' stopifnot(all(splat$TRA$contig_tbl$chain == 'TRA'))
#' stopifnot(all(splat$TRB$contig_tbl$chain == 'TRB'))
split_cdb = function(ccdb, fields, tbl = 'contig_tbl', drop = FALSE, equalize = TRUE){
    thetbl = access_cdb(ccdb, tbl)
    if(!is.character(fields)){
        stop('`field` must be a character naming fields in `tbl`')
    }
    if(length(missing <- setdiff(fields, names(thetbl)))>0){
        stop("The following fields are missing from ", tbl, ': ',
             paste0(missing, collapse = ', '), '.')
    }
    split_tbl = split(thetbl, thetbl[fields], drop = drop)
    out = purrr::map(split_tbl, function(tt){
        tmp = replace_cdb(ccdb, tbl, tt)
        if(equalize) equalize_ccdb(tmp)
    })
    out
}

#' @export
rbind.ContigCellDB <- function(..., deparse.level=1)
{
    objects <- list(...)
    .bind_rows_ccdb(objects[[1L]], objects[-1L])
}

#' Combine `ContigCellDB` along rows (contigs, cells or clusters).
#'
#' The union of the rows in each of the objects is taken,
#'  thus removing any rows that has an exact duplicate.  This
#'  includes all fields, not just the primary key for that table.
#' The union of the various primary keys is taken.
#' @param ... [ContigCellDB()]
#' @param deparse.level ignored
#' @return [ContigCellDB()]
#' @aliases rbind.ContigCellDB
#' @export
#' @importFrom S4Vectors rbind
#' @examples
#' data(ccdb_ex)
#' splat = split_cdb(ccdb_ex, 'chain', 'contig_tbl')
#' unite = rbind(splat$TRA, splat$TRB)
#' stopifnot(all.equal(unite, ccdb_ex))
#'
setMethod("rbind", "ContigCellDB",
          function (..., deparse.level = 1)
          {
              rbind.ContigCellDB(..., deparse.level = deparse.level)
          }
)


.bind_rows_ccdb = function(o1, objects, .id = NULL){
    all_objs = c(o1, objects)
    if(!all(purrr::map_lgl(all_objs, inherits, 'ContigCellDB')))
        stop("Can't rbind heterogenous objects.")

    tbls = list()
    for(tt in c('contig_tbl', 'cell_tbl', 'cluster_tbl')){
        tbls[[tt]] = unique(purrr::map_dfr(all_objs, access_cdb, name = tt, .id = .id))
    }

    pks = list()
    for(tt in c('contig_pk', 'cell_pk', 'cluster_pk')){
        pks[[tt]] = purrr::reduce(purrr::map(all_objs, access_cdb, name = tt), union)
    }
    ContigCellDB(tbls$contig_tbl, pks$contig_pk, tbls$cell_tbl, pks$cell_pk,
                 tbls$cluster_tbl, pks$cluster_pk)
}
