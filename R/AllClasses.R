setClass(
    'FineCluster',
    contains = 'list',
    slots = c(
        cluster = 'ANY',
        distance = 'matrix',
        medoid = 'integer',
        homology = 'numeric'
    )
)

setClass(
    "ContigCellDB",
    slots = c(
        contig_tbl = 'data.frame',
        contig_pk = 'character',
        cell_tbl = 'data.frame',
        cell_pk = 'character',
        cluster_pk = 'character',
        cluster_tbl = 'data.frame',
        cluster_type = 'character',
        equalized = 'logical'
    )
)


valid_ContigCellDB = function(object){
    if( length(missing_fields <- setdiff(object$cell_pk, names(object$contig_tbl))) > 0){
        return(sprintf("%s fields were named as primary keys in `cell_tbl`
                       but were missing from `contig_tbl`",
                       paste(missing_fields, collapse = ', ')))
    }
    if( length(missing_fields <- setdiff(object$cluster_pk, names(object$contig_tbl))) > 0){
        return(sprintf("%s fields were named as primary keys in `cluster_pk`
                       but were missing from `contig_tbl`",
                       paste(missing_fields, collapse = ', ')))
    }
    return(TRUE)
}

setValidity('ContigCellDB', valid_ContigCellDB)

