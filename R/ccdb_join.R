#' Join dataframe or SingleCellExperiment object with ContigCellDB object
#'
#' @param template data.frame or SingleCellExperiment object to be joined with ccdb.
#' @param ccdb A ContigCellDB object.
#' @param join_fun Function used for the join operation.
#' @param by A character vector of variables to join by.
#'
#' @return [ContigCellDB()]
#' @export
#' @examples
#' data(ccdb_ex)
#' to_join = dplyr::bind_rows(ccdb_ex$cell_tbl[1:10,],
#' dplyr::tibble(barcode = c('extra1', 'extra2'), sample = LETTERS[1:2],
#' pop = LETTERS[1:2]))
#' ccdb_join(to_join, ccdb_ex)
ccdb_join = function(template, ccdb, join_fun = dplyr::left_join, by = ccdb$cell_pk){
    if(!inherits(ccdb,  "ContigCellDB")) stop('ccdb must have class CellaRepertorium')

    #join
    if(inherits(template,  "SingleCellExperiment")) {
        #check if all keys ccdb keys are in template
        if(!all(by %in% colnames(SingleCellExperiment::colData(template)))) stop('Not all ccdb keys present in template')
        ccdb$cell_tbl = join_fun(as.data.frame(SingleCellExperiment::colData(template)), ccdb$cell_tbl, by = by)
    }
    else if(inherits(template, "data.frame") || inherits(template, "DataFrame")) {
        #check if all keys ccdb keys are in template
        if(!all(by %in% colnames(template))) stop('Not all ccdb keys present in template')
        ccdb$cell_tbl = join_fun(as.data.frame(template), ccdb$cell_tbl, by = by)
    }
    else stop('Template must inherit from `SingleCellExperiment` or `data.frame`')

    return(ccdb)
}
