#' Join dataframe or SingleCellExperiment object with ContigCellDB object
#'
#' @param template Dataframe or SingleCellExperiment object to be joined with ccdb.
#' @param ccdb A ContigCellDB object.
#' @param join_fun Function used for the join operation.
#' @param by A character vector of variables to join by.
#'
#' @return [ContigCellDB()]
#' @export
#'
ccdb_join = function(template, ccdb, join_fun = dplyr::left_join, by = ccdb$cell_pk){
    if(class(ccdb) == "CellaRepertorium") stop('ccdb must have class CellaRepertorium')
    
    #join
    if(class(template) == "SingleCellExperiment") {
        #check if all keys ccdb keys are in template
        if(!all(by %in% colnames(SingleCellExperiment::colData(template)))) stop('Not all ccdb keys present in template')
        ccdb$cell_tbl = join_fun(as.data.frame(SingleCellExperiment::colData(template)), ccdb$cell_tbl, by = by)
    }
    else if(class(template) == "data.frame") {
        #check if all keys ccdb keys are in template
        if(!all(by %in% colnames(template))) stop('Not all ccdb keys present in template')
        ccdb$cell_tbl = join_fun(template, ccdb$cell_tbl, by = by)
    }
    else stop('Template must have class SingleCellExperiment or data.frame')
    
    return(ccdb)
}
