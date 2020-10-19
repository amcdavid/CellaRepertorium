check_plot_infra = function(){
    if(!requireNamespace('ggplot2')) stop("Install ggplot2 >= 3.0.0.")
    if(!requireNamespace('cowplot')) stop("Install cowplot.")
}

#' Make a plot showing properties of the clustering
#'
#' The number of elements per cluster and the average distance between the medoid and other elements are plotted.
#' @param cdb A `fine_clustering` `ContigCellDB` object
#' @param return_plotlist should  a list of `ggplot2` plots be returned.  If FALSE, a `cowplot` composite is retuned.
#'
#' @return  a `cowplot` composite or a list of plots.
#' @export
#'
#' @example inst/examples/small_cluster_example.R
#' @examples
#' cluster_plot(ccdb_ex_small)
cluster_plot = function(cdb, return_plotlist = FALSE){
    check_plot_infra()
    if(!has_fineclustering(cdb)) stop("Run `cdhit_cdb(cdb)` and/or `fine_clustering(cdb)` first.")
    dist_expanded = dplyr::filter(cdb$cluster_tbl, .data$n_cluster>1)
    n_cluster = cdb$cluster_tbl
    plts = list(
        ggplot2::ggplot(dist_expanded, ggplot2::aes(x = .data$avg_distance)) + ggplot2::geom_histogram() + ggplot2::xlab('Intra distance') + ggplot2::ggtitle(' Distance (non-singletons)', subtitle = cdb$cluster_pk),
        ggplot2::ggplot(cdb$cluster_tbl, ggplot2::aes(x = .data$n_cluster)) + ggplot2::geom_histogram() + ggplot2::xlab('Number of members') + ggplot2::ggtitle('Cluster sizes', subtitle =  cdb$cluster_pk))
    if(return_plotlist) return(plts)

    cowplot::plot_grid(plotlist = plts)
}
