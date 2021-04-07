require_igraph = function(){
  if(!requireNamespace('igraph') || !requireNamespace('ggraph')) stop('Install `igraph` and `ggraph`.')
  TRUE
}

#' Visualization of pairs of cluster factor
#'
#' With `factors`, a pair of variables present in the `contig_tbl` and the `cluster_tbl`,
#' generate and plot cross-tabs of the number of contigs, or its pearson residual.
#' @param ccdb A ContigCellDB object.
#' @param factors `character` length 2 of fields present
#' @param type Type of visualization, a heatmap or a node-edge network plot
#' @param statistic Cluster characteristics visualized by pearson residuals or raw contig counts
#' @param ncluster `integer`.  Omit factors that occur less than `nclusters`. For clarity of visualization.
#' @param chaintype Character in ccdb$contig_tbl$chain. If passed will subset contigs belonging to specified chain (IGH,IGK,IGL,TRA,TRB)
#' @seealso canonicalize_cluster to "roll-up" additional contig variables into the `cluster_tbl``
#' @return A ggraph object if type == 'network', and a ggplot object if type == 'heatmap'
#' @export
#' @examples
#' library(ggraph)
#' data(ccdb_ex)
#' ccdb_germline_ex = cluster_germline(ccdb_ex, segment_keys = c('v_gene', 'j_gene', 'chain'),
#' cluster_pk = 'segment_idx')
#' ccdb_germline_ex = fine_clustering(ccdb_germline_ex, sequence_key = 'cdr3_nt', type = 'DNA')
#' plot_cluster_factors(ccdb_germline_ex,factors = c('v_gene','j_gene'),
#' statistic = 'pearson', type = 'network' ,ncluster = 10, chaintype = 'TRB')
#' plot_cluster_factors(ccdb_germline_ex,factors = c('v_gene','j_gene'),
#' statistic = 'contigs', type = 'heatmap')
#' plot_cluster_factors(ccdb_germline_ex,factors = c('v_gene','j_gene'),
#' statistic = 'contigs', type = 'network', ncluster = 10)
plot_cluster_factors = function(ccdb,
                           factors,
                           type = c('heatmap', 'network'),
                           statistic = c('pearson', 'contigs'),
                           ncluster = 0,
                           chaintype) {
  if (!inherits(ccdb,  "ContigCellDB") || !has_fineclustering(ccdb))
    stop('ccdb must be  `ContigCellDB` with `fine_clustering`.')
  check_plot_infra()
  statistic = match.arg(statistic, c('pearson', 'contigs'))
  type = match.arg(type, c('heatmap', 'network'))
  if (length(factors) != 2 ||
      !is.character(factors))
    stop('`factors` must be character vector that specifies exactly two columns from ccdb')

  cluster_tbl =  dplyr::filter(ccdb$cluster_tbl, n_cluster > ncluster)


  if (!missing(chaintype)) {
    cluster_tbl = cluster_tbl %>% filter(chain %in% chaintype)
  }

  if(!nrow(cluster_tbl) > 0) stop("No clusters remain after filtering.")

  if (statistic == 'pearson') {
    if (!missing(chaintype)) {
      ccdb$contig_tbl = ccdb$contig_tbl %>% filter(chain %in% chaintype)
    }
    if(!nrow(ccdb$contig_tbl)>0) stop("No contigs left after filtering.")
    resids = stats::chisq.test(table(
      ccdb$contig_tbl[[factors[1]]],
      ccdb$contig_tbl[[factors[2]]]
    ))
    resid_tab = as.data.frame(resids$residuals)
    colnames(resid_tab) = c(factors, 'residuals')
    cluster_tbl = left_join_warn(cluster_tbl, resid_tab, by = factors)
    fillby = 'residuals'
  } else {
    fillby = 'n_cluster'
  }
  cluster_tbl = cluster_tbl %>% dplyr::relocate(factors)

  if (type == 'network') {
    require_igraph()
    nodes = data.frame(names = c(unique(cluster_tbl[[factors[1]]]),unique(cluster_tbl[[factors[2]]])))
    graph = igraph::graph_from_data_frame(cluster_tbl, vertices = nodes)

    plot = ggraph::ggraph(graph, layout = 'circle') + ggraph::geom_node_point() + ggraph::geom_edge_fan(ggplot2::aes_string(colour = fillby)) + ggraph::geom_node_label(ggplot2::aes(label = igraph::V(graph)$name), size = 2) + ggplot2::ggtitle(
      paste(
        'Network Plot for factors',
        factors[1],
        'and',
        factors[2],
        'using',
        statistic
      )
    )

    if (statistic == 'pearson') {
      plot = plot + ggraph::scale_edge_color_gradient2()
    }

    return(plot)
  }

  else {
    plot = ggplot2::ggplot(cluster_tbl) + ggplot2::aes_string(x = factors[1], y = factors[2], fill = fillby) + ggplot2::geom_tile() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +  ggplot2::ggtitle(
      paste(
        'Heatmap for factors',
        factors[1],
        'and',
        factors[2],
        'using',
        statistic
      )
    )

    if (statistic == 'pearson') {
      plot = plot + ggplot2::scale_fill_gradient2()
    }

    return(plot)
  }
}

