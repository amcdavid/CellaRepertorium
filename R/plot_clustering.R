#' Visualization of clustering characteristics for a ContigCellDB object with inherited fine_clustering object
#' 
#' @param ccdb A ContigCellDB object.
#' @param cluster_factors A character vector of variables in ccdb$contig table to visualize clusters by
#' @param type Type of visualization, a heatmap or a node-edge network plot
#' @param statistic Cluster characteristics visualized by pearson residuals or raw contig counts
#' @param ncluster Optional integer. Only plot characteristics that have at least nclusters (for clarity of visualization)
#' @param chaintype Character in ccdb$contig_tbl$chain. If passed will subset contigs belonging to specified chain (IGH,IGK,IGL,TRA,TRB)
#'
#' @return A ggraph object if type == 'network', and a ggplot object if type == 'heatmap'
#' @examples
#' data(ccdb_ex)
#' ccdb_germline_ex = cluster_germline(ccdb_ex, segment_keys = c('v_gene', 'j_gene', 'chain'), cluster_pk = 'segment_idx')
#' ccdb_germline_ex = fine_clustering(germline_cluster, sequence_key = 'cdr3_nt', type = 'DNA')
#' plot_clustering(ccdb_germline_ex,cluster_factors = c('v_gene','j_gene'), statistic = 'pearson', type = 'network' ,ncluster = 50, chaintype = 'IGH')

  
plot_clustering = function(ccdb,cluster_factors,type = c('heatmap','network'),statistic = c('pearson','contigs'),ncluster,chaintype) {

if(!inherits(ccdb,  "ContigCellDB")) stop('ccdb must have class ContigCellDb')

if(length(cluster_factors) != 2 | !is.character(cluster_factors)) stop('cluster_factors must be character vector that specifies exactly two columns from ccdb')
  
if(!missing(ncluster)) {tmp = ccdb$cluster_tbl %>% filter(n_cluster > ncluster)} else {tmp = ccdb$cluster_tbl}

if(!missing(chaintype)) {tmp = tmp %>% filter(chain == chaintype)} else {tmp = tmp}
  
if(statistic == 'pearson') {

  if(!missing(chaintype)) {ccdb$contig_tbl = ccdb$contig_tbl %>% filter(chain == chaintype)}
  
  resids = chisq.test(table(ccdb$contig_tbl %>% pull(cluster_factors[1]), ccdb$contig_tbl %>% pull(cluster_factors[2])))
  
  resid_tab = as.data.frame(resids$residuals)
  
  colnames(resid_tab) = c(cluster_factors,'residuals')
  
  tmp = left_join(tmp,resid_tab,by = cluster_factors)
  
  fillby = 'residuals'
  
}
  
else {fillby = 'n_cluster'}
  
tmp = tmp %>% relocate(cluster_factors)

if (type == 'network') {
  
  nodes = data.frame(names = c(unique(tmp %>% pull(cluster_factors[1])),unique(tmp %>% pull(cluster_factors[2]))))
  
  graph = graph_from_data_frame(tmp,vertices = nodes)

  plot = ggraph(graph,layout = 'circle') + geom_node_point() + geom_edge_fan(aes_string(colour = fillby)) + geom_node_label(aes(label = V(graph)$name), size = 2) + ggtitle(paste('Network Plot for factors', cluster_factors[1],'and',cluster_factors[2], 'using',statistic))
  
  if(statistic == 'pearson') {plot = plot + scale_edge_color_gradient2()}
  
  return(plot) }
  
else {
  plot = ggplot(tmp) + aes_string(x = cluster_factors[1],y = cluster_factors[2],fill = fillby) + geom_tile() + theme(axis.text.x = element_text(angle = 90)) +  ggtitle(paste('Heatmap for factors', cluster_factors[1],'and',cluster_factors[2], 'using',statistic))
  
  if (statistic == 'pearson') {plot = plot + scale_fill_gradient2()}

return(plot) }
}
  
