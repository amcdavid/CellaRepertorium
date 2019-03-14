#' Generate a legible name for a series of contigs
#'
#' @param contig_frame An `all_contig_annotations.csv` file, output from VDJ Cell ranger.  Importantly, this should contain columns `chain`, `v_gene`, `d_gene`, `j_gene`
#' @param prefix an optional prefix added to each contig, eg, possibly a sample id.
#' @return \code{character}
#' @importFrom dplyr %>%
#' @importFrom stringr str_replace_all
#' @importFrom methods as
#' @importFrom stats as.dist dist hclust na.fail sd
#' @importFrom utils data
#' @export
#'
#' @examples
#' library(dplyr)
#' contig_anno_path = system.file('extdata', 'cellranger_contig_annotation.csv', package = 'CellaRepertorium')
#' contig_anno = readr::read_csv(contig_anno_path)
#' contig_anno = contig_anno %>% mutate(fancy_name = fancy_name_contigs(., prefix = paste(sample, pop, sep = '_')))
#' stopifnot(!any(duplicated(contig_anno$fancy_name)))
fancy_name_contigs = function(contig_frame, prefix){
    notin = setdiff(c('chain', 'v_gene', 'd_gene', 'j_gene'), names(contig_frame))
    if(length(notin>0)) stop("`contig_frame must contain all of ", paste(notin, collapse = ','))
    chain = contig_frame$chain
paste(contig_frame$v_gene, contig_frame$d_gene, contig_frame$j_gene, sep = ':') %>% str_replace_all('None', '') %>% str_replace_all('IG[KLH]|TR[ABDG]', '') %>% paste(chain, ., sep = ':') %>% paste(prefix, ., sep = '_') %>% make.unique()

}

variable_genes = function(ref = '10X'){
    # This should return genes that lie in the VDJ region for a given annotation
    # And should know something about the chromium reference DB.

}


#' Generate a list of tables representing clusters paired in cells
#'
#' A contingency table of every combination of `cluster_idx` up to `table_order` is generated.
#' Combinations that are found in at least `min_expansion` number of cells are reported.  All cells that have these combinations are returned, as well as cells that only have `orphan_level` of matching `cluster_idx`.
#'
#' For example, if `table_order=2` and `min_expansion=2` then heavy/light or alpha/beta pairs found two or more times will be returned (as well as alpha-alpha pairs, etc, if those are present).
#' If `orphan_level=1` then all cells that share just a single chain with an expanded clone will be returned.
#'
#' The `cluster_idx.1_fct` and `cluster_idx.2_fct` fields in `cell_tbl`, `idx1_tbl`, `idx2_tbl` are cast to factors and ordered such that pairings will tend to occur along the diagonal when they are cross-tabulated.
#' This facilitates plotting.
#'
#' @section Caveats and warnings:
#'  The cell_idx -> cluster_idx map is generally one-to-many.  For `table_order>1`, few collisions are expected as most cells will contain no more than 2 clusters.  Any collisions are resolved by returning the most prevalent cluster, across samples.  When two clusters are tied for most prevalent within a cell, the `cluster_idx` returned is arbitrary. Therefore, when `table_order=1`, it is strongly recommended to subset the `cluster_tbl` to just a single chain.
#' @param cluster_tbl a table with all combinations of clusters in all cells
#' @param cell_identifiers character vector naming fields that key a cell
#' @param cluster_idx character naming a single field IDing the clusters
#' @param min_expansion the minimal number of times a pairing needs to occur for it to be reported
#' @param cell_tbl optional, ancillary table with additional cell features.  Must also be keyed by `cell_identifiers`
#' @param feature_tbl optional, ancillary table with additional cluster features.  Must also be keyed by `cluster_idx`
#' @param table_order Integer larger than 1. What order of cluster_idx will be paired, eg, order = 2 means that the most common and second most common cluster_idx will be sought for each cell
#' @param orphan_level Integer larger than 0 and less than or equal to `table_order`.  Given that at least `min_expansion` cells are found that have `table_order` chains identical, how many `cluster_idx` pairs will we match on to select other cells.  Example: `ophan_level=1` means that cells that share just a single chain with the
#' @param cluster_whitelist a table of "cluster_idx" that should always be reported.  In contrast to the `cluster_tbl`, here the clusters must be named "cluster_idx.1", "cluster_idx.2" (if order-2 pairs are being selected).
#' @param cluster_blacklist a table of "cluster_idx" that will never be reported.  Must be named as per `cluster_whitelist`.
#'
#' @return list of tables.  The `cell_tbl` is keyed by the `cell_identifiers`, with fields "cluster_idx.1", "cluster_idx.2", etc, IDing the contigs present in each cell. "cluster_idx.1_fct" and "cluster_idx.2_fct" cast these fields to factors and are reordered to maximize the number of pairs along the diagonal. The `idx1_tbl` and `idx2_tbl` report information (passed in about the `cluster_idx` by `feature_tbl`.)  The `cluster_pair_tbl` reports all pairings found of contigs, and the number of times observed.
#' @export
#'
#' @seealso get_canonical_chain
#' @importFrom tibble as_data_frame
#' @importFrom dplyr bind_rows left_join ungroup summarize anti_join
#' @importFrom stringr str_length str_c
#' @importFrom rlang sym syms :=
#' @examples
#' library(dplyr)
#' cluster_tbl = data_frame(clust_idx = gl(3, 2), cell_idx = rep(1:3, times = 2))
#' # no pairs found twice
#' pairing_tables(cluster_tbl, 'cell_idx', 'clust_idx')
#' # all pairs found, found once.
#' pairing_tables(cluster_tbl, 'cell_idx', 'clust_idx', min_expansion = 1)
#' cluster_tbl2 = bind_rows(cluster_tbl, cluster_tbl %>% mutate(cell_idx = rep(4:6, times = 2)))
#' #all pairs found twice
#' pairing_tables(cluster_tbl2, 'cell_idx', 'clust_idx', min_expansion = 1)
#' # A warning if `table_order = 1` -- best to stratify by chain
#' # (could be accomplished by setting the `cell_identifiers` to include a field that identifies the chain).
#' pairing_tables(cluster_tbl2, 'cell_idx', 'clust_idx', min_expansion = 1, table_order = 1)
pairing_tables = function(cluster_tbl, cell_identifiers = 'barcode', cluster_idx = 'cluster_idx', table_order = 2, min_expansion = 2,  orphan_level = 1, cluster_whitelist = NULL, cluster_blacklist = NULL, cell_tbl = NULL, feature_tbl = NULL ){

    if(orphan_level > table_order) stop('`ophan_level` must be less than or equal to `table_order`')
    if(table_order < 1) stop('Table order must be at least 1')

    # get `table_order` most common clusters for each cell
    # forcibly rename cluster_idx -> "cluster_idx"
    bar_chain_tbls = purrr::map(seq_len(table_order), function(i) get_canonical_chain(cluster_tbl, cell_identifiers = cell_identifiers, cluster_idx = cluster_idx, order = i) %>% dplyr::rename( !!paste0('cluster_idx.', i) := !!cluster_idx))
    # for each cell, what clusters are present
    oligo_cluster_pairs = purrr::reduce(bar_chain_tbls, left_join, by = cell_identifiers)

    # Set up cluster_ids for indexing into the pairing tables
    cluster_ids = str_c('cluster_idx.', seq_len(table_order))
    cluster_ids_to_select = cluster_ids[seq_len(orphan_level)]

    # In how many cells do each cluster pairing appear?
    cluster_pair_tbl = oligo_cluster_pairs %>% group_by(!!!syms(cluster_ids)) %>% summarize(n_clone_pairs = n())
    # which clusters are expanded
    expanded_cluster = cluster_pair_tbl %>% dplyr::filter(n_clone_pairs >= min_expansion) %>% dplyr::filter_at(.vars = cluster_ids_to_select, .vars_predicate = all_vars(!is.na(.)))
    expanded_cluster = ungroup(expanded_cluster) %>% dplyr::select(!!!syms(cluster_ids_to_select), max_pairs = n_clone_pairs)
    if(!is.null(cluster_whitelist)){
        expanded_cluster = bind_rows(expanded_cluster, cluster_whitelist)
    }
    if(!is.null(cluster_blacklist)) expanded_cluster = anti_join(expanded_cluster, cluster_blacklist)

    # Could have duplicated cluster_ids after binding to the whitelist or from considering orphans
    expanded_cluster = expanded_cluster[!duplicated(expanded_cluster %>% dplyr::select(-max_pairs)),]
    expanded_c1 = oligo_cluster_pairs %>% right_join(expanded_cluster, by = cluster_ids_to_select)
    if(anyDuplicated(expanded_c1 %>% dplyr::select(!!!syms(cell_identifiers)))) stop("Ruhoh, duplicated cell identifiers, this is a bug!")

    # cross-tab pairings in order to figure out how to order the cluster_idx to put most common pairing on diagonal
    if(nrow(expanded_c1)>0 && table_order > 1){
        expanded_counts = reshape2::dcast(expanded_c1, cluster_idx.1 ~ cluster_idx.2, fun.aggregate = length, value.var = 'max_pairs')
        mat = expanded_counts[,-1]
        #rownames(mat) = rowid[['cluster_idx']]
        ro = hclust(dist(mat))$order
        co = hclust(dist(t(mat)))$order
    } else if (nrow(expanded_c1) == 0) {
       warning('No pairs found')
        return()
    } else{
            expanded_counts = reshape2::dcast(expanded_c1, cluster_idx.1 ~ ., fun.aggregate = length, value.var = 'max_pairs')
            ro = order(expanded_counts[[2]])
            co = NA_integer_
    }

    ci_class = class(cluster_tbl[[cluster_idx]])
    as_method = if(ci_class == 'factor') as.factor else function(x) as(x, ci_class)
    rowid = data_frame(cluster_idx = expanded_counts$cluster_idx.1 %>% as_method, plot_order = ro)
    colid = suppressWarnings(data_frame(cluster_idx = colnames(expanded_counts)[-1] %>% as_method, plot_order = co))
    rowid[['cluster_idx.1_fct']] = factor(rowid[['cluster_idx']], levels = rowid[['cluster_idx']][ro])
    colid[['cluster_idx.2_fct']] = factor(colid[['cluster_idx']], levels = colid[['cluster_idx']][co])


    # also fix levels of this tbl
    expanded_c1 = expanded_c1 %>% mutate(cluster_idx.1_fct = factor(cluster_idx.1, levels = levels(rowid[['cluster_idx.1_fct']])))
    if(table_order>1) expanded_c1 = expanded_c1 %>% mutate(cluster_idx.2_fct = factor(cluster_idx.2, levels = levels(colid[['cluster_idx.2_fct']])))

    if(!is.null(cell_tbl)){
        if(anyDuplicated(cell_tbl %>% select(!!!syms(cell_identifiers)))) stop('`cell_tbl` must not have duplicate `cell_identifiers`.')
        expanded_c1 = left_join(expanded_c1, cell_tbl, by = cell_identifiers)
    }

    idx1_tbl = rowid %>% dplyr::rename(!!cluster_idx := cluster_idx)
    idx2_tbl = colid %>% dplyr::rename(!!cluster_idx := cluster_idx)
    if(!is.null(feature_tbl)){
        if(anyDuplicated(feature_tbl %>% select(!!sym(cluster_idx)))) stop('`feature_tbl` must not have duplicate `cluster_idx`.')
        idx1_tbl = left_join(idx1_tbl, feature_tbl, by = cluster_idx)
        idx2_tbl = left_join(idx2_tbl, feature_tbl, by = cluster_idx)
    }

    list(cell_tbl = expanded_c1, idx1_tbl = idx1_tbl, idx2_tbl = idx2_tbl, cluster_pair_tbl = cluster_pair_tbl)

}

plot_pairing = function(pairing_list, color_labels_by){
    pl = pairing_list
    pairs_plt = ggplot(pairing_list$cell_tbl, aes(x = cluster_idx.1_fct, y = cluster_idx.2_fct, color = sample, shape = pop)) + geom_jitter(width = .3, height = .3)

    ylab = data_frame(!!color_labels_by :=  ggplot_build(pairs_plt)$layout$panel_params[[1]]$y.label) %>% left_join(feature_tbl) %>% mutate(class_color = ifelse(is.na(class_color), '#E41A1C', class_color))

    xlab = data_frame(!!color_labels_by :=  ggplot_build(pairs_plt)$layout$panel_params[[1]]$x.label) %>% left_join(feature_tbl) %>% mutate(class_color = ifelse(is.na(class_color), '#E41A1C', class_color))

    pairs_plt = pairs_plt + theme(axis.text.x = element_text(angle = 90, color = xlab$class_color, size = 8), axis.text.y = element_text(color = ylab$class_color, size = 8))


}

cleanup_annotations = function(json){
    anno = json[['annotations']]
    atomic = purrr::map(anno, ~ .[setdiff(names(.), 'mismatches')])
    names(atomic) = json[['contig_name']]
    as_data_frame(bind_rows(atomic, .id = 'contig_name'))
}

get_gaps = function(ca){
    vregion = which(ca[['feature.region_type']] == 'L-REGION+V-REGION')
    dregion = which(ca[['feature.region_type']] == 'D-REGION')
    jregion = which(ca[['feature.region_type']] == 'J-REGION')
    vj_gap <- dj_gap <- vd_gap <- NA_integer_
    if(length(vregion) > 0 & length(dregion) > 0) vd_gap = ca[[dregion, 'contig_match_start']] - ca[[vregion,'contig_match_end']]
    if(length(dregion) > 0 & length(jregion) > 0) dj_gap = ca[[jregion, 'contig_match_start']] - ca[[dregion,'contig_match_end']]
    if(length(vregion) > 0 & length(jregion) > 0) vj_gap = ca[[jregion,'contig_match_start']] - ca[[vregion,'contig_match_end']]
    data_frame(vd_gap, dj_gap, vj_gap)
}

read_contig_json = function(file, seq_cols = c('quals', 'aa_sequence', '')){
    jsn = fromJSON(file(anno_file), flatten = TRUE)
    # contig sequences and gaps
}
