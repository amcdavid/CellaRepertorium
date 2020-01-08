#' A filtration of clusters
#'
#' Return clusters that match all provided conditions
#' @param min_number `integer` At least this many cells
#' @param min_freq `numeric` At least this frequency
#' @param white_list `data.frame` keyed by cluster_pk that must match
#'
#' @return object representing the filtration (currently a list)
#' @export
#'
#' @examples
#' cluster_filterset(min_number = 1, min_freq = 0)
cluster_filterset = function(min_number = 0, min_freq = 0, white_list = NULL){
    if(!is.numeric(min_number)||length(min_number) > 1 || floor(min_number) != min_number) stop("`min_number` must be scalar integer.")
    if(!is.numeric(min_freq)|| length(min_freq) > 1 || min_freq < 0 || min_freq > 1) stop("`min_freq` must be scalar between 0 and 1.")
    if(!is.null(white_list) && !inherits(white_list, 'data.frame')) stop("`white_list` must be data.frame.")

    invisible(list(min_number = min_number, min_freq = min_freq, white_list = white_list))
}

.execute_filter = function(ccdb, filterset, tie_break_keys = character(), contig_filter_args = TRUE){
    # canonicalize
    contig_fields = union(tie_break_keys, ccdb$cluster_pk)
    # should this warning be hushed?
    canon = canonicalize_cell(ccdb,
                              contig_filter_args = !!rlang::enexpr(contig_filter_args),
                              tie_break_keys, overwrite = TRUE, contig_fields = contig_fields)
    count = canon$cell_tbl %>% group_by(!!!syms(ccdb$cluster_pk)) %>% summarize(n = n(), freq = n/nrow(canon$cell_tbl))
    if(!is.null(filterset$white_list)) count = semi_join(count, filterset$white_list, by = ccdb$cluster_pk)
    filter(count, n >= filterset$min_number, freq >= filterset$min_freq)
}

#' @describeIn cluster_logistic_test split `ccdb` and conduct tests within strata
#' @param ... passed to `cluster_logistic_test`
#' @inheritParams split_cdb
#' @export
cluster_test_by = function(ccdb, fields  = 'chain', tbl = 'cluster_tbl', ...){
    splat = split_cdb(ccdb, fields = fields, tbl = tbl, drop = TRUE, equalize = TRUE)
    res = purrr::map_dfr(splat, function(x) cluster_logistic_test(...,  ccdb = x), .id = paste0(fields, collapse = '.'))
    # in case .id was included in something else we joined to
    res[,unique(names(res))]
}

#' Test clusters for differential usage
#'
#' Typically one will want to stratify by chain by calling `cluster_test_by`, as this will calculate the number of cell "trials" separately depending on the chain recovered.
#' @param formula the **right-hand side** of a glmer or glm-style formula.
#' @param ccdb [ContigCellDB()]
#' @param keep_fit `logical` as to whether the fit objects should be returned as a list column
#' @param fitter a function taking arguments `formula`, `data`, `is_mixed` and `keep_fit` that is run on each cluster.  Should return a `tibble` or `data.frame`
#' @param silent `logical`. Should warnings from fitting functions should be suppressed?
#' @param filterset a call to [cluster_filterset()] that will be used to subset clusters.
#' @param add_cluster_tbl `logical` should the output be joined to the `cluster_tbl`?
#' @inheritParams canonicalize_cell
#' @return table with one row per cluster/term.
#' @export
#' @importFrom stats as.formula
#'
#' @examples
#' library(dplyr)
#' data(ccdb_ex)
#' ccdb_ex = cluster_germline(ccdb_ex)
#' trav1 = filter(ccdb_ex$cluster_tbl, v_gene == 'TRAV1')
#' cluster_logistic_test(~pop + (1|sample), ccdb_ex, cluster_whitelist = trav1)
#' # Fixed effect analysis of each cluster, by chain
#' prev4 = ccdb_ex$contig_tbl %>% group_by(cluster_idx) %>%
#' summarize(n()) %>% filter(`n()`>= 4)
#' cluster_test_by(ccdb = ccdb_ex, fields = 'chain',
#' tbl = 'cluster_tbl', formula = ~ pop, cluster_whitelist = prev4)
cluster_logistic_test = function(formula, ccdb, filterset = cluster_filterset(),
                                 contig_filter_args = TRUE, tie_break_keys = c('umis', 'reads'), add_cluster_tbl = FALSE,
                                 keep_fit = FALSE, fitter = glm_glmer, silent = FALSE){
    if(length(ccdb$cluster_pk) == 0 || nrow(ccdb$cluster_tbl) == 0) stop("No clusters to test.")

    # canonicalize
    contig_fields = union(tie_break_keys, ccdb$cluster_pk)
    # should this warning be hushed?
    canon = canonicalize_cell(ccdb,
                              contig_filter_args = !!rlang::enexpr(contig_filter_args),
                              tie_break_keys, overwrite = TRUE, contig_fields = contig_fields)

    keep = .execute_filter(ccdb, filterset, tie_break_keys, !!rlang::enexpr(contig_filter_args))

    # Only test clusters that are present
    cluster_whitelist = semi_join(keep, canon$cluster_tbl, by = canon$cluster_pk)
    # make a uniquely named column that enumerates clusters
    cluster_idx_nm = make.unique(c('clidx', union(names(cluster_whitelist), names(canon$cell_tbl))))[1]
    n_cluster = nrow(cluster_whitelist)
    cluster_whitelist[[cluster_idx_nm]] = seq_len(nrow(cluster_whitelist))
    canon$cell_tbl = left_join_warn(canon$cell_tbl, cluster_whitelist[union(canon$cluster_pk, cluster_idx_nm)], by = canon$cluster_pk)
    canon$cell_tbl[[cluster_idx_nm]] = ifelse(is.na(canon$cell_tbl[[cluster_idx_nm]]), -1, canon$cell_tbl[[cluster_idx_nm]])
    x_nm = make.unique(c('x_', names(canon$cell_tbl)))[1]

    # formula/method munging
    rhs = paste0(as.character(formula), collapse = '')
    is_mixed = stringr::str_detect(rhs, stringr::fixed('|'))
    message("Fitting ", if(is_mixed) 'mixed logistic ' else 'fixed logistic ',
            "models to ", n_cluster, ' clusters.')

    formula = as.formula(stringr::str_c(x_nm, rhs))

    if(identical(fitter, glm_glmer)  && !requireNamespace('broom')){
        stop("Please install broom.")
    } else if(is_mixed && !requireNamespace('lme4')){
        stop("Please install lme4.")
    }
    safe_fit = purrr::possibly(fitter, tibble())
    pb = progress::progress_bar$new(total = n_cluster)
    res = purrr::map_dfr(cluster_whitelist[[cluster_idx_nm]], function(cluster_idx){
        pb$tick()
        canon$cell_tbl[[x_nm]] = canon$cell_tbl[[cluster_idx_nm]] == cluster_idx
        if(silent) suppressWarnings(safe_fit(formula, canon$cell_tbl, is_mixed, keep_fit)) else safe_fit(formula, canon$cell_tbl, is_mixed, keep_fit)
    }, .id = cluster_idx_nm)
    res[[cluster_idx_nm]] = as.integer(res[[cluster_idx_nm]])
    res = left_join(res, cluster_whitelist[c(canon$cluster_pk, cluster_idx_nm)], by = cluster_idx_nm) %>% select(-!!sym(cluster_idx_nm))
    if(add_cluster_tbl) res = left_join_warn(res, ccdb$cluster_tbl, by = ccdb$cluster_pk)
    return(res)
}

glm_glmer = function(formula, data, is_mixed, keep_fit){
    func = if(is_mixed) lme4::glmer else stats::glm
    fit = func(formula, data = data, family = 'binomial')
    td = broom::tidy(fit, effects = 'fixed')
    if(keep_fit) td[['fit']] = list(fit)
    td
}
