#' @describeIn cluster_logistic_test split `ccdb` and conduct tests within strata
#' @param ... passed to `cluster_logistic_test`
#' @inheritParams split_cdb
#' @param ccdb [ContigCellDB()]
#' @export
cluster_test_by = function(ccdb, fields  = 'chain', tbl = 'cluster_tbl', ...){
    splat = split_cdb(ccdb, fields = fields, tbl = tbl, drop = TRUE)
    purrr::map_dfr(splat, function(x) cluster_logistic_test(...,  ccdb = x), .id = paste0(fields, collapse = '.'))
}

#' Test clusters for differential usage
#'
#' Typically one will want to stratify by chain by calling `cluster_test_by`, as this will calculate the number of cell "trials" separately depending on the chain recovered.
#' @param formula the **right-hand side** of a glmer or glm-style formula.
#' @param ccdb [ContigCellDB()]
#' @param replicate_keys keys naming columns in `ccdb$cell_tbl`
#' @param keep_fit `logical` as to whether the fit objects should be returned as a list column
#' @param fitter a function taking arguments `formula`, `data`, `is_mixed` and `keep_fit` that is run on each cluster.  Should return a `tibble` or `data.frame`
#' @param cluster_whitelist a table, keyed by `ccdb$cluster_pk` specifying the clusters to test. It does not alter which cells are included or how  canonicalization is performed. If omitted, all clusters will be tested and reported.
#' @inheritParams canonicalize_cell
#' @return table with one row per cluster/term.
#' @export
#'
#' @examples
#' data(ccdb_ex)
#' ccdb_ex = cluster_germline(ccdb_ex)
#' trav1 = dplyr::filter(ccdb_ex$cluster_tbl, v_gene == 'TRAV1')
#' cluster_logistic_test(~pop + (1|sample), ccdb_ex, 'sample', cluster_whitelist = trav1)
#' # Fixed effect analysis of each cluster, by chain
#' cluster_test_by(ccdb = ccdb_ex, fields = 'chain', tbl = 'cluster_tbl', formula = ~ pop, replicate_keys = 'sample')
cluster_logistic_test = function(formula, ccdb, replicate_keys, cluster_whitelist, contig_filter_args = TRUE, tie_break_keys = c('umis', 'reads'), keep_fit = FALSE, fitter = glm_glmer){
    if(!all(replicate_keys %in% names(ccdb$cell_tbl))) stop('Replicate must identify columns in `cell_tbl`.')
    if(length(ccdb$cluster_pk) == 0 || nrow(ccdb$cluster_tbl) == 0) stop("No clusters to test.")
    replicate_vec = do.call(interaction, ccdb$cell_tbl[replicate_keys])

    # canonicalize
    contig_fields = union(tie_break_keys, ccdb$cluster_pk)
    canon = canonicalize_cell(ccdb, contig_filter_args = !!rlang::enexpr(contig_filter_args), tie_break_keys, overwrite = TRUE, contig_fields = contig_fields)
    canon$cell_tbl[['_replicate_']] = replicate_vec

    # count clusters within replicate
    cluster_keys = union(ccdb$cluster_pk, replicate_keys)
    canon$cell_tbl = canon$cell_tbl %>% group_by(!!!syms(cluster_keys)) %>%
        mutate(x_ = dplyr::n()) %>% group_by(!!!syms(replicate_keys)) %>%
        mutate(n_ = dplyr::n())
    if(!missing(cluster_whitelist)){
        valid_KeyedTbl(cluster_whitelist, canon$cluster_pk)
        canon$cluster_tbl = semi_join(canon$cluster_tbl, cluster_whitelist, by = canon$cluster_pk)
    }
    n_cluster = nrow(canon$cluster_tbl)

    # formula/method munging
    rhs = paste0(as.character(formula), collapse = '')
    is_mixed = stringr::str_detect(rhs, stringr::fixed('|'))
    message("Fitting ", if(is_mixed) 'mixed logistic ' else 'fixed logistic ',
            "models to ", n_cluster, ' clusters.')

    formula = as.formula(stringr::str_c("cbind(x_, n_ - x_)", rhs))
    if(identical(fitter, glm_glmer)  && !requireNamespace('broom')){
        stop("Please install broom.")
    } else if(is_mixed && !requireNamespace('lme4')){
        stop("Please install lme4.")
    }
    safe_fit = purrr::possibly(fitter, tibble())
    canon$cell_tbl %>% group_by(!!!syms(ccdb$cluster_pk)) %>% do(safe_fit(formula, data = ., is_mixed, keep_fit))
}

glm_glmer = function(formula, data, is_mixed, keep_fit){
    func = if(is_mixed) lme4::glmer else stats::glm
    fit = func(formula, data = data, family = 'binomial')
    td = broom::tidy(fit, effects = 'fixed')
    if(keep_fit) td[['fit']] = list(fit)
    td
}
