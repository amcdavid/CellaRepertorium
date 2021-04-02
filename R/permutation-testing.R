.as.vector.tibble = function(x){
    if(inherits(x, 'data.frame')){
        if(ncol(x)>1) stop('Only supported for one column df')
        x[[1]]
    } else{
        as.vector(x)
    }
}


#' Calculate number of cluster-subject singletons for the purposes of permutation testing
#'
#' @param cluster_idx factor-like cluster variable
#' @param subject factor-like subject
#'
#' @return average number of singletons
#' @export
#' @examples
#' message("see example(cluster_permute_test)")
#' @seealso [cluster_permute_test()]
purity = function(cluster_idx, subject) {
    cluster_idx = .as.vector.tibble(cluster_idx)
    subject = .as.vector.tibble(subject)
    n_label_cluster = dplyr::tibble(cluster_idx = cluster_idx, subject = subject) %>%
        group_by(cluster_idx, subject) %>% summarize(n = dplyr::n()) %>% ungroup()
    #Average number of singleton clusters for each subject
    singletons = mean(n_label_cluster$n == 1)
    singletons
}


#' Tests for independence between labels and covariates using permutation of cells
#'
#' This tests a statistic for association between labels (for instance, cluster/clonal ID) and covariates (for instance, subject or treatment) by permuting the link between the two.
#' Each observation represents a cell.
#' \code{statistic} is any function of `labels`
#'
#' @param ccdb `ContigCellDB`
#' @param cell_covariate_keys  `character` naming fields in `ccdb$cell_tbl`
#' @param cell_label_key `character` naming a single field in `ccdb$cell_tbl`
#' @param cell_stratify_keys optional `character` naming fields in `ccdb$cell_tbl` under which permutations of `cell_label_key` will occur.
#' This means that the test will occur conditional on these covariates.
#' Must be disjoint from  `cell_covariate_keys`.
#' @param sanity_check_strata `logical`, should `cell_stratify_keys` be checked for sanity?
#' @param ... passed to \code{statistic}
#' @inheritParams .cluster_permute_test
#' @return a list containing the observed value of the statistic, its expectation (under independence), a p-value, and the Monte Carlo standard error (of the expected value).
#' @seealso [purity()]
#' @export
#'
#' @examples
#' library(dplyr)
#' # covariate should name one or more columns in `cell_tbl`
#'
#' cluster_idx = c(1, 1, 1, 2, 2, 3, 3)
#' subject = c('A', 'A', 'B', 'B', 'B', 'C', 'C')
#' contig_tbl = tibble(contig_pk = seq_along(cluster_idx), cluster_idx, subject)
#' ccdb_test = ContigCellDB(contig_tbl = contig_tbl, contig_pk = 'contig_pk',
#' cell_pk = c('contig_pk', 'subject', 'cluster_idx'), cluster_pk = 'cluster_idx')
#' ccdb_test$cell_tbl
#'
#' cluster_permute_test(ccdb_test, 'subject', 'cluster_idx',
#' statistic = purity, n_perm  = 50)
#'

cluster_permute_test = function(ccdb, cell_covariate_keys,
                                cell_label_key = ccdb$cluster_pk,
                                cell_stratify_keys, statistic,
                                n_perm,
                                alternative = c("two.sided", "less","greater"),
                                sanity_check_strata = TRUE,
                                ...){
    if(!is.character(cell_label_key) || length(cell_label_key) > 1)
        stop('`cell_label_key` must name a single column in ccdb$cell_tbl')
    label = ccdb$cell_tbl[[cell_label_key]]
    if( length(na <- which(is.na(label))) > 0){
        warning('Excluding ', length(na), ' cells with missing labels.')
    }
    covariates = ccdb$cell_tbl[cell_covariate_keys]
    if(!missing(cell_stratify_keys)){
        if(length(overlap <- intersect(cell_stratify_keys, cell_covariate_keys))>0)
            stop("`cell_stratify_keys` must be disjoint from `cell_covariate_keys`.")
        strata = do.call(interaction, ccdb$cell_tbl[cell_stratify_keys])
        if(sanity_check_strata && length(unique(strata)) > .75*length(unique(label))) stop("Lots of strata, set `sanity_check_strata = FALSE` if you are sure.")
    } else{
        strata  = NULL
    }
    .cluster_permute_test(label, covariates, strata, statistic, n_perm, alternative, ...)
}


#' Cell permutation tests (internal)
#'
#' @param labels \code{factor} of length n
#' @param covariates \code{data.frame} of length n
#' @param strata \code{factor}
#' @param statistic function of label (vector) and covariate (data.frame). Must return a scalar
#' @param alternative `character` naming the direction `statistic` should be fall under the alternative hypothesis
#' @param n_perm number of permutations to run
#' @param ... passed along to `statistic`
#' @return a list containing the observed value of the statistic, its expectation (under independence), a p-value, and the Monte Carlo standard error (of the expected value).
.cluster_permute_test = function(labels, covariates, strata, statistic, n_perm, alternative,  ...){
    permp = rep(NA, n_perm)
    alternative = match.arg(alternative[1], c("two.sided", "less", "greater"), several.ok = FALSE)

    observed = statistic(labels, covariates, ...)
    pb = progress::progress_bar$new(total = n_perm)
    if(is.null(strata)){
        sampler = sample
    } else{
        sampler = function(labels){
                unsplit(lapply(split(labels, strata), sample), strata)
        }
    }
    for(i in seq_len(n_perm)){
        pb$tick()
        ci = sampler(labels)
        permp[i] = statistic(ci, covariates, ...)
    }
    if(alternative == 'two.sided'){
        p_bound = min(mean(observed < permp), mean(observed > permp))*2
    } else{
        p_bound = if(alternative == 'less') mean(observed < permp) else mean(observed > permp)
    }
    list(observed = observed, expected = mean(permp),
         p.value = max(1/n_perm, p_bound), mc.se = sd(permp)/sqrt(n_perm))
}
