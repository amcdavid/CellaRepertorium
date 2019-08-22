
#' Calculate number of cluster-subject singletons for the purposes of permutation testing
#'
#' @param cluster_idx factor-like cluster variable
#' @param subject factor-like subject
#'
#' @return average number of singletons
#' @export
#'
#' @seealso [cluster_permute_test()]
purity = function(cluster_idx, subject) {
    n_label_cluster = dplyr::bind_cols(cluster_idx = cluster_idx, subject = subject) %>%
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
#' @param statistic function of \code{covariates} (a `data.frame`), \code{labels} (a `vector`) that returns a \code{numeric} of length 1.
#' @param n_perm number of permutations.
#' @param ... passed to \code{statistic}
#' @param alternative `character` naming the direction `statistic` should be fall under the alternative hypothesis
#'
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
#' purity_ex = function(cluster_idx, subject) {
#' n_label_cluster = tibble(cluster_idx = cluster_idx, subject = subject) %>%
#' group_by(cluster_idx, subject) %>% summarize(n = dplyr::n()) %>% ungroup()
#' #Average number of singleton clusters for each subject
#' singletons = mean(n_label_cluster$n == 1)
#' singletons
#' }

#' cluster_permute_test(ccdb_test, 'subject', 'cluster_idx',
#' statistic = purity, n_perm  = 50)
#'

cluster_permute_test = function(ccdb, cell_covariate_keys,
                                cell_label_key = ccdb$cluster_pk, statistic,
                                n_perm,
                                alternative = c("two.sided", "less","greater"),
                                ...){
    if(!is.character(cell_label_key) || length(cell_label_key) > 1) stop('`cell_label_key` must name a single column in ccdb$cell_tbl')
    label = ccdb$cell_tbl[[cell_label_key]]
    if( length(na <- which(is.na(label))) > 0){
        warning('Excluding ', length(na), ' cells with missing labels.')
    }
    covariates = ccdb$cell_tbl[cell_covariate_keys]
    .cluster_permute_test(label, covariates, statistic, n_perm, alternative, ...)
}


#' Cell permutation tests (internal)
#'
#' @param labels \code{factor} of length n
#' @param covariates \code{factor} of length n
#' @param statistic function of label and covariate
#' @inheritParams cluster_permute_test
.cluster_permute_test = function(labels, covariates, statistic, n_perm, alternative,  ...){
    permp = rep(NA, n_perm)
    alternative = match.arg(alternative[1], c("two.sided", "less", "greater"), several.ok = FALSE)

    # checck that covariate nested within label
    observed = statistic(labels, covariates, ...)
    pb = progress::progress_bar$new(total = n_perm)
    for(i in seq_len(n_perm)){
        pb$tick()
        ci = sample(labels)
        permp[i] = statistic(ci, covariates, ...)
    }
    if(alternative == 'two.sided'){
        p_bound = min(mean(observed < permp), mean(observed > permp))*2
    } else{
        p_bound = if(alternative == 'less') mean(observed < permp) else mean(observed > permp)
    }
    list(observed = observed, expected = mean(permp), p.value = max(1/n_perm, p_bound), mc.se = sd(permp)/sqrt(n_perm))
}
