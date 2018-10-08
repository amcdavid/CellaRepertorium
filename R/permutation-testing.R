#' Tests for independence between clusters and labels using permutation
#'
#' This tests a statistic for association between clusters and labels by permuting the link between the two
#' Both clusters and labels should be factors.  Each observation represents a cell.
#' For example, \code{clusters} could be a clonal ID, while \code{labels} is a subject ID.
#' \code{statistic} is any
#' @param clusters \code{factor} of length n
#' @param labels \code{factor} of length n
#' @param statistic function of \code{clusters} and \code{labels} that returns a \code{numeric} of length 1.  Currently, a two-tailed test is always run.
#' @param n_perm number of permutations.
#' @param ... passed to \code{statistic}
#'
#' @return a list containing the observed value of the statistic, its expectation (under independence), a two-sided p-value, and the monte carlo standard error (of the expected value).
#' @export
#'
#' @examples
#' library(dplyr)
#' purity = function(clusters, labels){
#' n_label_cluster = data_frame(labels = labels, clusters = clusters) %>% group_by(clusters, si) %>% summarize(n = n()) %>% ungroup()
#' singletons = mean(n_label_cluster$n == 1)
#' singletons
#' }
#'
#' clusters = c(1, 1, 1, 2, 2, 3, 3)
#' cluster_permute_test
#'
cluster_permute_test = function(clusters, labels, statistic, n_perm = 1e3, ...){
    permp= rep(NA, n_perm)
    observed = statistic(clusters, labels, ...)
    for(i in seq_len(n_perm)){
        si = sample(labels)
        permp[i] = statistic(clusters, si, ...)
    }
    list(observed=observed, expected=mean(permp), p.value=max(1/n_perm, min(mean(observed<permp), mean(observed>permp))*2), mc.se=sd(permp)/sqrt(n_perm))
}
