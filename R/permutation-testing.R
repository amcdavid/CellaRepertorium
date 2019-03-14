#' Tests for independence between clusters and labels using permutation
#'
#' This tests a statistic for association between clusters and labels by permuting the link between the two
#' Both clusters and labels should be factors.  Each observation represents a cell.
#' For example, \code{clusters} could be a clonal ID, while \code{labels} is a subject ID.
#' \code{statistic} is any
#'
#' @param clusters \code{factor} of length n
#' @param labels \code{factor} of length n
#' @param statistic function of \code{clusters}, \code{labels} and optionally \code{covariates} that returns a \code{numeric} of length 1.
#' @param n_perm number of permutations.
#' @param ... passed to \code{statistic}
#' @param covariates optional \code{factor} of length n.
#' @param alternative `character` naming the direction `statistic` should be fall under the alternative hypothesis
#'
#' @return a list containing the observed value of the statistic, its expectation (under independence), a p-value, and the monte carlo standard error (of the expected value).
#' @export
#'
#' @examples
#' library(dplyr)
#' purity = function(clusters, labels, covariates){
#' n_label_cluster = data_frame(labels = labels, clusters = clusters) %>% group_by(clusters, labels) %>% summarize(n = n()) %>% ungroup()
#' singletons = mean(n_label_cluster$n == 1)
#' singletons
#' }
#'
#' clusters = c(1, 1, 1, 2, 2, 3, 3)
#' labels = c('A', 'A', 'B', 'B', 'B', 'C', 'C')
#' covariates = c('X', 'X', 'Y', 'Y', 'Y', 'Y', 'Y')
#' cluster_permute_test(clusters, labels, statistic = purity, n_perm  = 50)
#'
cluster_permute_test = function(clusters, labels, covariates, statistic, n_perm = 1e3, alternative = c("two.sided", "less",
                                                                                                       "greater"),  ...){
    permp= rep(NA, n_perm)
    alternative = match.arg(alternative[1], c("two.sided", "less", "greater"), several.ok = FALSE)

    # checck that covariate nested within label
    if(!missing(covariates)){
        observed = statistic(clusters, labels, covariates, ...)
    } else{
        observed = statistic(clusters, labels, ...)
    }
    for(i in seq_len(n_perm)){
        ci = sample(clusters)
        permp[i] = if(!missing(covariates)) statistic(ci, labels, covariates,  ...) else statistic(ci, labels, covariates)
    }
    if(alternative == 'two.sided'){
        p_bound = min(mean(observed<permp), mean(observed>permp))*2
    } else{
        p_bound = if(alternative == 'less') mean(observed<permp) else mean(observed>permp)
    }
    list(observed = observed, expected = mean(permp), p.value = max(1/n_perm, p_bound), mc.se = sd(permp)/sqrt(n_perm))
}
