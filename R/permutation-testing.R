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
#' @return a list containing the observed value of the statistic, the permuted values of the statistic, its expectation (under independence), a p-value, and the Monte Carlo standard error (of the expected value).
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
#' clust_test = cluster_permute_test(ccdb_test, 'subject', 'cluster_idx',
#' statistic = purity, n_perm  = 50)
#' library(ggplot2)
#' plot_permute_test(perm_test = clust_test)
cluster_permute_test = function(ccdb, cell_covariate_keys,
                                cell_label_key = ccdb$cluster_pk,
                                cell_stratify_keys, statistic,
                                contrasts = NULL, n_perm,
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
    .cluster_permute_test(label, covariates, strata, statistic, contrasts, n_perm, alternative, ...)
}

pairs = function (n, names) {
  m = matrix(0, nrow = n*(n-1)/2, ncol = n)
  ii = 1
  for(i in seq_len(n - 1)){
    for(j in seq(from = i + 1, to = n)){
      m[ii, c(i,j)] = c(1, -1)
      ii = ii + 1
    }
  }
  colnames(m) = names
  m
}

#' Cell permutation tests (internal)
#'
#' @param labels \code{factor} of length n
#' @param covariates \code{data.frame} of length n
#' @param strata \code{factor}
#' @param statistic function of label (vector) and covariate (data.frame). 
#' If this returns a vector, then by default each level will be compared against each other, pairwise, but see the next section.
#' @param contrasts an optional list of numeric vectors. 
#' Each will be dotted with the statistic, or optionally a matrix provided in which case each **row** would be tested one-by-one.
#' @param alternative `character` naming the direction `statistic` should be fall under the alternative hypothesis
#' @param n_perm number of permutations to run
#' @param ... passed along to `statistic`
#' @return a list containing the observed value of the statistic, the permuted values of the statistic, its expectation (under independence), a p-value, and the Monte Carlo standard error (of the expected value).
.cluster_permute_test = function(labels, covariates, strata, statistic, contrasts, n_perm, alternative,  ...){
    cl = match.call()
    observed = statistic(labels, covariates, ...)
    len_stat = length(observed)
    permp = matrix(NA_real_, nrow = n_perm, ncol = len_stat)
    
    # Wrangle contrasts
    # empty
    if(is.null(contrasts)) contrasts = pairs(len_stat, names = names(observed))
    # list -> matrix
    if(is.list(contrasts)) contrasts = do.call(rbind, contrasts)
    # coerce numeric non-matrices
    if(is.numeric(contrasts) && !is.matrix(contrasts)) dim(contrasts) = c(1, len_stat)
    # check for correct structure
    if(!is.matrix(contrasts) || !is.numeric(contrasts) || ncol(contrasts) != len_stat){
      stop("`contrasts` must be numeric matrix with ", len_stat, " columns.")
    }
    
    # finally back to a list 
    # (this seemed like the easiest way to check that a user-provided list has correct shape..?)
    contrasts = split(contrasts, gl(nrow(contrasts), 1, length(contrasts)))
    
    alternative = match.arg(alternative[1], c("two.sided", "less", "greater"), several.ok = FALSE)
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
        permp[i,] = statistic(ci, covariates, ...)
    }
    
    permute_tail_probs = function(x, obs, contrast) {
      if(alternative == 'two.sided'){
        p_bound = min(mean(obs < x), mean(obs > x))*2
      } else{
        p_bound = if(alternative == 'less') mean(obs < x) else mean(obs > x)
      }
      out = list(observed = obs, expected = mean(x),
                 p.value = max(1/n_perm, p_bound), 
                 mc.se = sd(x)/sqrt(n_perm),statistics = x,
                 call = cl, contrast = contrast)
      class(out) = 'PermuteTest'
      return(out)
    }
    
    if(len_stat == 1 && length(contrasts) == 0){
      return(permute_tail_probs(as.vector(permp), 
                                 observed, contrast = 'identity'))
    }
    contrasts_out = list()
    for(j in seq_len(len_stat)){
      contrasts_out[[j]] = permute_tail_probs(as.vector(permp %*% as.matrix(contrasts[[j]])), 
                                              as.vector(crossprod(observed, contrasts[[j]])), contrast = contrasts[[j]])
    }
    class(contrasts_out) = 'PermuteTestList'
    return(contrasts_out)
}

#' @param x `PermuteTestList`
#' @importFrom generics tidy
#' @describeIn cluster_permute_test return a permutations run using a sequence of contrasts as a `tibble`
#' @export
tidy.PermuteTestList = function(x, ...){
  scalar_vars = c('observed', 'expected', 'p.value', 'mc.se')
  purrr::map_dfr(x, function(y){
    tbl = as_tibble(y['statistics'])
    tbl[scalar_vars] = y[scalar_vars]
    tbl = mutate(tbl, contrast = str_c(y[['contrast']], collapse = ', '))
  })
}

#' @param perm_test \code{PermuteTest} output from \code{cluster_permute_test()}
#' @return A ggplot2 plot
#' @describeIn cluster_permute_test Plot a histogram of permuted vs observed test statistic
#' @export
plot_permute_test = function(perm_test) {
    check_plot_infra()
    if(!inherits(perm_test, 'PermuteTest') && !inherits(perm_test, 'PermuteTestList')) {
        stop('run cluster_permute_test(cdb) first')
    }
  if(inherits(perm_test, 'PermuteTestList')){
    stat_df = tidy.PermuteTestList(perm_test)
  } else{
    stat_df = data.frame(statistics = perm_test$statistics,
                         observed = perm_test$observed)
  }
    plt = ggplot2::ggplot(data = stat_df, ggplot2::aes(x = .data$statistics)) + 
        ggplot2::geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
        ggplot2::geom_vline(xintercept = perm_test$observed, col = 'red') +
        ggplot2::ggtitle("Permuted and observed test statistics")
    if(inherits(perm_test, 'PermuteTestList')){
      plt = plt + ggplot2::facet_wrap(~.data$contrast)
    }
    return(plt)
}
