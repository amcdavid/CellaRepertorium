globalVariables('cluster_idx')
#' R interface to CDHIT/CDHITest
#'
#' CDHIT is a greedy algorithm to cluster amino acid or DNA sequences based on a minimum identity.
#' By default, in this package it is configured perform ungapped, global alignments with no clipping at start or end.
#' The `identity` is the number of identical characters in alignment
#' divided by the full length of the shorter sequence.
#' Set `s` < 1 to change the minimum coverage of the shorter sequence, which will allow clipping at start or end.
#' Changing `G` = 0 changes the meaning of the `identity` to be the number of
#' identical characters in the alignment divided by the length of the alignment.
#' In this case, you must also set the alignment coverage controls `aL`, `AL`, `aS`, `AS`.
#'
#' CDHit is by Fu, Niu, Zhu, Wu and Li (2012).  The R interface is originally by
#' Thomas Lin Pedersen and was transcribed here because it is not exported from the package FindMyFriends, which is orphaned.
#'
#' @param seqs \code{AAseq} or \code{DNAseq}
#' @param identity minimum proportion identity
#' @param kmerSize word size.  If NULL, it will be chosen automatically based on the identity.
#' You may need to lower it below 5 for AAseq with identity less than .7.
#' @param min_length Minimum length for sequences to be clustered.  An error if something smaller is passed.
#' @param s fraction of shorter sequence covered by alignment.
#' @param name program name (?)
#' @param showProgress show a status bar
#' @param only_index if TRUE only return the integer cluster indices, otherwise return a tibble.
#' @param ... other arguments that can be passed to cdhit, see https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHIT for details.  These will override any default values.
#'
#' @useDynLib CellaRepertorium
#' @return vector of \code{integer} of length \code{seqs} providing the cluster ID for each sequence, or a `tibble`.  See details.
#' @export
#' @importFrom tibble data_frame
#' @importFrom dplyr group_by mutate filter
#'
#' @examples
#' fasta_path = system.file('extdata', 'demo.fasta', package='CellaRepertorium')
#' aaseq = Biostrings::readAAStringSet(fasta_path)
#' # 100% identity, global alignment
#' cdhit(aaseq, identity = 1, only_index = TRUE)[1:10]
#' # 100% identity, local alignment with no padding of endpoints
#' cdhit(aaseq,identity = 1, G = 0, aL = 1, aS = 1,  only_index = TRUE)[1:10]
#' # 100% identity, local alignment with .9 padding of endpoints
#' cdhit(aaseq,identity = 1, G = 0, aL = .9, aS = .9,  only_index = TRUE)[1:10]
#' # a tibble
#' tbl = cdhit(aaseq, identity = 1, G = 0, aL = .9, aS = .9, only_index = FALSE)
cdhit = function(seqs, identity = NULL, kmerSize = NULL, min_length = 6, s = 1, name = 'CD-Hit', only_index = FALSE, showProgress = interactive(), ...) {
    if(any(width(seqs) < min_length)) stop("Some sequences shorter than `min_length`; remove these or decrease min_length")
    uopts = list(...)
    options = list()
    options$i <- tempfile()
    options$s = s
    writeXStringSet(seqs, options$i)
    on.exit(unlink(options$i))
    type = switch(
        class(seqs),
        AAStringSet = 'cdhitC',
        DNAStringSet = 'cdhitestC',
        stop('seqs must be either AAStringSet or DNAStringSet')
    )
    if(type == 'cdhitestC'){ #DNA
        kmerSize = case_when(identity < .8 ~ 4, identity < .85 ~ 5, identity < .88 ~ 6, identity < .9 ~ 7, identity < .95 ~ 9, identity < 1 ~ 10, TRUE ~ 11)
        options = c(options, list(ap = 1, r = 0))
    } else{
        kmerSize = 5
    }
    options$n = kmerSize
    options$c = identity
    options$l = min_length - 1
    options = c(uopts, options)
    options = options[!duplicated(names(options))]
    options = lapply(options, as.character)
    out = switch(
        class(seqs),
        AAStringSet = cdhitC(options, name, showProgress) + 1,
        DNAStringSet = cdhitestC(options, name, showProgress) + 1,
        stop('seqs must be either AAStringSet or DNAStringSet')
    )
    if(only_index) return(out)
    dplyr::data_frame(query_name = names(seqs), seq = as.character(seqs), cluster_idx = out) %>%
        dplyr::group_by(cluster_idx) %>% dplyr::mutate(n_cluster = dplyr::n())
}

##' @describeIn cdhit Run `cdhit` on `ClusterContigDB` object
##' @param object An object of class `ClusterContigDB`
##' @param sequence_key `character` naming the column in the `contig_tbl` containing the sequence to be clustered
##' @param type one of 'DNA' or 'AA'
##' @param cluster_tbl_name What index should the clustering be stored in?  By default, a new, unnamed cluster is added.
##' @export
cdhit_ccdb = function(object, sequence_key, type = c('DNA', 'AA'), cluster_name = 'cluster_idx', ...){
    seqs = object$contig_tbl[[sequence_key]]
    if(length(seqs) < 1) stop("No sequences were provided")
    type = match.arg(type, c('DNA', 'AA'))
    if(type == 'DNA'){
        seqset = DNAStringSet(seqs)
    } else {
        seqset = AAStringSet(seqs)
    }
    cluster_idx = cdhit(seqset, ..., only_index = TRUE)
    contig_tbl = dplyr::mutate(object$contig_tbl, !!sym(cluster_name) := cluster_idx)
    cluster_tbl = as_tibble(unique(contig_tbl[cluster_name]))
    replace_cluster_tbl(object, cluster_tbl, contig_tbl, cluster_pk = cluster_name)
}
