globalVariables('cluster_idx')
#' R interface to CDHIT/CDHITest
#'
#' CDHIT is a greedy algorithm to cluster amino acid or DNA sequences
#' by Fu, Niu, Zhu, Wu and Li (2012).  The R interface is originally by
#' Thomas Lin Pedersen and was transcribed here because it is not exported from the package FindMyFriends, which is orphaned.
#'
#' @param seqs \code{AAseq} or \code{DNAseq} (untested..)
#' @param identity minimum proportion identity
#' @param kmerSize word size.  Set to 5 for 70\%-90\% identity.  Lower for lesser identity.
#' @param min_length Minimum length for sequences to be clustered.  An error if something smaller is passed.
#' @param name program name (?)
#' @param showProgress show a status bar
#' @param only_index if TRUE only return the integer cluster indices, otherwise return a tibble.
#' @param ... other arguments that can be passed to cdhit, see https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHIT for details
#'
#' @useDynLib CellaRepertorium
#' @return vector of \code{integer} of length \code{seqs} providing the cluster ID for each sequence, or a tibble.  See details.
#' @export
#' @importFrom tibble data_frame
#' @importFrom dplyr group_by mutate
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
cdhit = function(seqs, identity = .9, kmerSize = 5, min_length = 6, name = 'CD-Hit', only_index = FALSE, showProgress = interactive(), ...) {
    if(any(width(seqs) < min_length)) stop("Some sequences shorter than `min_length`; remove these or decrease min_length")
    options = list(...)
    options$i <- tempfile()
    writeXStringSet(seqs, options$i)
    on.exit(unlink(options$i))
    options$n = kmerSize
    options$c = identity
    options$l = min_length - 1
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
