#' Calculate distances and perform hierarchical clustering on a set of sequences
#'
#' The distances between AA sequences is defined to be 1-score/max(score) times the median length of the input sequences.
#' The distances between nucleotide sequences is defined to be edit_distance/max(edit_distance) times the median length of input sequences.
#' @param seqs character vector, DNAStringSet or AAStringSet
#' @param type character either `AA` or `DNA` specifying type of `seqs`
#' @param big_memory_brute attempt to cluster more than 4000 sequences?  Clustering is quadratic, so this will take a long time and might exhaust memory
#' @param substitution_matrix a character vector naming a substition matrix used to weight
#' @param cluster_method character passed to `hclust`
#'
#' @seealso hclust, stringDist
#' @return dendrogram of class `hclust`
#' @export
#' @import Biostrings
#' @examples
#' fasta_path = system.file('extdata', 'demo.fasta', package='CellaRepertorium')
#' aaseq = Biostrings::readAAStringSet(fasta_path)[1:100]
#' cls = fine_cluster(aaseq)
#' plot(cls)
fine_cluster = function(seqs, type = 'AA', big_memory_brute = FALSE, substitution_matrix = 'BLOSUM100', cluster = 'hclust', cluster_method = 'complete'){
    if(length(seqs) > 4000 & !big_memory_brute) stop("Won't try to cluster ", length(seqs), " sequences unless `big_memory_brute` = TRUE.  (Make sure you actually have lots of memory!)")
    type = match.arg(type, choices = c('AA', 'DNA'))
    cluster = match.arg(cluster, c('hclust', 'none'))
    if(type == 'AA'){
        biostrings_data_avail = data(package = 'Biostrings')$results
        if(substitution_matrix %in% biostrings_data_avail[,'Item']){
            e = environment()
            data(list = substitution_matrix, envir = e, package = 'Biostrings')
            substitution_matrix = e[[substitution_matrix]]
        }
        if(is.character(seqs)){
            ss = AAStringSet(seqs)
        } else if(inherits(seq, 'AAStringSet')){
            ss = seq
        } else{
            stop("Must be character or AAStringSet")
        }
        sd = as.matrix(stringDist(ss, method = 'substitutionMatrix', substitutionMatrix = substitution_matrix))
        # alignment score along diagonal
        pw = pairwiseAlignment(ss, ss, substitutionMatrix = substitution_matrix, scoreOnly = TRUE)
        diag(sd) = pw

        if(length(ss)>1) sd = as.dist((1 - sd/diag(sd))*median(nchar(ss)))
    } else{
        if(is.character(seqs)){
            ss = DNAStringSet(seqs)
        } else if(inherits(seq, 'DNAStringSet')){
            ss = seq
        } else{
            stop("Must be character or DNAStringSet")
        }
        sd = stringDist(ss, method = 'levenshtein')
        if(length(ss)>1) sd = sd/max(sd)*median(nchar(ss))
    }

    if(cluster == 'hclust'){
        if(length(seqs)>1){
    hc <- stats::hclust(sd, method = cluster_method)
    hc$labels = names(ss)
    return(hc)
        } else{
         return(NA)
        }
    } else if(cluster =='none'){
     return(sd)
    }
}

#' Calculate the entropy of a vector
#'
#' @param v categorical vector
#' @param pseudo_count number of pseudo counts to add on, to stablize empty categories
#' @param na.action how to handle NA values
#'
#' @return the sample entropy
#' @export
#'
#' @examples
#' v2 = gl(2, 4)
#' v4 = gl(4, 4)
#' stopifnot(entropy(v2) < entropy(v4))
#' v_empty = v2[1:4] #empty level 2
#' stopifnot(is.finite(entropy(v_empty)) # pseudo_count
#'
#' np(v4, p = .2, pseudo_count = 0)
#' np(v4, p = .25, pseudo_count = 0)
#' np(v4, p = .25, pseudo_count = .0001)
#'
#' modal_category(v4)
#' modal_category(v4[-1])
entropy = function(v, pseudo_count = length(v)/1000, na.action = na.fail){
    v = na.action(v)
    pr = table(v)/length(v) + pseudo_count
    e = -sum(log(pr)*pr)
    if(e<0) warning("Negative entropy; adjust pseudo_counts?")
    e
}


#' @export
#' @describeIn entropy The number of categories exceeding `p` proportion of the total
np = function(v, p = .05, pseudo_count = p/5, na.action = na.fail){
    v = na.action(v)
    if(length(v) == 0) return(0)
    pr = table(v)/length(v) + pseudo_count
    e = sum(pr>p)
    e
}

#' @export
#' @describeIn entropy The modal category of v.  Ties are broken by lexicographic order of the factor levels.
modal_category = function(v, na.action = na.fail){
    v = na.action(v)
    if(length(v)>0) names(sort(-table(v)))[1] else NA_character_
}
