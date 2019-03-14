#' For each cell, return a single, canonical chain-cluster
#'
#' In single cell data, multiple chains (heavy-light or alpha-beta) are expected.  In some cases, there could be more than two (eg multiple alpha alleles for T cells).
#' This picks a cluster id for each cell based on the overall prevalence of cluster ids over all cells in `tbl`.
#' If order = 1 then the canonical chain-cluster will be the most prevalent, and if order = 2, it will be the 2nd most prevalent, and so on.  Ties are broken arbitrarily (possibly by lexicographic order of `cluster_idx`).
#' @param tbl `data.frame` containing columns specified in `cell_identifiers` and `cluster_idx`
#' @param cell_identifiers `character` vector specifying columns in `tbl` that identify a cell
#' @param cluster_idx `character` specifying the column in `tbl` that identifies a
#' @param order return the 1st, 2nd, 3rd, etc, most common chain-cluster
#'
#'
#' @return `data.frame` with columns from `cell_identifiers` and a single `cluster_idx` for each cell
#' @export
get_canonical_chain = function(tbl, cell_identifiers = 'barcode', cluster_idx = 'cluster_idx', order = 1){
    prevalence = tbl %>% group_by(!!sym(cluster_idx)) %>% summarize(prev = n())
    tbl_prevalence = left_join(tbl, prevalence, by = cluster_idx)
    tbl_order = tbl_prevalence %>% group_by(!!!syms(cell_identifiers)) %>% summarize(!!cluster_idx := dplyr::nth(!!sym(cluster_idx), -order, prev))
    ungroup(tbl_order)
}

# Define canonical chain types per cell

#' Given a family of similar sequences, return a "representative"
#'
#' @param seqs character vector
#' @param medoid_idx optional index into seqs
#' @param warn_if_distinct Should a warning be emitted if there are distinct elements in seqs?
#'
#' If medoid_idx is supplied, the medoid sequence is returned, otherwise the longest
#' sequence is returned
#'
#' @return character vector
#' @export
#'
#' @examples
#' get_canonical_representative(c('apple', 'manzana', 'pomme'))
get_canonical_representative = function(seqs, medoid_idx, warn_if_distinct = FALSE){
    if(!missing(medoid_idx)){
        if(!is.integer(medoid_idx) || medoid_idx < 1 || medoid_idx > length(seqs)) stop("Illegal `medoid_idx`")
        ii = medoid_idx
    } else{
        len = str_length(seqs)
        ii = which.max(len)
    }
    rep = seqs[ii]

    if(warn_if_distinct && length(nomatch <- which(seqs != rep)) > 0) warning("At indices ", nomatch, " sequences did not match the representative")

    return(rep)
}





fine_cluster_by = function(seqs, by, max_dist = NULL, ...){
    seq_list = split(seqs, by)
    fc_list = lapply(seq_list, fine_cluster, ...)
    if(is.null(max_dist)){
        max_max = max(map_dbl(fc_list, 'max_dist'))
    } else {
        max_max = max_dist
    }
    affinities = map(fc_list, ~ max_max - .[['distance']])
    aff_matrix = Matrix::.bdiag(affinities)
    homologies = unlist(map(fc_list, 'homology'), use.names = FALSE)
    medoid = map_dbl(fc_list, 'medoid')
    list(affinities, aff_matrix, homologies, medoid)
}

#' Calculate distances and perform hierarchical clustering on a set of sequences
#'
#' The distances between AA sequences is defined to be 1-score/max(score) times the median length of the input sequences.
#' The distances between nucleotide sequences is defined to be edit_distance/max(edit_distance) times the median length of input sequences.
#' @param seqs character vector, DNAStringSet or AAStringSet
#' @param type character either `AA` or `DNA` specifying type of `seqs`
#' @param big_memory_brute attempt to cluster more than 4000 sequences?  Clustering is quadratic, so this will take a long time and might exhaust memory
#' @param method one of 'substitutionMatrix' or 'levenshtein'
#' @param substitution_matrix a character vector naming a substition matrix available in Biostrings, or a substitution matrix itself
#' @param cluster_method character passed to `hclust`
#' @seealso hclust, stringDist
#' @return dendrogram of class `hclust`
#' @export
#' @import Biostrings
#' @examples
#' fasta_path = system.file('extdata', 'demo.fasta', package='CellaRepertorium')
#' aaseq = Biostrings::readAAStringSet(fasta_path)[1:100]
#' cls = fine_cluster(aaseq)
#' plot(cls$cluster)
fine_cluster = function(seqs, type = 'AA', big_memory_brute = FALSE, method = 'levenshtein', substitution_matrix = 'BLOSUM100', cluster = 'hclust', cluster_method = 'complete'){
    if(length(seqs) > 4000 & !big_memory_brute) stop("Won't try to cluster ", length(seqs), " sequences unless `big_memory_brute` = TRUE.  (Make sure you actually have lots of memory!)")
    type = match.arg(type, choices = c('AA', 'DNA'))
    cluster = match.arg(cluster, c('hclust', 'none'))
    method = match.arg(method, c('substitutionMatrix', 'levenshtein'))

    if(is.character(seqs)){
        if(type == 'AA'){
            ss = AAStringSet(seqs)
            } else {
            ss = DNAStringSet(seqs)
            }
    } else if(inherits(seqs, 'AAStringSet') || inherits(seqs, 'DNAStringSet')){
        ss = seqs
    } else{
        stop("Must be character, AAStringSet or DNAStringSet")
    }

    if(method == 'substitutionMatrix'){
        biostrings_data_avail = data(package = 'Biostrings')$results
        if(substitution_matrix %in% biostrings_data_avail[,'Item']){
            e = environment()
            data(list = substitution_matrix, envir = e, package = 'Biostrings')
            substitution_matrix = e[[substitution_matrix]]
        }
        sd = as.matrix(Biostrings::stringDist(ss, method = method, substitutionMatrix = substitution_matrix))
        # alignment score along diagonal
        pw = Biostrings::pairwiseAlignment(ss, ss, substitutionMatrix = substitution_matrix, scoreOnly = TRUE)
        diag(sd) = pw
        sd  = if(length(ss)>1) (diag(sd) - sd) else matrix(0, nrow = 1, ncol = 1)
    } else{ #levenshtein
        sd = as.matrix(Biostrings::stringDist(ss, method = 'levenshtein'))
        #if(length(ss)>1) sd = sd/max(sd)*median(Biostrings::nchar(ss))
    }

    if(length(seqs)>1 && cluster != 'none'){
        hc = stats::hclust(as.dist(sd), method = cluster_method)
        hc$labels = names(ss)
    } else{
        hc = NULL
    }

    medoid = which.min(colMeans(sd))
    homology = sd[medoid,]
    list(cluster = hc, distance = sd, homology = homology, medoid = medoid, max_dist = max(sd))
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
#' stopifnot(is.finite(entropy(v_empty))) # pseudo_count
#'
#' np(v4, p = .2, pseudo_count = 0)
#' np(v4, p = .25, pseudo_count = 0)
#' np(v4, p = .25, pseudo_count = .0001)
#'
#' modal_category(v4)
#' modal_category(v4[-1])
#'
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
