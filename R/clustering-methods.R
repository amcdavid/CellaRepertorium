
cluster_germline = function(ccdb, segment_keys = c('v_gene', 'j_gene', 'chain'), cluster_name = 'cluster_idx'){
    contig_tbl = ccdb$contig_tbl
    seg_types = contig_tbl %>% group_by(!!!syms(segment_keys)) %>% summarize() %>% ungroup() %>% mutate(!!sym(cluster_name) := seq_len(nrow(.)))
    cl_con_tbl = left_join_warn(seg_types, contig_tbl, by = segment_keys)
    cluster_tbl = as_tibble(unique(cl_con_tbl[union(cluster_name, segment_keys)]))
    replace_cluster_tbl(ccdb, cluster_tbl, cl_con_tbl, cluster_pk = cluster_name)
}


# Also canonicalize..
#' Perform additional clustering of sequences within groups
#'
#' @param ccdb A ContigCellDB object
#' @param sequence_key `character` naming column in `contig_tbl` with sequence
#' @param type 'AA' or 'DNA'
#' @param max_affinity `numeric` naming the maximal affinity for the sparse affinity matrix that is constructed.  Not currently used.
#' @param keep_clustering_details `logical` -- should output of `fine_cluster_seqs` be kept as a list column
#' @param ... passed to `fine_cluster_seqs`
#' @importFrom dplyr select bind_cols
#'
#' @examples 
#' ccdb_ex = CellaRepertorium:::cdhit_ccdb(ccdb_ex, 'cdr3_nt', type = 'DNA', cluster_name = 'DNA97', identity = .965, min_length = 12, G = 1)
#' ccdb_ex = fine_clustering(ccdb_ex, sequence_key = 'cdr3_nt', type = 'DNA') 
#' @return `Clustering` object with updated `contig_tbl` and `cluster_tbl`
#' @export
fine_clustering = function(ccdb, sequence_key, type, max_affinity = NULL, keep_clustering_details = FALSE, ...){
    contig_tbl = ccdb$contig_tbl
    message('Calculating intradistances on ', nrow(ccdb$cluster_tbl), ' clusters.')
    # run `fine_cluster_seqs` within each cluster_pk
    cluster_tbl = contig_tbl %>% group_by(!!!syms(ccdb$cluster_pk)) %>% summarize(fc = list(fine_cluster_seqs(!!sym(sequence_key), type = type, ...)), n_cluster = n())
    message('Summarizing')
    contig_by_cluster = contig_tbl[union(ccdb$contig_pk, ccdb$cluster_pk)] %>% tidyr::nest(!!!syms(ccdb$contig_pk)) %>%
        right_join(cluster_tbl %>% dplyr::select(!!!syms(ccdb$cluster_pk)), by=ccdb$cluster_pk) # need to make sure these are in the same order!

    if(is.null(max_affinity)){
        max_max = max(purrr::map_dbl(cluster_tbl$fc, 'max_dist'))
    } else {
        max_max = max_affinity
    }
    affinities = purrr::map(cluster_tbl$fc, ~ max_max - .[['distance']])
    aff_matrix = Matrix::.bdiag(affinities)
    d_medoid = purrr::map2_dfr(cluster_tbl$fc, contig_by_cluster$data, function(.x, .y){
        is_medoid = seq_len(nrow(.y)) == .x$medoid
        dplyr::bind_cols(.y, `d(medoid)` = .x$homology, is_medoid = is_medoid)
    })
    contig_tbl = left_join_warn(d_medoid, contig_tbl, by = ccdb$contig_pk, overwrite = TRUE)
    avg_homology = contig_tbl %>% group_by(!!!syms(ccdb$cluster_pk)) %>% summarize(avg_homology = mean(`d(medoid)`))
    cluster_tbl = avg_homology %>% left_join_warn(cluster_tbl, by = ccdb$cluster_pk, overwrite = TRUE)
    if(!keep_clustering_details) cluster_tbl = cluster_tbl %>% select(-fc)
    replace_cluster_tbl(ccdb, cluster_tbl, contig_tbl)

}

right_join_warn = function(...) left_join_warn(..., join = right_join)

left_join_warn = function(x, y, by, overwrite = FALSE, join = left_join, ...){
    if(missing(by)) stop('by must be provided')
    nx = setdiff(names(x), by)
    ny = setdiff(names(y), by)
    if(length(inter <- intersect(nx, ny))>0){
        if(overwrite){
            warning("Overwriting fields ",
                    paste(inter, collapse = ', '),
                    ' in table ', deparse(substitute(y)))
            y = y[setdiff(names(y), inter)]
        } else{
            warning('Tables share fields ',
                    paste(inter, collapse = ', '),
                    ', which will gain a suffix in ',
                    deparse(substitute(y)))
        }
    }
    join(x  = x, y = y, by = by, suffix = c('', '.y'), ...)
}

canonicalize_cluster = function(ccdb, representative, contig_fields = c('cdr3', 'cdr3_nt', 'chain', 'v_gene', 'd_gene', 'j_gene')){
    cluster_tbl = ccdb$contig_tbl %>% filter(is_medoid) %>% select(!!!syms(union(ccdb$cluster_pk, contig_fields)))
    cluster_tbl = cluster_tbl %>% left_join_warn(ccdb$cluster_tbl, by = ccdb$cluster_pk) %>% mutate(representative = make.unique(!!sym(representative)), representative = forcats::fct_reorder(representative, n_cluster))
    contig_tbl = left_join_warn(ccdb$contig_tbl, cluster_tbl %>% select(!!!syms(ccdb$cluster_pk), representative), by = ccdb$cluster_pk)
    replace_cluster_tbl(ccdb, cluster_tbl, contig_tbl)
}

#' Calculate distances and perform hierarchical clustering on a set of sequences
#'
#' The distances between AA sequences is defined to be 1-score/max(score) times the median length of the input sequences.
#' The distances between nucleotide sequences is defined to be edit_distance/max(edit_distance) times the median length of input sequences.
#'
#' @param seqs character vector, DNAStringSet or AAStringSet
#' @param type character either `AA` or `DNA` specifying type of `seqs`
#' @param big_memory_brute attempt to cluster more than 4000 sequences?  Clustering is quadratic, so this will take a long time and might exhaust memory
#' @param method one of 'substitutionMatrix' or 'levenshtein'
#' @param substitution_matrix a character vector naming a substition matrix available in Biostrings, or a substitution matrix itself
#' @param cluster_fun `character`, one of "hclust" or "none", determining if distance matrices should also be clustered with `hclust`
#' @param cluster_method character passed to `hclust`
#'
#' @seealso hclust, stringDist
#' @return `list` containing
#' @export
#' @import Biostrings
#' @examples
#' fasta_path = system.file('extdata', 'demo.fasta', package='CellaRepertorium')
#' aaseq = Biostrings::readAAStringSet(fasta_path)[1:100]
#' cls = fine_cluster_seqs(aaseq)
#' plot(cls$cluster)
fine_cluster_seqs = function(seqs, type = 'AA', big_memory_brute = FALSE, method = 'levenshtein', substitution_matrix = 'BLOSUM100', cluster_fun = 'hclust', cluster_method = 'complete'){
    if(length(seqs) > 4000 & !big_memory_brute) stop("Won't try to cluster ", length(seqs), " sequences unless `big_memory_brute` = TRUE.  (Make sure you actually have lots of memory!)")
    type = match.arg(type, choices = c('AA', 'DNA'))
    cluster_fun = match.arg(cluster_fun, c('hclust', 'none'))
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

    if(length(seqs)>1 && cluster_fun != 'none'){
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
#' @param p proportion threshold
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
