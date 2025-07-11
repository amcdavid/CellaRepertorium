globalVariables('ngrp')

#' Cluster contigs by germline properties
#'
#' @param ccdb [ContigCellDB()]
#' @param segment_keys fields in `contig_tbl` that identify a cluster
#' @param cluster_pk name of cluster to be added to `cluster_tbl`
#'
#' @return [ContigCellDB()]
#' @export
#'
#' @examples
#' data(ccdb_ex)
#' ccdb_ex = cluster_germline(ccdb_ex)
#' ccdb_ex$cluster_tbl
cluster_germline = function(ccdb, segment_keys = c('v_gene', 'j_gene', 'chain'), cluster_pk = 'cluster_idx'){
    contig_tbl = ccdb$contig_tbl
    seg_types = contig_tbl %>% group_by(!!!syms(segment_keys)) %>% summarize() %>% ungroup()
    seg_types[[cluster_pk]] = seq_len(nrow(seg_types))
    cl_con_tbl = left_join_warn(seg_types, contig_tbl, by = segment_keys)
    cluster_tbl = as_tibble(unique(cl_con_tbl[union(cluster_pk, segment_keys)]))
    replace_cluster_tbl(ccdb, cluster_tbl, cl_con_tbl, cluster_pk = cluster_pk)
}

globalVariables(c('fc', 'd(medoid)', 'is_medoid', 'n_cluster'))
# attempt to detect if we need nest_legacy or regular nest
nest_pre1.0 = if('nest_legacy' %in% getNamespaceExports('tidyr')) tidyr::nest_legacy else tidyr::nest

# Also canonicalize..
#' Perform additional clustering of sequences within groups
#'
#' @param ccdb A [ContigCellDB()] object
#' @param sequence_key `character` naming column in `contig_tbl` with sequence
#' @param type 'AA' or 'DNA'
#' @param max_affinity `numeric` naming the maximal affinity for the sparse affinity matrix that is constructed.  Not currently used.
#' @param keep_clustering_details `logical` -- should output of `fine_cluster_seqs` be kept as a list column
#' @inheritDotParams fine_cluster_seqs -seqs -type -cluster_fun -cluster_method
#' @importFrom dplyr select bind_cols
#'
#' @example inst/examples/small_cluster_example.R
#' @return [ContigCellDB()] object with updated `contig_tbl` and `cluster_tbl`
#' @export
fine_clustering = function(ccdb, sequence_key, type, max_affinity = NULL, keep_clustering_details = FALSE, ...){
    contig_tbl = ccdb$contig_tbl
    message('Calculating intradistances on ', nrow(ccdb$cluster_tbl), ' clusters.')
    # run `fine_cluster_seqs` within each cluster_pk
    cluster_tbl = contig_tbl %>% group_by(!!!syms(ccdb$cluster_pk)) %>% summarize(fc = list(fine_cluster_seqs(!!sym(sequence_key), type = type, ...)), n_cluster = dplyr::n())
    cluster_tbl_orig = ccdb$cluster_tbl
    cluster_tbl = left_join_warn(cluster_tbl, cluster_tbl_orig, by = ccdb$cluster_pk, overwrite = TRUE)
    message('Summarizing')
    contig_by_cluster = contig_tbl[union(ccdb$contig_pk, ccdb$cluster_pk)]
    contig_by_cluster = nest_pre1.0(contig_by_cluster, !!!syms(ccdb$contig_pk))
    contig_by_cluster = left_join(cluster_tbl[ccdb$cluster_pk], contig_by_cluster, by=ccdb$cluster_pk)

    if(is.null(max_affinity)){
        max_max = max(purrr::map_dbl(cluster_tbl$fc, 'max_dist'))
    } else {
        max_max = max_affinity
    }
    affinities = purrr::map(cluster_tbl$fc, ~ max_max - .[['distance']])
    aff_matrix = Matrix::.bdiag(affinities)
    d_medoid = purrr::map2_dfr(cluster_tbl$fc, contig_by_cluster$data, function(.x, .y){
        is_medoid = seq_len(nrow(.y)) == .x$medoid
        #if(nrow(.y) != length(.x$distance)) browser()
        dplyr::bind_cols(.y, `d(medoid)` = .x$distance, is_medoid = is_medoid)
    })
    contig_tbl = left_join_warn(d_medoid, contig_tbl, by = ccdb$contig_pk, overwrite = TRUE)
    avg_distance = contig_tbl %>% group_by(!!!syms(ccdb$cluster_pk)) %>% summarize(avg_distance = mean(`d(medoid)`))
    cluster_tbl = avg_distance %>% left_join_warn(cluster_tbl, by = ccdb$cluster_pk, overwrite = TRUE)
    if(!keep_clustering_details) cluster_tbl = cluster_tbl %>% dplyr::select(-fc)
    replace_cluster_tbl(ccdb, cluster_tbl, contig_tbl)

}

#' @describeIn left_join_warn perform a `dplyr::right_join()`
#' @export
right_join_warn = function(...) left_join_warn(..., join = right_join)

#' Perform a `dplyr::left_join()` but check for non-key overlapping fields
#'
#' Perform a `dplyr` join, but either warn if the two tables
#' share non-key fields  If `overwrite = TRUE`, then shared columns will pull from `x` otherwise a suffix will be added to `y`.  To perform this check, `by` must be specified, and it is an error if it is not.
#' @inheritParams dplyr::left_join
#' @param overwrite `logical` -- should non-key fields in y be overwritten using x, or should a suffix (".y") be added
#' @param by `character` specifying columns in `x` and `y` to key on.
#' @param join function giving the type of join to perform, eg, left, right, inner, outer.
#' @param ... passed to joining function
#'
#' @return `data.frame` or `tibble`
#' @export
#'
#' @examples
#' left_join_warn(mtcars, mtcars, by  = 'mpg')
#' left_join_warn(mtcars, mtcars, by = 'mpg', overwrite = TRUE)
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

#' Cluster "And" intersection
#'
#' For each contig present in both `X` and `Y`,
#' a new cluster is defined that combines cluster identities in both `X` and `Y`.
#' In the resulting `ContigCellDB`, two contigs are in the same cluster if they are in the
#' same cluster in `X` and the same cluster
#'  in `Y`.  `X` and `Y` must have matching `contig_pk`.
#'  The `contig_tbl` has fields from `X` for contigs present in both `X` and `Y`.
#'  The  `cell_tbl` from `X` is carried forward unchanged, while the `cluster_tbl`
#'   in the result contains the mapping between the ancestral clustering, and the derived.
#' @param X `ContigCellDB`
#' @param Y `ContigCellDB`
#' @param new_pk optional `character` naming the new pk.
#' @export
#' @examples
#' data(ccdb_ex)
#' ccdb_germ = cluster_germline(ccdb_ex, cluster_pk = 'germline_idx')
#' ccdb_cdr3 = cdhit_ccdb(ccdb_ex, 'cdr3_nt', type = 'DNA',
#' cluster_name = 'DNA97', identity = .965, min_length = 12, G = 1)
#' ccdb_cdr3 = cland(ccdb_cdr3, ccdb_germ)
#'
cland = function(X, Y, new_pk){
    if(!inherits(X, 'ContigCellDB') || !inherits(Y, 'ContigCellDB') || !all(X$contig_pk == Y$contig_pk))
        stop("`X` and `Y` must be `ContigCellDB` with common `contig_pk`.")

    if(all(X$cluster_pk == Y$cluster_pk)){
        oypk = Y$cluster_pk
        nypk = str_c(oypk, '.y')
        Y = replace_cluster_tbl(Y,
                                dplyr::rename(Y$cluster_tbl, !!nypk := !!oypk),
                                dplyr::rename(Y$contig_tbl, !!nypk := !!oypk),
                                nypk, warn = FALSE)
        message('Renaming `Y$cluster_pk`: ', oypk, ' -> ', nypk)
    }
    if(missing(new_pk)) new_pk = str_c(X$cluster_pk, '_and_', Y$cluster_pk)
    contigs = left_join_warn(X$contig_tbl,
                            Y$contig_tbl[c(Y$cluster_pk, Y$contig_pk)], join = dplyr::inner_join, by = X$contig_pk, overwrite = FALSE)
    contigs[[new_pk]] = interaction(contigs[[X$cluster_pk]], contigs[[Y$cluster_pk]])
    # I guess we have mostly been an integer, otherwise ?
    contigs[[new_pk]] = as.integer(factor(contigs[[new_pk]]))
    cluster_tblX = unique(contigs[c(new_pk, X$cluster_pk, Y$cluster_pk)]) %>% left_join_warn(X$cluster_tbl, by = X$cluster_pk, overwrite =  TRUE) %>% left_join_warn(Y$cluster_tbl, by = Y$cluster_pk, overwrite = TRUE)
    cdb = replace_cluster_tbl(X, cluster_tblX, contigs, new_pk, warn = FALSE)
    if (has_fineclustering(cdb)) warning("Rerun `fine_clustering` to update results.")
    cdb
}


has_fineclustering = function(ccdb){
    ('is_medoid' %in% names(ccdb$contig_tbl))
}

#' Find a canonical contig to represent a cluster
#'
#' @param ccdb [ContigCellDB()]
#' @param contig_filter_args an expression passed to [dplyr::filter()].
#' Unlike `filter`, multiple criteria must be `&` together, rather than using
#' commas to separate. These act on `ccdb$contig_tbl`
#' @param tie_break_keys (optional) `character` naming fields in `contig_tbl`
#' that are used sort the contig table in descending order.
#' Used to break ties if `contig_filter_args` does not return a unique contig
#' for each cluster
#' @param order The rank order of the contig, based on `tie_break_keys`
#' to return.  If `tie_break_keys` included an ordered factor (such as chain)
#' this could be used to return the second chain.
#' @param representative an optional field from `contig_tbl` that will be made
#' unique. Serve as a surrogate `cluster_pk`.
#' @param contig_fields Optional fields from `contig_tbl` that will be copied into
#' the `cluster_tbl` from the canonical contig.
#'
#' @inheritParams left_join_warn
#' @return [ContigCellDB()] with some number of clusters/contigs/cells but with "canonical" values copied into `cluster_tbl`
#' @export
#' @seealso [canonicalize_cell()] [left_join_warn()]
#' @example inst/examples/small_cluster_example.R
canonicalize_cluster = function(ccdb, contig_filter_args,
tie_break_keys = character(), order = 1, representative = ccdb$cluster_pk[1], contig_fields = c('cdr3', 'cdr3_nt', 'chain', 'v_gene', 'd_gene', 'j_gene'), overwrite = TRUE){
    if(missing(contig_filter_args) && 'is_medoid' %in% names(ccdb$contig_tbl)){
        message("Filtering `contig_tbl` by `is_medoid`, override by setting `contig_filter_args == TRUE`")
        contig_filter_args = quote(is_medoid)
    }
    if(!has_fineclustering(ccdb)){
        stop('Run `fine_clustering(ccdb)` first.')
    }
    req_contig_fields = unique(c(contig_fields, representative, tie_break_keys))
    if(length(missing_contig <- setdiff(req_contig_fields, names(ccdb$contig_tbl))) > 0) stop('`contig_tbl` is missing fields, ', paste(missing_contig, collapse = ', '), '.')
    sub_contig_tbl = filter(.data = ccdb$contig_tbl, !!rlang::enexpr(contig_filter_args))

    if(nrow(sub_contig_tbl) != nrow(ccdb$cluster_tbl)){
        message('Subset of `contig_tbl` has ', nrow(sub_contig_tbl) , ' rows for ',
                nrow(ccdb$cluster_tbl), ' clusters. Filling missing values and breaking ties ', appendLF = FALSE)
        if(length(tie_break_keys) == 0) message('arbitrarily by first contig.')
        else message('with ', paste(tie_break_keys, sep = ', '), '.')
    }
    # setup quosures to arrange the data
    arranging = purrr::map(tie_break_keys, ~ rlang::quo(desc(!!sym(.x))))

    # take first row of each cluster
    cluster_tbl = sub_contig_tbl %>% group_by(!!!syms(ccdb$cluster_pk)) %>% dplyr::arrange(!!!arranging)
    idx = cluster_tbl %>% dplyr::transmute(ngrp = dplyr::n(), idx = seq_along(ngrp))
    cluster_tbl = cluster_tbl[idx$idx==order,,drop = FALSE]
    cluster_tbl = cluster_tbl %>% dplyr::select(!!!syms(unique(c(ccdb$cluster_pk, contig_fields, representative))))

    # fill any missing clusters after the filtering
    cluster_tbl = cluster_tbl  %>% ungroup() %>% right_join_warn(ccdb$cluster_tbl, by = ccdb$cluster_pk, overwrite = overwrite)
    cluster_tbl = cluster_tbl %>% mutate(representative = make.unique(as.character(!!sym(representative))), representative = forcats::fct_reorder(representative, n_cluster))
    contig_tbl = left_join_warn(ccdb$contig_tbl, cluster_tbl %>% dplyr::select(!!!syms(ccdb$cluster_pk), representative), by = ccdb$cluster_pk)
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
#' @param substitution_matrix a character vector naming a substitution matrix available in Biostrings, or a  substitution matrix itself
#' @param cluster_fun `character`, one of "hclust" or "none", determining if distance matrices should also be clustered with `hclust`
#' @param cluster_method character passed to `hclust`
#'
#' @seealso [hclust()], [pwalign::stringDist()]
#' @return `list`
#' @export
#' @import Biostrings
#' @import pwalign
#' @examples
#' fasta_path = system.file('extdata', 'demo.fasta', package='CellaRepertorium')
#' aaseq = Biostrings::readAAStringSet(fasta_path)[1:100]
#' cls = fine_cluster_seqs(aaseq, cluster_fun = 'hclust')
#' plot(cls$cluster)
fine_cluster_seqs = function(seqs, type = 'AA', big_memory_brute = FALSE, method = 'levenshtein', substitution_matrix = 'BLOSUM100', cluster_fun = 'none', cluster_method = 'complete'){
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
        sd = as.matrix(pwalign::stringDist(ss, method = method, substitutionMatrix = substitution_matrix))
        # alignment score along diagonal
        pw = pwalign::pairwiseAlignment(ss, ss, substitutionMatrix = substitution_matrix, scoreOnly = TRUE)
        diag(sd) = pw
        sd  = if(length(ss)>1) (diag(sd) - sd) else matrix(0, nrow = 1, ncol = 1)
    } else{ #levenshtein
        sd = as.matrix(pwalign::stringDist(ss, method = 'levenshtein'))
        #if(length(ss)>1) sd = sd/max(sd)*median(Biostrings::nchar(ss))
    }

    if(length(seqs)>1 && cluster_fun != 'none'){
        hc = stats::hclust(as.dist(sd), method = cluster_method)
        hc$labels = names(ss)
    } else{
        hc = NULL
    }

    medoid = which.min(colMeans(sd))
    distance = sd[medoid,]
    list(cluster = hc, distance_mat = sd, distance = distance, medoid = medoid, max_dist = max(sd))
}

#' Calculate the entropy of a vector
#'
#' @param v categorical vector
#' @param pseudo_count number of pseudo counts to add on, to stabilize empty categories
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
