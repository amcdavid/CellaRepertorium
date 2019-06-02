#' Generate a legible name for a series of contigs
#'
#' @param contig_frame An `all_contig_annotations.csv` file, output from VDJ Cell ranger.  Importantly, this should contain columns `chain`, `v_gene`, `d_gene`, `j_gene`
#' @param prefix an optional prefix added to each contig, eg, possibly a sample id.
#' @return \code{character}
#' @importFrom dplyr %>%
#' @importFrom stringr str_replace_all
#' @importFrom methods as
#' @importFrom stats as.dist dist hclust na.fail sd
#' @importFrom utils data
#' @export
#'
#' @examples
#' library(dplyr)
#' contig_anno_path = system.file('extdata', 'all_contig_annotations_balbc_1.csv.xz',
#'     package = 'CellaRepertorium')
#' contig_anno = readr::read_csv(contig_anno_path)
#' contig_anno = contig_anno %>% mutate(fancy_name =
#'     fancy_name_contigs(., prefix = 'b6_1'))
#' stopifnot(!any(duplicated(contig_anno$fancy_name)))
fancy_name_contigs = function(contig_tbl, prefix){
    contig_frame = contig_tbl
    notin = setdiff(c('chain', 'v_gene', 'd_gene', 'j_gene'), names(contig_frame))
    if(length(notin>0)) stop("`contig_tbl must contain all of ", paste(notin, collapse = ','))
    chain = contig_frame$chain
paste(contig_frame$v_gene, contig_frame$d_gene, contig_frame$j_gene, sep = ':') %>% str_replace_all('None', '') %>% str_replace_all('IG[KLH]|TR[ABDG]', '') %>% paste(chain, ., sep = ':') %>% paste(prefix, ., sep = '_') %>% make.unique()
}

variable_genes = function(ref = '10X'){
    # This should return genes that lie in the VDJ region for a given annotation
    # And should know something about the chromium reference DB.

}



cleanup_annotations = function(json){
    anno = json[['annotations']]
    atomic = purrr::map(anno, ~ .[setdiff(names(.), 'mismatches')])
    names(atomic) = json[['contig_name']]
    as_data_frame(bind_rows(atomic, .id = 'contig_name'))
}

get_gaps = function(ca){
    vregion = which(ca[['feature.region_type']] == 'L-REGION+V-REGION')
    dregion = which(ca[['feature.region_type']] == 'D-REGION')
    jregion = which(ca[['feature.region_type']] == 'J-REGION')
    vj_gap <- dj_gap <- vd_gap <- NA_integer_
    if(length(vregion) > 0 & length(dregion) > 0) vd_gap = ca[[dregion, 'contig_match_start']] - ca[[vregion,'contig_match_end']]
    if(length(dregion) > 0 & length(jregion) > 0) dj_gap = ca[[jregion, 'contig_match_start']] - ca[[dregion,'contig_match_end']]
    if(length(vregion) > 0 & length(jregion) > 0) vj_gap = ca[[jregion,'contig_match_start']] - ca[[vregion,'contig_match_end']]
    data_frame(vd_gap, dj_gap, vj_gap)
}

read_contig_json = function(file, seq_cols = c('quals', 'aa_sequence', '')){
    jsn = fromJSON(file(anno_file), flatten = TRUE)
    # contig sequences and gaps
}
