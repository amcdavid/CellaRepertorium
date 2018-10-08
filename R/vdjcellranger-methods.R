#' Generate a legible name for a series of contigs
#'
#' @param contig_frame An `all_contig_annotations.csv` file, output from VDJ Cell ranger.  Importantly, this should contain columns `chain`, `v_gene`, `d_gene`, `j_gene`
#' @param prefix an optional prefix added to each contig, eg, possibly a sample id.
#' @return \code{character}
#' @export
#'
#' @examples
#' contig_anno_path = system.file('extdata', 'cellranger_contig_annotation.csv', package = 'CellaRepertorium')
#' contig_anno = readr::read_csv(contig_anno_path)
#' contig_anno = contig_anno %>% mutate(fancy_name = fancy_name_contigs(., prefix = paste(sample, pop, sep = '_')))
#' stopifnot(any(duplicated(contig_anno$fancy_name)))
fancy_name_contigs = function(contig_frame, prefix){
    chain = contig_frame$chain
paste(contig_frame$v_gene, contig_frame$d_gene, contig_frame$j_gene, sep = ':') %>% str_replace_all('None', '') %>% str_replace_all('IG[KLH]|TR[ABDG]', '') %>% paste(chain, ., sep = ':') %>% paste(prefix, ., sep = '_') %>% make.unique()

}
