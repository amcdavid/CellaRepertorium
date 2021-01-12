#' Concatenate fasta files and rename sequence headers
#'
#' This facilitates submission to that they can be submitted to highV-quest enmass.
#' @param file_tbl optional `tibble()` mapping between filenames and name
#' @param file_field column in `file_tbl` containing the input filenames
#' @param id_field column in `file_tbl` containing the identifier
#' @param output_name name for output fastas
#' @param sequence_name_prefix ??
#'
#' @importFrom dplyr rowwise transmute do
concat_fasta = function(file_tbl, file_field, id_field, output_name = tempfile(pattern = 'contigs', tmpdir = getwd()), sequence_name_prefix = 'contig_id'){
    if(missing(file_tbl)){
        file_tbl = tibble(file = file_field, id = id_field)
        file_field = 'file'
        id_field = 'id'
    } else if (!inherits('data.frame', file_tbl) || !is.character(file_field) || !(length(file_field) == 1) || !(file_field %in% names(file_field))){
        stop('`file_tbl` must be a data.frame containing a column indexed by `file_field`, a length-1 `character`')
    } else if(!is.character(id_field)){
        id_field = setdiff(names(file_tbl), file_field)
    }
    fasta_tbl = file_tbl %>% rowwise() %>% transmute(fasta = Biostrings::readDNAStringSet(!!sym(file_field)))
    file_decorate = paste0(id_field, '=')
    fasta_tbl = fasta_tbl %>% cbind(file_tbl %>% rowwise() %>% transmute(suffix = paste0(file_decorate, c(!!!syms(id_field)), collapse = ',')))
    fasta_tbl %>% rowwise() %>% do({
        seq_nm = str_c(sequence_name_prefix, '=', names(fasta), ',', suffix)
        names(fasta) = seq_nm
        Biostrings::writeXStringSet(fasta, output_name, append = TRUE)
        tibble()
    })
    invisible(fasta_tbl)
}
