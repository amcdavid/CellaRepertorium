# Download original data from 10X website
# PBMCs from C57BL/6 mice - TCR enrichment from amplified cDNA
# V3 chemistry
# https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_mm_c57bl6_pbmc_t
download.file('http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_c57bl6_pbmc_t/vdj_v1_mm_c57bl6_pbmc_t_all_contig_annotations.csv',
              destfile = 'data-raw/vdj_v1_mm_c57bl6_pbmc_t_all_contig_annotations.csv')
download.file('http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_c57bl6_pbmc_t/vdj_v1_mm_c57bl6_pbmc_t_all_contig_annotations.json',
              destfile = 'data-raw/vdj_v1_mm_c57bl6_pbmc_t_all_contig_annotations.json')
# PBMCs from BALB/c mice - TCR enrichment from amplified cDNA
# https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_mm_balbc_pbmc_t
download.file('http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_balbc_pbmc_t/vdj_v1_mm_balbc_pbmc_t_all_contig_annotations.csv',
              destfile = 'data-raw/vdj_v1_mm_balbc_pbmc_t_all_contig_annotations.csv')
download.file('http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_balbc_pbmc_t/vdj_v1_mm_balbc_pbmc_t_all_contig_annotations.json',
              destfile = 'data-raw/vdj_v1_mm_balbc_pbmc_t_all_contig_annotations.json')

# Subsample / split
library(readr)
library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)

csvs = list.files('data-raw', pattern = '.+?t_all_contig.+?\\.csv$', full.names = TRUE)
jsns = list.files('data-raw', pattern = '.+?t_all_contig.+?\\.json$', full.names = TRUE)

# Pull out sample and population names
contig_map = tibble(file = csvs, strain = factor(c('balbc', 'b6'), levels = c('balbc', 'b6')))
jsn_map = tibble(file = jsns, strain = factor(c('balbc', 'b6'), levels = c('balbc', 'b6')))

# read in CSV
contig_map = contig_map %>% rowwise() %>% mutate(contigs = list(read_csv(file)))

# Split into (potentially overlapping) subsets of 200 barcodes
set.seed(1234)
subsample_split = function(contigs, ncells = 200, nsplit = 3){
    cells = contigs %>% group_by(barcode) %>% summarize()
    contigs_splt = map(seq_len(nsplit), function(i){
          cell_samp = cells[sample(nrow(cells), ncells),]
          contigs %>% semi_join(cell_samp, by = 'barcode') %>% ungroup()
    })
}

contig_map = contig_map %>% rowwise() %>% mutate(subsample_list = list(subsample_split(contigs))) %>% unnest(subsample_list) %>% mutate(sample_idx = seq_along(.$strain)) %>% ungroup()

contig_map %>% rowwise() %>% do({
    fname = sprintf('inst/extdata/all_contig_annotations_%s_%d.csv', .$strain, .$sample_idx)
    write_csv(.$subsample_list, path = fname)
    system2('xz', fname)
    tibble()
})

# Derive json
jsn_map = jsn_map %>% rowwise() %>% mutate(json = list(fromJSON(file(file), flatten = FALSE) %>% as_tibble))

# barcodes present in subsamples
ss_bc = contig_map %>% unnest(subsample_list) %>% select(strain, barcode, sample_idx) %>% split(., .$strain)


jsn_ss = map2_dfr(jsn_map$json, ss_bc, function(.x, .y){
    j = right_join(.x, .y, by = 'barcode')
    nest(j, -strain, -sample_idx, .key = 'jsn')
})

jsn_ss %>% rowwise() %>% do({
    code = toJSON(.$jsn)
    fname = sprintf('inst/extdata/all_contig_annotations_%s_%d.json', .$strain, .$sample_idx)
    write_lines(code, path = fname)
    system2('xz', fname)
    tibble()
})

knitr::purl('vignettes/mouse_tcell_qc.Rmd', output = 'data-raw/mouse_tcells_qc.R')
source('data-raw/mouse_tcells_qc.R')
ccdb_ex = ContigCellDB_10XVDJ(contigs_qc, contig_pk = c('pop',   'sample', 'barcode', 'contig_id'), cell_pk = c('pop',   'sample', 'barcode'))
use_data(contigs_qc, overwrite = TRUE)
use_data(ccdb_ex, overwrite = TRUE)
unlink('data-raw/mouse_tcells_qc.R')
