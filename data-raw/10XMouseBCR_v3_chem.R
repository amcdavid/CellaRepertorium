# Download original data from 10X website
# PBMCs from C57BL/6 mice - BCR and 5' from amplified cDNA
# V3 chemistry
# https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_mm_c57bl6_pbmc_t
download.file('http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_c57bl6_pbmc_5gex/vdj_v1_mm_c57bl6_pbmc_5gex_raw_feature_bc_matrix.h5',
              destfile = 'data-raw/vdj_v1_mm_c57bl6_pbmc_5gex_raw_feature_bc_matrix.h5')
download.file('http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_c57bl6_pbmc_b/vdj_v1_mm_c57bl6_pbmc_b_all_contig_annotations.csv',
              destfile = 'data-raw/vdj_v1_mm_c57bl6_pbmc_b_all_contig_annotations.csv')

# PBMCs from BALB/c mice - 5' from amplified cDNA
# V3 chemistry
# https://support.10xgenomics.com/single-cell-vdj/datasets/3.0.0/vdj_v1_mm_balbc_pbmc_5gex
#

download.file('http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_balbc_pbmc_5gex/vdj_v1_mm_balbc_pbmc_5gex_raw_feature_bc_matrix.h5',
              destfile = 'data-raw/vdj_v1_mm_balbc_pbmc_5gex_raw_feature_bc_matrix.h5')

download.file('http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_balbc_pbmc_b/vdj_v1_mm_balbc_pbmc_b_all_contig_annotations.csv',
              destfile = 'data-raw/vdj_v1_mm_balbc_pbmc_b_all_contig_annotations.csv')

# Subsample / split
library(readr)
library(dplyr)
library(purrr)
library(tidyr)

csvs = list.files('data-raw', pattern = '.+?pbmc_b_all.+?\\.csv$', full.names = TRUE)

# Pull out sample and population names
contig_map = tibble(file = csvs, pop = factor(c('balbc', 'b6'), levels = c('balbc', 'b6')))

# read in CSV
all_contigs = contig_map %>% rowwise() %>% mutate(contigs = list(read_csv(file))) %>% unnest()

ccdb_bcell = ContigCellDB_10XVDJ(all_contigs, contig_pk = c('pop',  'barcode', 'contig_id'), cell_pk = c('pop', 'barcode'))
use_data(ccdb_bcell, overwrite = TRUE)
