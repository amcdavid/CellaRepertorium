context("pairing-methods")

barcode <- as.character(c(1:3,
             rep(4:5, each = 2),
             rep(6:7,each = 3),
             8:10, rep(11:12, each = 2), rep(13:14,each = 3)))
chain <- c(rep('TRA',3),
           rep('TRA',2), rep('TRB',2), 'TRA','TRA','TRB','TRA','TRB','TRA',
           rep('TRB',3),rep(c('TRA','TRB'),2), 'TRB','TRB','TRA','TRB','TRA','TRB')
contig_id <- make.unique(barcode)

contig_tbl = tibble(barcode = barcode, contig_id = contig_id, chain = chain)
cdb = ContigCellDB(contig_tbl, contig_pk = c('barcode','contig_id'), cell_pk = 'barcode')
test_that('Chain Pairings',
          {
            tbl <- enumerate_pairing(cdb,chain_recode_fun = 'guess')
            expect_is(tbl, 'data.frame')
            expect_equal(nrow(tbl),length(unique(barcode)))
            expect_equal(sum(tbl$pairing == 'paired'),6)
            expect_equal('pairing' %in% names(enumerate_pairing(cdb)), FALSE)
            expect_equal('canonical' %in% names(enumerate_pairing(cdb)), FALSE)
            expect_error(enumerate_pairing(cdb, chain_key = 'pop'))
            expect_equal(sum(tbl$canonical == 'double-alpha'),3)
            expect_equal(sum(tbl$canonical == 'classical'),8)
          })

context('Canonicalize cells')
data("ccdb_ex")
canon1 = canonicalize_cell(ccdb_ex, contig_filter_args = TRUE, contig_fields = 'umis')
test_that('Canonicalization preserves number and order of cells', {
    canon2 = canonicalize_cell(ccdb_ex, contig_filter_args = FALSE, contig_fields = 'umis')
    expect_equivalent(canon2$cell_tbl[ccdb_ex$cell_pk], ccdb_ex$cell_tbl[ccdb_ex$cell_pk])
    expect_true(all(is.na(canon2$cell_tbl$umis)))
})

test_that('Overwrites fields appropriately', {
    expect_warning(canon3 <- canonicalize_cell(canon1, contig_filter_args = FALSE, contig_fields = 'umis', overwrite = TRUE),  'Overwriting')
    expect_true(all(is.na(canon3$cell_tbl$umis)))
    expect_warning(canon4 <- canonicalize_cell(canon1, contig_filter_args = FALSE, contig_fields = 'umis', overwrite = FALSE), 'suffix')
    expect_equal(canon4$cell_tbl$umis.y, canon1$cell_tbl$umis)
})


context('Pairing Tables')
library(dplyr)
tbl <- tibble(clust_idx = gl(3, 2), cell_idx = rep(1:3, times = 2), contig_idx = 1:6)
ccdb <- ContigCellDB(tbl, contig_pk = c('cell_idx', 'contig_idx'), cell_pk = 'cell_idx', cluster_pk = 'clust_idx')

tbl2 <- bind_rows(tbl, tbl %>% mutate(cell_idx = rep(4:6, times = 2)))
ccdb2 <- ContigCellDB(tbl2, contig_pk = c('cell_idx', 'contig_idx'), cell_pk = 'cell_idx', cluster_pk = 'clust_idx')



test_that('Generate a list of tables representing clusters paired in cells',{
  pt2 <- pairing_tables(ccdb, canonicalize_by_prevalence, min_expansion = 1)
  expect_known_value(pt2$cell_tbl, "out/pairing1.rda")

  pt3 <- pairing_tables(ccdb2, canonicalize_by_prevalence, min_expansion = 1)
  expect_known_value(pt3$cell_tbl,  "out/pairing2.rda")

  ccdb2$contig_tbl = ccdb2$contig_tbl %>% mutate(umis = 1, reads = 1, chain = rep(c('TRA', 'TRB'), times = 6))
  pt4 <- pairing_tables(ccdb2, canonicalize_by_chain, min_expansion = 1, table_order = 2)
  expect_known_value(pt4$cell_tbl, "out/pairing3.rda")

})

