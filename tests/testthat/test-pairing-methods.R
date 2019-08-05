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
    expect_equal(canon2$cell_tbl[ccdb_ex$cell_pk], ccdb_ex$cell_tbl[ccdb_ex$cell_pk])
    expect_true(all(is.na(canon2$cell_tbl$umis)))
})

test_that('Overwrites fields appropriately', {
    expect_warning(canon3 <- canonicalize_cell(canon1, contig_filter_args = FALSE, contig_fields = 'umis', overwrite = TRUE),  'Overwriting')
    expect_true(all(is.na(canon3$cell_tbl$umis)))
    expect_warning(canon4 <- canonicalize_cell(canon1, contig_filter_args = FALSE, contig_fields = 'umis', overwrite = FALSE), 'suffix')
    expect_equal(canon4$cell_tbl$umis.y, canon1$cell_tbl$umis)
})

