context("ContigDB-methods")

contig_tbl = tibble(contig_key = 1:5, cell_key = c(1, 1, 1, 2, 3))
test_that('Can construct ContigCellDB',
          {
      cdb = ContigCellDB(contig_tbl, 'contig_key', cell_pk = 'cell_key')
      expect_is(cdb, 'ContigCellDB')
      expect_equal(cdb$contig_tbl, contig_tbl)
      expect_equal(nrow(cdb$cell_tbl), 3)
      expect_error(ContigCellDB(contig_tbl, contig_pk = 'cell_key', cell_pk = 'cell_key'))
          })
