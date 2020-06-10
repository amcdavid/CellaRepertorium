context("ContigDB-methods")

contig_tbl = tibble(contig_key = 1:5, cell_key = c(1, 1, 1, 2, 3))

test_that('Can construct ContigCellDB',
          {
      cdb = ContigCellDB(contig_tbl, 'contig_key', cell_pk = 'cell_key')
      expect_is(cdb, 'ContigCellDB')
      expect_equal(cdb$contig_tbl, contig_tbl)
      expect_equal(nrow(cdb$cell_tbl), 3)
      expect_error(ContigCellDB(contig_tbl, contig_pk = 'cell_key', cell_pk = 'cell_key'))
      expect_equal(access_cdb(cdb,"contig_tbl"), contig_tbl)
      expect_error(replace_cdb(cdb,'cell_pk', 'cell'))
          })


test_that('Can split', {
   data(ccdb_ex)
   splat = split_cdb(ccdb_ex, 'chain', 'contig_tbl', equalize = TRUE)
   expect_equal(length(splat), 2)
   expect_equal(splat$TRA$contig_tbl$chain, rep('TRA', sum(ccdb_ex$contig_tbl$chain == 'TRA')))

   # Drop missing cells
   expect_equal(nrow(splat$TRB$cell_tbl), nrow(unique(dplyr::filter(ccdb_ex$contig_tbl, chain =='TRB')[ccdb_ex$cell_pk])))

   splat_cell = split_cdb(ccdb_ex, c('sample', 'pop'), 'cell_tbl')
})

test_that('Can rbind', {
   data(ccdb_ex)
   splat = split_cdb(ccdb_ex, 'chain', 'contig_tbl')
   unite = equalize_ccdb(rbind(splat$TRA, splat$TRB),sort = TRUE)
   expect_is(unite, 'ContigCellDB')

   expect_equal(unite, ccdb_ex)
   ccdb_ex = cluster_germline(ccdb_ex, segment_keys = 'sample')
   splat_cell = split_cdb(ccdb_ex, c('sample', 'pop'), 'cell_tbl', drop = TRUE, equalize = TRUE)
   unite_cell = do.call(rbind, splat_cell)
   expect_equal(equalize_ccdb(unite_cell, sort = TRUE), equalize_ccdb(ccdb_ex, sort = TRUE))

   splat_cell[[1]] = suppressWarnings(replace_cluster_tbl(splat_cell[[1]], cluster_tbl = tibble(), cluster_pk = character()))
   unite_cluster = do.call(rbind, splat_cell)
   expect_equal(nrow(unite_cluster$cluster_tbl), nrow(ccdb_ex$cluster_tbl) - 1)
})
