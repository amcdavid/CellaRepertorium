context('ccdb_join')
library(SingleCellExperiment)
contig_tbl = tibble(contig_key = 1:5, cell_key = c(1, 1, 1, 2, 3))
cdb = ContigCellDB(contig_tbl, 'contig_key', cell_pk = 'cell_key')

coldata = DataFrame(cell_key = c(1, 2, 3), feature = c(4, 9, 1))
sce = SingleCellExperiment::SingleCellExperiment(colData = coldata)

test_that('sce join',
          {
              expect_s4_class(ccdb_join(sce, cdb), 'ContigCellDB')
              expect_equal(colnames(ccdb_join(sce, cdb)), c('cell_key', 'feature'))
              expect_equal(ccdb_join(sce, cdb)$feature, c(4, 9, 1))
          })

test_that('data.frame join',
          {
              expect_s4_class(ccdb_join(coldata, cdb), 'ContigCellDB')
              expect_equal(colnames(ccdb_join(coldata, cdb)), c('cell_key', 'feature'))
              expect_equal(ccdb_join(coldata, cdb)$feature, c(4, 9, 1))
          })
