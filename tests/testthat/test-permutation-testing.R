context('Permutation Testing')
library(dplyr)
cluster_idx = c(1, 1, 1, 2, 2, 3, 3)
subject = c('A', 'A', 'B', 'B', 'B', 'C', 'C')
contig_tbl = tibble(contig_pk = seq_along(cluster_idx), cluster_idx, subject)
ccdb_test = ContigCellDB(contig_tbl = contig_tbl, contig_pk = 'contig_pk',
                         cell_pk = c('contig_pk', 'subject', 'cluster_idx'), cluster_pk = 'cluster_idx')
test_that('Tests for independence between labels and covariates using permutation of cells',{
  set.seed(123)
  x <- cluster_permute_test(ccdb_test, 'subject', 'cluster_idx', statistic = purity, n_perm  = 50)
  expect_known_value(x,'out/permutation1.rda')
  expect_is(x,'list')
})
