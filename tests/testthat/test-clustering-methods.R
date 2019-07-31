context("clustering-methods")

data("ccdb_ex")
test_that("Cluster contigs by germline properties",{
  cdb <- cluster_germline(ccdb_ex)
  expect_is(cdb,'ContigCellDB')
  expect_gt(ncol(cdb$cluster_tbl),ncol(ccdb_ex$cluster_tbl))
  expect_equal(names(cdb$cluster_tbl),c('cluster_idx','v_gene','j_gene','chain'))
  cdb1 <- cluster_germline(ccdb_ex,segment_keys = c('v_gene', 'j_gene'),cluster_name = 'clusterID')
  expect_gt(ncol(cdb$cluster_tbl),ncol(cdb1$cluster_tbl))
  expect_equal(names(cdb1$cluster_tbl)[1],'clusterID')
  expect_error(cluster_germline(ccdb_ex,segment_keys = c('v_gene', 'm_gene')))
})

ccdb_ex_small <- ccdb_ex
ccdb_ex_small$cell_tbl <- ccdb_ex_small$cell_tbl[1:200,]
ccdb_ex_small <- cdhit_ccdb(ccdb_ex_small,
                           sequence_key = 'cdr3_nt', type = 'DNA', cluster_name = 'DNA97',
                           identity = .965, min_length = 12, G = 1)
ccdb_ex_small_fine <- fine_clustering(ccdb_ex_small, sequence_key = 'cdr3_nt', type = 'DNA')

test_that("additional clustering of sequences within groups",{
  expect_is(ccdb_ex_small_fine,'ContigCellDB')
  expect_gt(ncol(ccdb_ex_small_fine$cluster_tbl),ncol(ccdb_ex_small$cluster_tbl))
  expect_equal(ncol(ccdb_ex_small$contig_tbl)+2,ncol(ccdb_ex_small_fine$contig_tbl))
  ccdb_ex_small_fine1 <- fine_clustering(ccdb_ex_small, sequence_key = 'cdr3_nt', type = 'DNA', keep_clustering_details = TRUE)
  expect_equal(ncol(ccdb_ex_small_fine$cluster_tbl)+1,ncol(ccdb_ex_small_fine1$cluster_tbl))
})

test_that("Find a canonical contig to represent a cluster",{
  ccdb_medoid = canonicalize_cluster(ccdb_ex_small_fine)
  expect_is(ccdb_medoid,'ContigCellDB')
  
})