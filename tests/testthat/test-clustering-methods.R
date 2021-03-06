context("clustering-methods")

data("ccdb_ex")
test_that("Cluster contigs by germline properties",{
  cdb <- cluster_germline(ccdb_ex)
  expect_is(cdb,'ContigCellDB')
  expect_gt(ncol(cdb$cluster_tbl),ncol(ccdb_ex$cluster_tbl))
  expect_equal(names(cdb$cluster_tbl),c('cluster_idx','v_gene','j_gene','chain'))
  cdb1 <- cluster_germline(ccdb_ex,segment_keys = c('v_gene', 'j_gene'),cluster_pk = 'clusterID')
  expect_gt(ncol(cdb$cluster_tbl),ncol(cdb1$cluster_tbl))
  expect_equal(names(cdb1$cluster_tbl)[1],'clusterID')
  expect_error(cluster_germline(ccdb_ex,segment_keys = c('v_gene', 'm_gene')))
})


ccdb_ex_small <- ccdb_ex
ccdb_ex_small$cell_tbl <- ccdb_ex_small$cell_tbl[1:200,]
ccdb_ex_small <- cdhit_ccdb(ccdb_ex_small,
                           sequence_key = 'cdr3_nt', type = 'DNA', cluster_pk = 'DNA97',
                           identity = .965, min_length = 12, G = 1)
ccdb_ex_small_fine <- fine_clustering(ccdb_ex_small, sequence_key = 'cdr3_nt', type = 'DNA')

contig_tbl <- tibble(contig_key = 1:12, cell_key = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
                    cluster_idx = rep(1:3, each = 4),
                    seq = c('AACCAA','AACCAA','AAAAAA','AAAAAA',
                            'AACCTTGG','ACTG','AACCTT','AACCGG',
                            'ACTGA','ACTG','AACTGA','AACTGA'))
cdb = ContigCellDB(contig_tbl, 'contig_key', cell_pk = 'cell_key', cluster_pk = 'cluster_idx')

test_that("additional clustering of sequences within groups",{
  expect_is(ccdb_ex_small_fine,'ContigCellDB')
  expect_gt(ncol(ccdb_ex_small_fine$cluster_tbl),ncol(ccdb_ex_small$cluster_tbl))
  expect_equal(ncol(ccdb_ex_small$contig_tbl)+2,ncol(ccdb_ex_small_fine$contig_tbl))
  ccdb_ex_small_fine1 <- fine_clustering(ccdb_ex_small, sequence_key = 'cdr3_nt', type = 'DNA', keep_clustering_details = TRUE)
  expect_equal(ncol(ccdb_ex_small_fine$cluster_tbl)+1,ncol(ccdb_ex_small_fine1$cluster_tbl))
})

cdb_fine <- fine_clustering(cdb,sequence_key = 'seq',type = 'DNA')
test_that("whether fine_cluster's calculation is correct",{
  expect_equal(dplyr::filter(cdb_fine$contig_tbl, is_medoid) %>% dplyr::pull(contig_key), c(1, 7, 9))
  avg_distance <- cdb_fine$cluster_tbl[['avg_distance']]
  expect_equal(avg_distance[1], 1)
  expect_equal(avg_distance[2], 1.75)
  expect_equal(avg_distance[3], 0.75)
})

test_that("Find a canonical contig to represent a cluster",{
  ccdb_medoid <- canonicalize_cluster(ccdb_ex_small_fine)
  expect_is(ccdb_medoid,'ContigCellDB')
  expect_error(cdb_canonicalize <- canonicalize_cluster(cdb_fine), 'missing fields')
  cdb_canonicalize <- canonicalize_cluster(cdb_fine,contig_fields = 'seq', representative = 'seq')
  expect_is(cdb_canonicalize$cluster_tbl$representative, 'factor')
 expect_equal(as.character(cdb_canonicalize$cluster_tbl$representative), dplyr::filter(cdb_fine$contig_tbl, is_medoid) %>% dplyr::pull(seq))
})

test_that("cland preserves contigs / cells", {
  cdb1 = cluster_germline(ccdb_ex, 'v_gene', cluster_pk = 'v_idx')
  cdb2 = cluster_germline(ccdb_ex, segment_keys = c('j_gene'), cluster_pk = 'j_idx')
  cdb3 = cland(cdb1, cdb2, new_pk = 'new')
  inames = intersect(names(cdb1$contig_tbl), names(cdb3$contig_tbl))
  expect_mapequal(cdb1$contig_tbl[inames], cdb3$contig_tbl[inames])
  #expect_mapequal(cdb1$contig_tbl, cdb3$contig_tbl)

  expect_equal(cdb1$cell_tbl, cdb3$cell_tbl)
  cdb4 = cluster_germline(ccdb_ex, segment_keys = c('v_gene', 'j_gene'), cluster_pk = 'cluster_idx3')
  j = hushWarning(left_join_warn(cdb4$contig_tbl, cdb3$contig_tbl, by = cdb3$contig_pk, overwrite = TRUE), 'Overwriting')
  ncomb = j %>% dplyr::group_by(new) %>% dplyr::summarize(ncomb = dplyr::n_distinct(cluster_idx3))
  expect_setequal(ncomb$ncomb, 1)
})

test_that("cland renames and is self-consistent", {
  cdb1 = cluster_germline(ccdb_ex, 'v_gene')
   hushWarning(expect_message(cdb2<- cland(cdb1, cdb1, new_pk = 'new'), 'Renaming'), 'Overwriting')
   expect_equal(cdb1$cluster_tbl$cluster_idx, cdb2$cluster_tbl$cluster_idx)
})
