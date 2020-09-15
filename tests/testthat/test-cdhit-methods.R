context("cdhit-methods")
library(Biostrings)
seqs = AAStringSet(c('CASSPGRGAYEQYF', 'CASSPGRGAYEQYF', 'CASSSGRGAYEQYF', 'CASSPGRGAYEQY', 'CASSPF', 'CASSP'))
test_that("cdhit returns appropriate lengths", {
    expect_equal(length(cdhit(seqs, min_length = 5, kmerSize = 3, only_index = TRUE)), length(seqs))
    expect_error(cdhit(seqs, min_length = 7, kmerSize = 3))
})

test_that("cdhit only clusters identical sequences with identity = 1, G= 0, aL = 1, aS = 1", {
    # Total equality
    res = cdhit(seqs = seqs, kmerSize = 3, min_length = 5, only_index = TRUE, identity = 1, G = 0, aL = 1, aS = 1)
    expect_equal(res, c(1, 1, 2, 3, 4, 5))
    # Equality except for 90% differences in ends
    res = cdhit(seqs = seqs, kmerSize = 3, min_length = 5, only_index = TRUE, identity = 1, G = 0, aL = .9, aS = .9, s = 0 )
    expect_equal(res, c(1, 1, 2, 1, 3, 4))
    # 90% equality, same length
    res = cdhit(seqs = seqs, kmerSize = 3, min_length = 5, only_index = TRUE, identity = .9, G = 0, aL = 1, aS = 1)
    expect_equal(res, c(1, 1, 1, 2, 3, 4))
})

test_that("Don't need to be sorted", {
    # if this gets too wild, the results are equivalent but re-numbered
    scram = c(2, 1, 3, 4, 5, 6)
    seqs_scramble = seqs[scram]
    res = cdhit(seqs, identity = 1, kmerSize = 3, min_length = 5, only_index = TRUE)
    res_scram = cdhit(seqs_scramble, identity = 1, kmerSize = 3, min_length = 5, only_index = TRUE)
    expect_equal(res_scram, res[scram])
})

data(ccdb_ex)

test_that("cdhit_ccdb do update", {
  res = cdhit_ccdb(ccdb_ex, 'cdr3_nt', type = 'DNA', cluster_name = 'DNA97', identity = .965, min_length = 12, G = 1)
  expect_is(res, 'ContigCellDB')
  expect_equal(res$cluster_pk, 'DNA97')
  expect_equal(ncol(ccdb_ex$contig_tbl)+1,ncol(res$contig_tbl))
})
