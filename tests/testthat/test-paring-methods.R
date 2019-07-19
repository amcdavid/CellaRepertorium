context("pairing-methods")

pop <- rep(c('a','b'), each = 13)
sample <- c(rep('1',3), rep('2',4), rep('3',6), rep('4', 3), rep('5', 4), rep('6',6))
barcode <- c('AABBCC','CCBBAA','AACCBB', rep(c('ABCABC','CBACBA'), each = 2), rep(c('AAABBB','AAACCC'),each = 3),
             'DDEEFF','FFEEDD','DDFFEE', rep(c('DEFDEF','FEDFED'), each = 2), rep(c('DDDEEE','EEEFFF'),each = 3))
chain <- c(rep('TRA',3),rep('TRA',2), rep('TRB',2), 'TRA','TRA','TRB','TRA','TRB','TRA',
           rep('TRB',3),rep(c('TRA','TRB'),2), 'TRB','TRB','TRA','TRB','TRA','TRB')
contig_id <- c('AABBCC_contig_1','CCBBAA_contig_1','AACCBB_contig_1',
               'ABCABC_contig_1','ABCABC_contig_2','CBACBA_contig_1','CBACBA_contig_2',
               'AAABBB_contig_1','AAABBB_contig_2','AAABBB_contig_3','AAACCC_contig_1','AAACCC_contig_2','AAACCC_contig_3',
               'DDEEFF_contig_1','FFEEDD_contig_1','DDFFEE_contig_1',
               'DEFDEF_contig_1','DEFDEF_contig_2','FEDFED_contig_1','FEDFED_contig_2',
               'DDDEEE_contig_1','DDDEEE_contig_2','DDDEEE_contig_3','EEEFFF_contig_1','EEEFFF_contig_2','EEEFFF_contig_3')

contig_tbl = tibble(barcode = barcode, pop = pop, sample = sample, contig_id = contig_id, chain = chain)
cdb = ContigCellDB(contig_tbl, contig_pk = c('barcode','pop','sample','contig_id'), cell_pk = c('barcode','pop','sample'))
test_that('Chain Pairings',
          {
            tbl <- enumerate_pairing(cdb,chain_recode_fun = 'guess')
            expect_is(tbl, 'tbl')
            expect_equal(nrow(tbl),14)
            expect_equal(sum(tbl$pairing == 'paired'),6)
            expect_equal('pairing' %in% names(enumerate_pairing(cdb)), FALSE)
            expect_equal('canonical' %in% names(enumerate_pairing(cdb)), FALSE)
            expect_error(enumerate_pairing(cdb, chain_key = 'pop'))
            expect_equal(sum(tbl$canonical == 'double-alpha'),3)
            expect_equal(sum(tbl$canonical == 'classical'),8)
          })