context('Cluster Logistic Testing')

data(ccdb_ex)
ccdb_ex = cluster_germline(ccdb_ex) %>% fine_clustering('cdr3_nt', 'DNA')
ccdb_ex = canonicalize_cluster(ccdb_ex)
ccdb_ex_trb = filter_cdb(ccdb_ex, chain == 'TRB') %>% canonicalize_cell(, contig_filter_args = chain == 'TRB', contig_fields = c('chain', 'v_gene', 'j_gene'))

theclust = filter(ccdb_ex$cluster_tbl, chain == "TRB", v_gene == "TRBV13-2", j_gene == "TRBJ2-7") %>% mutate( x_ = 1)

manual_cell = left_join(ccdb_ex_trb$cell_tbl, theclust, by = c("chain", "v_gene", "j_gene")) %>% mutate(x_ = ifelse(is.na(x_), 0, 1))

#manual_contig = left_join(ccdb_ex$contig_tbl, theclust)

theform = ~ pop + (1|sample)
test_that('Manual glmer testing equals cluster_logistic_test', {
    if(!requireNamespace('lme4')) skip('Install lme4')
    library(lme4)
    # test against a manually data -- TRUE if the cell was in the cluster, false otherwise
    manual_glmer = glmer(update(theform, x_ ~ .) , data = manual_cell,  family = 'binomial')
    automatic1 = cluster_logistic_test(theform, ccdb = ccdb_ex_trb, cluster_whitelist = theclust, keep_fit = TRUE, silent = TRUE)
    expect_equal(fixef(automatic1$fit[[1]]), fixef(manual_glmer))

    wl2 = tibble(cluster_idx = c(theclust$cluster_idx, ccdb_ex_trb$cluster_tbl$cluster_idx[1:3]))
    automatic2 = cluster_logistic_test(theform, ccdb = ccdb_ex_trb, cluster_whitelist = wl2, keep_fit = TRUE, silent = TRUE)
    expect_equal(unique(automatic2$cluster_idx), wl2$cluster_idx)
    expect_equal(semi_join(automatic2, theclust, by = 'cluster_idx') %>% select(term:p.value), automatic1 %>% select(term:p.value))

    automatic3 = cluster_test_by(ccdb_ex, formula = theform, cluster_whitelist = wl2, silent = TRUE)
    expect_false(any(automatic3$chain == 'TRA'))
    expect_equal(semi_join(automatic3, theclust, by = 'cluster_idx') %>% select(term:p.value), automatic1 %>% select(term:p.value))

    ccdb_ex_trb2 = ccdb_ex_trb
    n_contig = nrow(ccdb_ex_trb$contig_tbl)
    ccdb_ex_trb2$contig_tbl$contig_id = str_c(ccdb_ex_trb2$contig_tbl$contig_id, '_2')
    ccdb_ex_trb2$contig_tbl = ccdb_ex_trb2$contig_tbl[sample(seq_len(n_contig)),]
    ccdb_ex_double = rbind(ccdb_ex_trb, ccdb_ex_trb2)
    expect_equal(nrow(ccdb_ex_double$contig_tbl), 2*nrow(ccdb_ex_trb$contig_tbl))
    expect_equal(ccdb_ex_double$cell_tbl, ccdb_ex_trb$cell_tbl)
    automatic4 = cluster_logistic_test(theform, ccdb = ccdb_ex_double, cluster_whitelist = wl2, silent = TRUE)
    expect_equal(automatic3 %>% select(term:p.value), automatic4 %>% select(term:p.value))
})



# Count cells, not contigs
# Results don't depend on what's in the whitelist
# Results don't depend on order of contigs/cells
#
