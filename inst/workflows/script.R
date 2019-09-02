data("ccdb_bcell")

cdb <- filter_cdb(ccdb_bcell, is_cell & high_confidence & productive == 'True' & chain %in% c('IGH', 'IGK', 'IGL'))
Y_kj <- cdb$contig_tbl %>% group_by(chain, cdr3, v_gene, j_gene) %>% summarize(total_umis = sum(umis))
Y_j <- cdb$contig_tbl %>% group_by(chain) %>% summarize(total_umis = sum(umis))
Y_j <- c(rep(Y_j$total_umis[1],sum(Y_kj$chain == 'IGH')), 
             rep(Y_j$total_umis[2],sum(Y_kj$chain=='IGK')), 
                 rep(Y_j$total_umis[3],sum(Y_kj$chain=='IGL')))
Y_kj <- Y_kj$total_umis

Z_kj <- cdb$contig_tbl %>% group_by(pop, barcode) %>% summarize(cells = dplyr::n())
pop <- Z_kj$pop
Z_j <- cdb$contig_tbl %>% group_by(pop) %>% summarize(cells = dplyr::n())
Z_j <- c(rep(Z_j$cells[1],sum(Z_kj$pop == 'balbc')),rep(Z_j$cells[2],sum(Z_kj$pop == 'b6')))
Z_kj <- Z_kj$cells
dat <- data.frame(Z_kj = Z_kj, Z_j = Z_j, pop = pop)
