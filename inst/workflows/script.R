data("ccdb_bcell")
ccdb_bcell =  filter_cdb(ccdb_bcell, is_cell & high_confidence & productive == 'True' & chain %in% c('IGH', 'IGK', 'IGL'))

by_chain = split_cdb(ccdb_bcell, 'chain', equalize = TRUE)

bulk_y <- generate_pseudobulk(by_chain$IGH,type = 'umi', class_keys = c('cdr3', 'v_gene', 'j_gene'), total_keys = 'pop')

bulk_z <- generate_pseudobulk(by_chain$IGH,type = 'cell', class_keys = c('cdr3', 'v_gene', 'j_gene'), total_keys = 'pop')

join = left_join(bulk_y, bulk_z, by = c('pop', 'cdr3', 'v_gene', 'j_gene'), suffix = c('_y', '_z'))

ggplot(data = join, aes(y=n_class_y/total_y, x = n_class_z/total_z, color = pop)) + geom_point() + geom_smooth(method = 'lm')

glms_y = bulk_y %>% group_by(v_gene, j_gene, cdr3) %>% do({
  fit = glm(cbind(n_class, total-n_class)~pop, data =., family = 'binomial')
  tidy(fit)
})

glms_z = bulk_z %>% group_by(v_gene, j_gene, cdr3) %>% do({
  fit = glm(cbind(n_class, total-n_class)~pop, data =., family = 'binomial')
  tidy(fit)
})

is_sig_y = filter(glms_y, term == 'popb6') %>% filter(p.value < .05)
is_sig_z =  filter(glms_z, term == 'popb6') %>% filter(p.value < .05)
