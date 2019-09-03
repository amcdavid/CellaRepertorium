data("ccdb_bcell")

generate_pseudobulk = function(ccdb, type = c('cell', 'umi')){
  cdb <- filter_cdb(ccdb, is_cell & high_confidence & productive == 'True' & chain %in% c('IGH', 'IGK', 'IGL'))
  if(type == 'umi'){
    bulk <- cdb$contig_tbl %>% group_by(pop, chain, cdr3, v_gene, j_gene) %>% summarize(n_umis = sum(umis))
    bulk <- bulk %>% group_by(pop, chain) %>% mutate(total_umis = sum(n_umis))
  } else if(type == 'cell'){
    bulk <- cdb$contig_tbl %>% group_by(pop, chain, cdr3, v_gene, j_gene) %>% summarize(n_bc = dplyr::n_distinct(barcode))
    bulk <- bulk %>% group_by(pop, chain) %>% mutate(total_bc = sum(n_bc))
  }
  
  return(bulk)
}

bulk_y <- generate_pseudobulk(ccdb_bcell,'umi')
bulk_z <- generate_pseudobulk(ccdb_bcell,'cell')

bulk_y %>% tidyr::complete(.,pop,fill=list(n_umis=0,total_umis=0)) %>% group_by(chain,v_gene,j_gene,cdr3)
