#' Parse a CD-hit results file as a data.frame
#'
#' Parses a cd-hit (\url{http://weizhong-cluster.ucsd.edu/cdhit_suite/cgi-bin/index.cgi?cmd=cd-hit}) cluster file and returns it as a \code{data.frame}.
#' The \code{data.frame} contains the columns `cluster_idx`: the cluster number, `query_idx`: the fasta identifier of the query sequence, `homology_pct`: the percent similarity to the exemplar sequence, `member_idx`, an arbitrary enumeration of each sequence in the a cluster, `len`: the length of the query sequence.
#'
#' @param file path to a file with cdhit results
#' @return \code{data.frame} with columns
#' @export
#'
#' @examples
#' demo_path = system.file('extdata', 'demo150.clstr.sorted', package='CellaRepertorium')
#' out = read_cdhit(demo_path)
#' @import stringr
read_cdhit = function(file){
    # Dang you CD-HIT, Dang you to heck.
    lns = readLines(file)
    cl_list = list()
    lii = 1
    li = 1
    ci = 1
    ci_name = ''
    while(li <= length(lns)){
        if(substr(lns[li], 1, 1) == '>' || li == length(lns) ){ # new cluster or EOF
            if(li>1){
                cl_list[[ci]] = lns[lii:(li-1)] %>% str_c(collapse = '\n') %>% str_replace_all('>|at|%|aa', '') %>% str_replace_all('(,[:space:])|(\\.{3,3}[:space:]*)', '\t') %>% str_c('\n') %>% str_replace_all('\\*', '100.00')
                names(cl_list)[ci] = ci_name
            }
            lii = li + 1
            ci_name = lns[li]
            ci = ci + 1

        }
        li = li + 1
    }
    tsvs = lapply(cl_list[-1], readr::read_tsv, col_names = c('member_idx', 'len', 'query_idx', 'homology_pct'), col_types = 'iicn')
    all_seqs = dplyr::bind_rows(tsvs, .id = 'cluster_idx') %>% dplyr::mutate(cluster_idx = str_replace_all(cluster_idx, fixed('>Cluster '), '') %>% as.integer)
    all_seqs[-nrow(all_seqs),] # kill last row, which is erroneous, and I can't count.
}

globalVariables('cluster_idx')

CDHIT_SERVER = "http://weizhong-lab.ucsd.edu/cdhit-web-server/cgi-bin/index.cgi?cmd=cd-hit"

#' Query the web interface to cdhit
#'
#' This sends a query to \url{http://weizhong-lab.ucsd.edu/cdhit-web-server/cgi-bin/index.cgi?cmd=cd-hit}
#' and returns a link to results.
#' @param sequences An object of class \code{AAStringSet}
#' @param results_path If non-null, a \code{character} specifying a path to a results.  Directories therein must exist. If null, the R temporary directory will be used.
#' @param identity_cutoff minimum identity to be clustered together
#' @param bandwidth see CDhit docs
#' @param results_timeout number of seconds to wait for result
#' @return path to results file
#' @export
#' @importFrom dplyr %>%
#' @import httr
#' @importFrom utils data download.file
#' @import Biostrings
#' @seealso AAStringSet
#' @examples
#' \dontrun{
#' fasta_path = system.file('extdata', 'demo.fasta', package='CellaRepertorium')
#' aaseq = Biostrings::readAAStringSet(fasta_path)
#' results_path = query_cdhit(aaseq[1:100])
#' read_cdhit(results_path)
#' }
query_cdhit = function(sequences, results_path = NULL, identity_cutoff = .9, bandwidth = 20, results_timeout = 120){
    RELOAD_WAIT = 15
    tmp_path = file.path(tempdir(), str_c(as.numeric(Sys.time())))
    if(!dir.create(tmp_path)) stop("Couldn't create directory within tmpdir()")
    query_file = file.path(tmp_path, 'query.fasta')
    writeXStringSet(sequences, query_file)
    body = list(program = 'cd-hit', SeqF = httr::upload_file(query_file), anno = 1, level = 1, lc1 = identity_cutoff,
                uG = 1, lg = 1, lb = bandwidth, lauL = 0, uAuL = 'unlimited',
                lauS = 0, uAuS = 'unlimited', ls = 0, uS = 'unlimited')
    resp = httr::POST(CDHIT_SERVER, body = body, encode = "multipart")
    check_link = resp$url
    t = 0
    success = FALSE
    jobid = str_extract(check_link, "(?<=JOBID=)[0-9]+")
    message(str_c('Waiting for results to post at ', check_link))
    pb = progress::progress_bar$new(total = results_timeout)
    pb$tick(0)
    while(t < results_timeout){
        pb$tick()
        if( (t %% RELOAD_WAIT) == 0){
            results_try = httr::GET(check_link)
        }
        if(results_try %>% as.character() %>% str_detect('is finished')){   #done!
            cluster_url = str_c('http://weizhong-cluster.ucsd.edu/cdhit_suite/output/', jobid, '/', jobid, '.fas.1.clstr')
            if(is.null(results_path)) results_path = file.path(tmp_path, 'results.clstr')
            if(download.file(cluster_url, results_path)) stop('Error downloading ', cluster_url)
            success = TRUE
            break
        } # try again
        t = t + 1
        Sys.sleep(1)
    }
    if(!success) stop("Timeout trying to get results at ", check_link)
    return(results_path)
}


#' R interface to CDHIT/CDHIT est
#'
#' CDHIT is a greedy algorithm to cluster amino acid or DNA sequences
#' by Fu, Niu, Zhu, Wu and Li (2012).  The R interface is by
#' Thomas Lin Pedersen and was transcribed here because it is not exported from the package FindMyFriends.
#' @param seqs \code{AAseq} or \code{DNAseq} (untested..)
#' @param identity minimum proportion identity
#' @param kmerSize word size.  Set to 5 for 70%-90% identity.  Lower for lesser identity.
#' @param name program name (?)
#' @param showProgress show a status bar
#' @param ... other arguments that can be passed to cdhit, see https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide#CDHIT for details
#' @useDynLib CellaRepertorium
#' @return vector of \code{integer} of length \code{seqs} providing the cluster ID for each sequence
#' @export
#'
#' @examples
#' fasta_path = system.file('extdata', 'demo.fasta', package='CellaRepertorium')
#' aaseq = Biostrings::readAAStringSet(fasta_path)
#' cdhit(aaseq)
cdhit = function(seqs, identity = .9, kmerSize = 5, name = 'CD-Hit', showProgress = interactive(), ...) {
    options = list(...)
    options$i <- tempfile()
    writeXStringSet(seqs, options$i)
    on.exit(unlink(options$i))
    options$n = kmerSize
    options$c = identity
    options = lapply(options, as.character)
    switch(
        class(seqs),
        AAStringSet = cdhitC(options, name, showProgress) + 1,
        DNAStringSet = cdhitestC(options, name, showProgress) + 1,
        stop('seqs must be either AAStringSet or DNAStringSet')
    )
}
