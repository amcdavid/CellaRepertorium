% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering-methods.R
\name{fine_clustering}
\alias{fine_clustering}
\title{Perform additional clustering of sequences within groups}
\usage{
fine_clustering(
  ccdb,
  sequence_key,
  type,
  max_affinity = NULL,
  keep_clustering_details = FALSE,
  ...
)
}
\arguments{
\item{ccdb}{A \code{\link[=ContigCellDB]{ContigCellDB()}} object}

\item{sequence_key}{\code{character} naming column in \code{contig_tbl} with sequence}

\item{type}{'AA' or 'DNA'}

\item{max_affinity}{\code{numeric} naming the maximal affinity for the sparse affinity matrix that is constructed.  Not currently used.}

\item{keep_clustering_details}{\code{logical} -- should output of \code{fine_cluster_seqs} be kept as a list column}

\item{...}{
  Arguments passed on to \code{\link[=fine_cluster_seqs]{fine_cluster_seqs}}
  \describe{
    \item{\code{big_memory_brute}}{attempt to cluster more than 4000 sequences?  Clustering is quadratic, so this will take a long time and might exhaust memory}
    \item{\code{method}}{one of 'substitutionMatrix' or 'levenshtein'}
    \item{\code{substitution_matrix}}{a character vector naming a substitution matrix available in Biostrings, or a  substitution matrix itself}
  }}
}
\value{
\code{\link[=ContigCellDB]{ContigCellDB()}} object with updated \code{contig_tbl} and \code{cluster_tbl}
}
\description{
Perform additional clustering of sequences within groups
}
\examples{
library(dplyr)
data(ccdb_ex)
ccdb_ex_small = ccdb_ex
ccdb_ex_small$cell_tbl = ccdb_ex_small$cell_tbl[1:200,]
ccdb_ex_small = cdhit_ccdb(ccdb_ex_small,
sequence_key = 'cdr3_nt', type = 'DNA', cluster_name = 'DNA97',
identity = .965, min_length = 12, G = 1)
ccdb_ex_small = fine_clustering(ccdb_ex_small, sequence_key = 'cdr3_nt', type = 'DNA')

# Canonicalize with the medoid contig is probably what is most common
ccdb_medoid = canonicalize_cluster(ccdb_ex_small)

# But there are other possibilities.
# To pass multiple "AND" filter arguments must use &
ccdb_umi = canonicalize_cluster(ccdb_ex_small,
contig_filter_args = chain == 'TRA' & length > 500, tie_break_keys = 'umis',
contig_fields = c('chain', 'length'))
ccdb_umi$cluster_tbl \%>\% dplyr::select(chain, length) \%>\% summary()
}
