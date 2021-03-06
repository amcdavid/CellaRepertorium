% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{cluster_plot}
\alias{cluster_plot}
\title{Make a plot showing properties of the clustering}
\usage{
cluster_plot(cdb, return_plotlist = FALSE)
}
\arguments{
\item{cdb}{A \code{fine_clustering} \code{ContigCellDB} object}

\item{return_plotlist}{should  a list of \code{ggplot2} plots be returned.  If FALSE, a \code{cowplot} composite is retuned.}
}
\value{
a \code{cowplot} composite or a list of plots.
}
\description{
The number of elements per cluster and the average distance between the medoid and other elements are plotted.
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
cluster_plot(ccdb_ex_small)
}
