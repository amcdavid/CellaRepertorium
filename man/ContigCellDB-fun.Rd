% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ContigCellDB-methods.R
\name{ContigCellDB}
\alias{ContigCellDB}
\alias{ContigCellDB_10XVDJ}
\title{Construct a ContigCellDB}
\usage{
ContigCellDB(
  contig_tbl,
  contig_pk,
  cell_tbl,
  cell_pk,
  cluster_tbl,
  cluster_pk = character(),
  equalize = TRUE
)

ContigCellDB_10XVDJ(
  contig_tbl,
  contig_pk = c("barcode", "contig_id"),
  cell_pk = "barcode",
  ...
)
}
\arguments{
\item{contig_tbl}{a data frame of contigs, and additional fields describing their properties}

\item{contig_pk}{character vector naming fields in \code{contig_tbl} that uniquely identify a row/contig}

\item{cell_tbl}{a data frame of cell barcodes, and (optional) additional fields describing their properties}

\item{cell_pk}{character vector naming fields in \code{cell_tbl} that uniquely identify a cell barcode}

\item{cluster_tbl}{A data frame that provide cluster assignments for each contig}

\item{cluster_pk}{If \code{cluster_tbl} was provided, a character vector naming fields in \code{cluster_tbl} that uniquely identify a cluster}

\item{equalize}{\code{logical}. Should the contig, cells and clusters be equalized by taking the intersection of their common keys?}

\item{...}{passed to \code{\link[=ContigCellDB]{ContigCellDB()}}}
}
\value{
\code{ContigCellDB}
}
\description{
Construct a ContigCellDB
}
\section{Functions}{
\itemize{
\item \code{ContigCellDB_10XVDJ}: provide defaults that correspond to identifiers in 10X VDJ data
}}

\section{Accessors/mutators}{

See \code{\link[=ContigCellDB-mutate]{$,ContigCellDB-method}} for more on how to access and mutate slots.
See \code{\link[=mutate_cdb]{mutate_cdb()}} and \code{\link[=filter_cdb]{filter_cdb()}} for endomorphic filtering/mutation methods
See \code{\link[=split_cdb]{split_cdb()}} to split into a list, and \code{\link[=rbind.ContigCellDB]{rbind.ContigCellDB()}} for the inverse operation.
}

\examples{
data(contigs_qc)
contigs_qc

cdb = ContigCellDB(contigs_qc, contig_pk = c('barcode', 'pop', 'sample', 'contig_id'),
 cell_pk = c('barcode', 'pop', 'sample'))
 cdb

 # everything that was in contigs_qc
 cdb$contig_tbl

 # Only the cell_pk are included by default (until clustering/canonicalization)
 cdb$cell_tbl

 # Empty, since no cluster_pk was specified
 cdb$cluster_tbl

 # Keys
 cdb$contig_pk
 cdb$cell_pk
 cdb$cluster_pk
}
\seealso{
\code{\link[=ContigCellDB-mutate]{$,ContigCellDB-method}}
}
