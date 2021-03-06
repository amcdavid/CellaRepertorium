% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ContigCellDB-methods.R
\name{rbind,ContigCellDB-method}
\alias{rbind,ContigCellDB-method}
\alias{rbind.ContigCellDB}
\title{Combine \code{ContigCellDB} along rows (contigs, cells or clusters).}
\usage{
\S4method{rbind}{ContigCellDB}(..., deparse.level = 1)
}
\arguments{
\item{...}{\code{\link[=ContigCellDB]{ContigCellDB()}}}

\item{deparse.level}{ignored}
}
\value{
\code{\link[=ContigCellDB]{ContigCellDB()}}
}
\description{
The union of the rows in each of the objects is taken,
thus removing any rows that has an exact duplicate.  This
includes all fields, not just the primary key for that table.
The union of the various primary keys is taken.
}
\examples{
data(ccdb_ex)
splat = split_cdb(ccdb_ex, 'chain', 'contig_tbl')
unite = equalize_ccdb(rbind(splat$TRA, splat$TRB), sort = TRUE)
stopifnot(all.equal(unite, ccdb_ex))

}
