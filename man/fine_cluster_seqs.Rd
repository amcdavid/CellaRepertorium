% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering-methods.R
\name{fine_cluster_seqs}
\alias{fine_cluster_seqs}
\title{Calculate distances and perform hierarchical clustering on a set of sequences}
\usage{
fine_cluster_seqs(
  seqs,
  type = "AA",
  big_memory_brute = FALSE,
  method = "levenshtein",
  substitution_matrix = "BLOSUM100",
  cluster_fun = "none",
  cluster_method = "complete"
)
}
\arguments{
\item{seqs}{character vector, DNAStringSet or AAStringSet}

\item{type}{character either \code{AA} or \code{DNA} specifying type of \code{seqs}}

\item{big_memory_brute}{attempt to cluster more than 4000 sequences?  Clustering is quadratic, so this will take a long time and might exhaust memory}

\item{method}{one of 'substitutionMatrix' or 'levenshtein'}

\item{substitution_matrix}{a character vector naming a substitution matrix available in Biostrings, or a  substitution matrix itself}

\item{cluster_fun}{\code{character}, one of "hclust" or "none", determining if distance matrices should also be clustered with \code{hclust}}

\item{cluster_method}{character passed to \code{hclust}}
}
\value{
\code{list}
}
\description{
The distances between AA sequences is defined to be 1-score/max(score) times the median length of the input sequences.
The distances between nucleotide sequences is defined to be edit_distance/max(edit_distance) times the median length of input sequences.
}
\examples{
fasta_path = system.file('extdata', 'demo.fasta', package='CellaRepertorium')
aaseq = Biostrings::readAAStringSet(fasta_path)[1:100]
cls = fine_cluster_seqs(aaseq, cluster_fun = 'hclust')
plot(cls$cluster)
}
\seealso{
\code{\link[=hclust]{hclust()}}, \code{\link[Biostrings:stringDist]{Biostrings::stringDist()}}
}
