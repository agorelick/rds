% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rds.R
\name{cluster_probability}
\alias{cluster_probability}
\title{cluster_probability}
\usage{
cluster_probability(k, m)
}
\arguments{
\item{k}{number of remaining cancer samples of the same phylogeny}

\item{m}{number of metastases samples}
}
\value{
probabilility of cluster by chance
}
\description{
Calculate the probability that a tree where all m metastases form a cluster (have one root) exists by chance
}
