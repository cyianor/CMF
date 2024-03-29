% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CMF.R
\name{p_check_sparsity}
\alias{p_check_sparsity}
\title{Internal function for checking whether the input is in the right format}
\usage{
p_check_sparsity(mat, max_row, max_col)
}
\arguments{
\item{mat}{An input matrix of class \code{matrix}}

\item{max_row}{Maximum row index for \code{mat}}

\item{max_col}{Maximum column index for \code{mat}}
}
\value{
\code{TRUE} if the input is in coordinate/triplet format.
\code{FALSE} otherwise.
}
\description{
Internal function for checking whether the input is in the right format
}
