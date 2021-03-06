% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMIT.R
\name{tscor}
\alias{tscor}
\title{OTU correlation matrix for each subject}
\usage{
tscor(ana, method = "kendall", subject_var, fill.na = 0)
}
\arguments{
\item{ana}{a phyloseq object of counts/relative abundance data.}

\item{method}{an option of the correlation method ("pearson","kendall","spearman"). The default method is "kendall".}

\item{subject_var}{a numeric vector of subject.}

\item{fill.na}{a number between 0 and 1 to fill the missing value. The default value is 0.}
}
\value{
a list of temporal correlation matrix for each subject.
}
\description{
OTU correlation matrix for each subject
}
