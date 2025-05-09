% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conjugate_class.R
\docType{class}
\name{conjugate-class}
\alias{conjugate-class}
\title{Class \code{conjugate} for output from the \code{pcvr::conjugate} function.}
\description{
Comparisons made by the \code{conjugate} function return objects of this class
containing parameters of the prior and posterior distributions, hypothesis tests,
ROPE tests, Bayes Factors, and plots of the posterior.
}
\details{
See \code{methods(class = "conjugate")} for an overview of available methods.
}
\section{Slots}{

\describe{
\item{\code{summary}}{Summary data frame of results}

\item{\code{posterior}}{Posterior distribution as a list of named lists}

\item{\code{prior}}{Prior distribution as a list of named lists}

\item{\code{plot}}{Optionally a plot of the distributions and their differences}

\item{\code{data}}{The data from s1 and s2 arguments to \link{conjugate}.}

\item{\code{call}}{Matched call to \link{conjugate}.}
}}

\seealso{
\code{\link{conjugate}}
}
