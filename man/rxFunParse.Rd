% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rudf.R
\name{rxFunParse}
\alias{rxFunParse}
\alias{rxRmFunParse}
\title{Add user function to rxode2}
\usage{
rxFunParse(name, args, cCode)

rxRmFunParse(name)
}
\arguments{
\item{name}{This gives the name of the user function}

\item{args}{This gives the arguments of the user function}

\item{cCode}{This is the C-code for the new function}
}
\value{
nothing
}
\description{
This adds a user function to rxode2 that can be called.  If needed,
these functions can be differentiated by numerical differences or
by adding the derivatives to rxode2's internal derivative table
with rxode2's \code{rxD} function
}
\author{
Matthew L. Fidler
}
\keyword{internal}
