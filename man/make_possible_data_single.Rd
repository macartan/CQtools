% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_possible_data.R
\name{make_possible_data_single}
\alias{make_possible_data_single}
\title{Make possible data for a single strategy step}
\usage{
make_possible_data_single(
  model,
  observed = NULL,
  N = 1,
  withins = FALSE,
  conditions = TRUE,
  vars = NULL
)
}
\arguments{
\item{model}{A \code{causal_model}. A model object generated by \code{\link[CausalQueries]{make_model}}}

\item{observed}{A data.frame with observations.}

\item{N}{Number of node to seek}

\item{withins}{logical Whether to seek node within existing data}

\item{conditions}{A list of character strings indicating for which cases data should be gathered. Options are: (i) to gather additional data on node specified via \code{vars} for any possible cases in the model ("any"), (ii) to gather data in all cases within an observed dataset ("within"), or (iii) to specify the subset of cases for which within-case data should be collected (e.g. "Y == 1").}

\item{vars}{Variables to be sought or NA. If NA \code{make_possible_data} gathers data on all node containing NA for the specified data strategy.}
}
\value{
A dataset
}
\description{
Creates a database of possible data from a data strategy.
Users can gather additional data on node specified via \code{vars} for any possible cases in the model ("any"). Or they can
gather data in all cases within an observed dataset ("withins"). Or they can specify  the subset of cases for which withins-case data should be collected (e.g. "Y == 1").
}
\examples{
model <- make_model("X->M->Y")  \%>\%
   set_restrictions(c("Y[M=1]<Y[M=0]", "M[X=1]<M[X=0]")) \%>\%
   set_parameter_matrix()
df <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1))
observed <- collapse_data(df, model, drop_family = TRUE)

make_possible_data_single(model, observed = observed, vars = "M", withins = TRUE, N = 2)
}
\keyword{internal}
