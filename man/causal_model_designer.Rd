% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CQ_designer.R
\name{causal_model_designer}
\alias{causal_model_designer}
\title{Integrate CausalQueries with DeclareDesign}
\usage{
causal_model_designer(
  reference_model = make_model("X -> Y"),
  param_type = "prior_draw",
  prior_data = NULL,
  analysis_model = reference_model,
  n = 1,
  query,
  ...
)
}
\arguments{
\item{reference_model}{A causal model as created by \code{\link[CausalQueries]{make_model}}}

\item{param_type}{parameter types to use for data draw; defaults to prior}

\item{prior_data}{if prior data available this gets added to new data}

\item{analysis_model}{model used for analysis. updating is based on priors from this model so any prior data based updating is communicate by re-providing prior data.}

\item{n}{Number of observations}

\item{query}{A character string. An expression defining nodal types to interrogate \code{\link[CausalQueries]{reveal_outcomes}}}

\item{...}{arguments passed to \code{query_model}}
}
\description{
Define a model, estimand, and answer strategy. will need @ import DeclareDesign
}
\examples{
require("DeclareDesign")

my_design <- causal_model_designer(
  reference_model = make_model("X -> Y") \%>\%
    set_priors(c(1,1, 1, 1, 7, 1)),
  analysis_model = make_model("X -> Y"),  # Prior model
  n = 100,
  query = list(ATE = "Y[X=1] - Y[X=0]"),
  given = "X==1 & Y==1")

df <- DeclareDesign::draw_data(my_design)
estimand <- draw_estimands(my_design)

# Estimation and diagnosis
# options(mc.cores = parallel::detectCores())
estimate <- draw_estimates(my_design)
diag <- diagnose_design(my_design, sims = 2)
diag
sim <- simulate_design(my_design, sims = 2)

}
