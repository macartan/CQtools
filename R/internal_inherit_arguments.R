#' Create argument documentation to inherit
#'
#' @param model A \code{causal_model}. A model object generated by \code{\link[CausalQueries]{make_model}}
#' @param reference_model A causal model as created by \code{\link[CausalQueries]{make_model}}
#' @param analysis_model A causal model as created by \code{\link[CausalQueries]{make_model}}
#' @param parameters A vector of real numbers in [0,1]. Values of parameters to specify (optional). By default, parameters is drawn from \code{model$parameters_df}
#' @param given A conditioning set as a character string that evaluates to a logical, for example 'Y==1'
#' @param sims Integer, number of estimand draws, defaults to max(sims, 4000)
#' @param iter Integer, passed to stan, defaults to 4000
#' @param chains Integer, passed to stan, defaults to 4
#' @param query A character string. An expression defining nodal types to interrogate \code{\link[CausalQueries]{reveal_outcomes}}
#' @param strategy A set of node to be sought
#' @param observed A data.frame with observations.


#'
#' @keywords internal
#'
CQtools_internal_inherit_params <- function(model,
																							 reference_model,
																							 analysis_model,
																							 parameters,
																							 given,
																							 sims,
																							 iter,
																							 chains,
																							 query,
																							 strategy,
																							 observed){

}
