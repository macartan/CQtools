#' Integrate CausalQueries with DeclareDesign
#'
#' Define a model, estimand, and answer strategy. will need @ import DeclareDesign
#'
#' @inheritParams CQtools_internal_inherit_params
#' @param n Number of observations
#' @param analysis_model model used for analysis. updating is based on priors from this model so any prior data based updating is communicate by re-providing prior data.
#' @param param_type parameter types to use for data draw; defaults to prior
#' @param prior_data if prior data available this gets added to new data 
#' @param ... arguments passed to \code{query_model}
#'
#' @importFrom dplyr %>%
#' @import Rcpp
#' @export
#' @examples
#' require("DeclareDesign")
#'
#' my_design <- causal_model_designer(
#'   reference_model = make_model("X -> Y") %>%
#'     set_priors(c(1,1, 1, 1, 7, 1)),
#'   analysis_model = make_model("X -> Y"),  # Prior model
#'   n = 100,
#'   query = list(ATE = "Y[X=1] - Y[X=0]"),
#'   given = "X==1 & Y==1")
#'
#' df <- DeclareDesign::draw_data(my_design)
#' estimand <- draw_estimands(my_design)
#'
#' # Estimation and diagnosis
#'# options(mc.cores = parallel::detectCores())
#' estimate <- draw_estimates(my_design)
#' diag <- diagnose_design(my_design, sims = 2)
#' diag
#' sim <- simulate_design(my_design, sims = 2)
#'

causal_model_designer <- function(
	reference_model = make_model("X -> Y"),
	param_type = "prior_draw", 
	prior_data = NULL,
	analysis_model  = reference_model,
	n = 1,
	query,
	...
) {

	args <- list(...)
	arg_names <- names(args)


	data_args  <- arg_names %in% c("vars", "probs", "n_steps", "probs", "subsets")
	query_args <- arg_names %in% c("given", "stats", "expand_grid")


	# Data step
	data_step <-
	    DeclareDesign::declare_population(data =
			{
			reference_model <- set_parameters(reference_model, param_type = "prior_draw")
			data <- do.call(
			    make_data, 
			    c(model = list(reference_model), n = n, args[data_args])
			    )
			if(!is.null(prior_data)) data <- bind_rows(data, prior_data)
			attr(data, "parameters") <- get_parameters(reference_model)
			data}
		)

	# Inquiry
	estimand <- DeclareDesign::declare_estimand(handler = function(data) {
		reference_model <- set_parameters(reference_model, parameters = attr(data, "parameters"))
		value <- do.call(query_model, 
		                 c(model = list(reference_model),
		                   query = list(query),
		                   using = "parameters",
		                   args[arg_names %in% query_args]))
		
		names(value)[names(value) == "mean"] <- "estimand"
		names(value)[names(value) == "Query"] <- "estimand_label"
		value <- value[c("estimand_label", "estimand")]
		value})



	# Estimator runs CQ updated model assuming answer-strategy model
	estimate  <- DeclareDesign::declare_estimator(handler = function(data) {
		updated <- update_model(model = analysis_model,  data = data)
		value   <- do.call(query_model, c(model = list(updated),
																		query = list(query),
																		using = "posteriors",
																		args[arg_names %in% query_args]))

		names(value)[names(value) == "mean"] <- "estimate"
		names(value)[names(value) == "Query"] <- "estimand_label"
		names(value)[names(value) == "sd"]  <-  "std.error"
		value$estimator_label <- paste0("est_", value$estimand_label)
		value
	})

	# Declare design
	design <- data_step + estimand + estimate

	attr(design, "reference_model") <- reference_model

	class(design) <- c(class(design), "causal_model_design")
	
	my_diagnosands <- DeclareDesign::declare_diagnosands(select = c(mean_estimate, sd_estimate, mean_estimand, bias),
																				MSE = mean((estimate - estimand)^2),
																				posterior_var = mean(std.error^2))

	DeclareDesign::set_diagnosands(design, diagnosands = my_diagnosands)

}


#' @export
plot.causal_model_design <- function(x, ...) {
    CausalQueries:::plot_dag(attr(x, "reference_model"))
}

