#' Integrate CausalQueries with DeclareDesign
#'
#' Define a model, estimand, and answer strategy. will need @ import DeclareDesign
#'
#' @inheritParams CQtools_internal_inherit_params
#' @param n Number of observations
#' @param ... arguments passed to \code{query_model}
#'
#' @importFrom dplyr %>%
#' @import Rcpp
#' @export
#' @examples
#' require("DeclareDesign")
#'
#' my_design <- CQ_designer(
#'   reference_model = make_model("X -> Y") %>%
#'     set_priors(c(1,1, 1, 1, 7, 1)),
#'   analysis_model = make_model("X -> Y"),  # Prior model
#'   n = 100,
#'   query = list(ATE = "Y[X=1] - Y[X=0]"),
#'   given = "X==1 & Y==1")
#'
#' df <- draw_data(my_design)
#' mand <- draw_estimands(my_design)
#'
#' # Estimation and diagnosis
#' options(mc.cores = parallel::detectCores())
#' mate <- draw_estimates(my_design)
#' diag <- diagnose_design(my_design, sims = 2)
#' diag
#' sim <- simulate_design(my_design, sims = 2)
#'

CQ_designer <- function(
	reference_model = make_model("X -> Y"),
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
		declare_population(data =
			{
			reference_model <- set_parameters(reference_model, param_type = "prior_draw")
			data <- do.call(make_data, c(model = list(reference_model), n = n,
																				args[data_args]))
			attr(data, "parameters") <- get_parameters(reference_model)
			data}
		)

	# Inquiry
	estimand <- declare_estimand(handler = function(data) {
		reference_model <- set_parameters(reference_model, parameters = attr(data, "parameters"))
		value <- do.call(query_model, c(model = list(reference_model),
																		query = list(query),
																		using = "parameters",
																		args[arg_names %in% query_args]))[,c(1,4)]

		names(value) <- c("estimand_label", "estimand")
		value})



	# Estimator runs CQ updated model assuming answer-strategy model
	estimate  <- declare_estimator(handler = function(data) {
		updated <- update_model(model = analysis_model,  data = data)
		value   <- do.call(query_model, c(model = list(updated),
																		query = list(query),
																		using = "posteriors",
																		args[arg_names %in% query_args]))

		names(value)[c(1,4)]<- c("estimand_label", "estimate")
		value$estimator_label <- paste0("est_", value$estimand_label)
		value
	})

	# Declare design
	design <- data_step + estimand + estimate

	my_diagnosands <- declare_diagnosands(select = c(mean_estimate, sd_estimate, mean_estimand, bias),
																				MSE = mean((estimate - estimand)^2),
																				posterior_var = mean(sd^2))

	set_diagnosands(design, diagnosands = my_diagnosands)

}

