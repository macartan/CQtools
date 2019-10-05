#' Diagnose a data strategy
#'
#' @param reference_model A causal model as created by \code{make_model}
#' @param analysis_model A causal model as created by \code{make_model}
#' @param given A data frame with existing data
#' @param queries queries
#' @export
#' @return A dataframe
#' @examples
#'
#' library(dplyr)
#' # Example with minimal arguments, assumes search for one case
#' diagnose_strategies(
#'   analysis_model = analysis_model,
#'   queries = queries,
#'   sims = 10)
#'
#'# Example comparing two data strategies with two queries
#'  rm(list = ls())
#'	if(!exists("fit")) fit  <- fitted_model()
#'	analysis_model <-
#'    make_model("X->M->Y")  %>%
#'    set_restrictions(c("Y[M=1]<Y[M=0]", "M[X=1]<M[X=0]")) %>%
#'    set_parameter_matrix()
#' given   <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1)) %>%
#'              collapse_data(analysis_model, remove_family = TRUE)
#' queries <- list(ATE = "Y[X=1]-Y[X=0]", PC = "Y[X=1]-Y[X=0]")
#' subsets <- list(TRUE, "Y==1 & X==1") # Subsets for queries
#' data_strategies <- list(strategy1 = list(N=1, within = TRUE, vars = "M", conditions = list("Y==1 & X==1")),
#' 											   strategy2 = list(N=1, within = TRUE, vars = "M", conditions = list("Y==X")))
#'
#' diagnosis <- diagnose_strategies(
#'   analysis_model = analysis_model,
#'   data_strategies = data_strategies,
#'   given = given,
#'   queries = queries,
#'   subsets = subsets,
#'   sims = 400)


diagnose_strategies <- function(reference_model = NULL,
																analysis_model,
																given = NULL,
																queries,
																subsets = TRUE, # Subsets for queries
																expand_grid = FALSE, # For queries
																data_strategies = list(strategy1 = list(N=1, within = TRUE, vars = NULL, conditions = list(TRUE))),
																sims = 1000,
																estimands_database = NULL,
																estimates_database = NULL,
																possible_data_list = NULL){
	# Housekeeping
	if(!exists("fit")) fit  <- fitted_model()

	# Add a strategy with no data at the top of the list
	data_strategies[["Prior"]] <- list(N=0, within = TRUE, vars = NULL, conditions = list(TRUE))
	data_strategies <- data_strategies[c(length(data_strategies), 1:(length(data_strategies) - 1))]

	#' make_possible_data(model, N = 2, within = FALSE)

	# 1. REFERENCE MODEL
	# If not provided, reference model should be the analysis model, updated
	if(is.null(reference_model)) {

		if(is.null(given)) {
		  reference_model <- analysis_model
		} else {
			data <- expand_data(given, analysis_model)
		  reference_model <- gbiqq(analysis_model, data, stan_model = fit)
	}}

	# 2. REFERENCE PARAMETERS DISTRIBUTION
	#########################################

	if(is.null(given) | is.null(reference_model$posterior_distribution)) {
		using <- "priors"

		if(!is.null(reference_model$prior_distribution)) {if(nrow(reference_model$prior_distribution) > sims) {
			reference_model$prior_distribution <- reference_model$prior_distribution[1:sims,]}}

		if(!is.null(reference_model$prior_distribution)) {if(nrow(model$prior_distribution) < sims) {
			stop("model contains priors distribution with too few parameter draws")}}

		if(is.null(reference_model$prior_distribution)) reference_model <-  set_prior_distribution(reference_model, sims)
		param_dist <- reference_model$prior_distribution


	} else {

		using <- "posteriors"
		param_dist <- (rstan::extract(reference_model$posterior, pars= "lambdas")$lambdas)[1:sims,]

	}


	# 3 ESTIMANDS DATABASE: If not provided, this gets generated using reference model and queries
	##################################################################################

	if(is.null(estimands_database)) {
	# Distribution of causal types
		type_distribution <- draw_type_prob_multiple(reference_model, using = using)

	# sims * length(queries) matrix of estimands
	estimands_database <- mapply(query_distribution,
															 model  = list(reference_model),
															 query  = queries,
															 subset = subsets,
															 using  = using,
															 type_distribution = list(type_distribution))
	}

	# 4 POSSIBLE DATA: Generate if not provided
	##################################################################################
	# list of possible_data given different strategies: Can handle multiple data strategies
		if(is.null(possible_data_list)){

			possible_data_list = lapply(data_strategies, function(ds)
				                     with(ds,
				                     		 make_possible_data(model = reference_model, given = given, N = N, within = within, condition = conditions, vars = vars)))
		}

  # 5 POSSIBLE DATA DISTRIBUTION: FOR EACH STRATEGY nrow(possible_data) * sims matrix of data probabilities
	##################################################################################################
  ## FLAG: THIS FUNCTION IS THE SLOWEST STEP: HOW TO SPEED UP?
	## Speed up possible by skipping for the priors strategy as it

	data_probabilities_list <-
		lapply(possible_data_list, function(possible_data) {

			# Arguments generated one prior to applying apply; FLAG strategy definition is nasty
			A_w <- (get_likelihood_helpers(reference_model)$A_w)[possible_data$event, ]
			strategy <- merge(possible_data[,1:2], collapse_data(expand_data(possible_data[, 1:2], reference_model), reference_model), by = "event")$strategy
			strategy_set <- unique(strategy)

		apply(param_dist, 1, function(pars) {  # FLAG NEEDS TO WORK WITH POSTERIOR ALSO
		make_data_probabilities(
			reference_model,
			pars = pars,
			possible_data = possible_data,
			A_w  = A_w,
			strategy = strategy,
			strategy_set = strategy_set)})})


	# 6 ESTIMATES DISTRIBUTION: ESTIMAND FOR EACH QUERY FROM UPDATED DATA
	##################################################################################################
	if(is.null(estimates_database)) {
		estimates_database <- lapply(possible_data_list, function(possible_data){

		        make_estimates_database(analysis_model, given = given, possible_data = possible_data,
														queries = queries, subsets = subsets, expand_grid = expand_grid)
		})
	}

  # 7 MAGIC: Calculate diagnostics for each data strategy
	##################################################################################################
	estimands <- apply(estimands_database, 2, mean)

diagnosis <-
	lapply(1:length(estimates_database), function(j) {

		mate <- estimates_database[[j]]
		# Error if each datatype observed  n_query * n_data_possibilities
		estimate <- sapply(1:length(mate), function(i)  mate[[i]]$mean)
		post_var <- sapply(1:length(mate), function(i) (mate[[i]]$sd)^2)
		prob     <- data_probabilities_list[[j]]

		# Return
		data.frame(
			strategy = names(data_strategies)[j],
			estimates_database[[1]][[1]][,1:2],
			estimand = estimands,
			estimate = apply(estimate%*%prob, 1, mean), # Double averaging: over parameters and data type draw
			MSE      = apply((estimate%*%prob - (t(estimands_database)[,1:sims]))^2, 1, mean),
			post_var = apply(post_var%*%prob, 1, mean))
	})


# Clean up and export
diagnosis <- do.call("rbind", diagnosis)

return_list <- list()
return_list$diagnoses_df        <- diagnosis
return_list$possible_data_list  <- possible_data_list
return_list$estimands_database  <- estimands_database
return_list$estimates_database  <- estimates_database
return_list$data_probabilities_list  <-  data_probabilities_list

return_list

class(return_list) <- c("strategy_diagnoses")

return_list

}

###

#' @export
print.strategy_diagnoses <- function(x, ...) {
	print(summary(x))
	invisible(x)
}


#' @export
summary.strategy_diagnoses <- function(object, ...) {
	structure(object, class = c("summary.strategy_diagnoses", "data.frame"))

}

#' @export
print.summary.strategy_diagnoses <- function(x, ...){
	print(x$diagnoses_df)
}



#' @export
print.strategy_diagnosis_single <- function(x, ...) {
	print(summary(x))
	invisible(x)
}


#' @export
summary.strategy_diagnosis_single <- function(object, ...) {
	structure(object, class = c("summary.strategy_diagnosis_single", "data.frame"))

}

#' @export
print.summary.strategy_diagnosis_single <- function(x, ...){
	print(x$diagnoses_df)

}


