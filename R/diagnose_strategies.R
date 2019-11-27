#' Diagnose a data strategy
#'
#' @param reference_model A causal model as created by \code{make_model}
#' @param analysis_model A causal model as created by \code{make_model}
#' @param given A data frame with existing data
#' @param queries Vector of causal statements characterizing queries
#' @param subsets Vector of statements refining queries
#' @param expand_grid Logical, expands grid over query arguments (queries, subsets)
#' @param data_strategies list containing arguments for data strategies.
#' For instance \code{list(strategy1 = list(N=1, withins = TRUE, vars = NULL, conditions = list(TRUE)))}
#' @param use_parameters Logical, defaulting to FALSE. If TRUE use  parameter vector rather than priors/posteriors.
#' Used in process tracing problems with a single case.
#' In this case the estimand is taken to be case level and is drawn  *conditional on case data*,
#' and the posterior variance is defined on the type, not the parameter (which is known, after all).
#' @param sims = 1000,
#' @param estimands_database Database of estimands, optional, for speed
#' @param estimates_database Database of estimates, optional, for speed
#' @param possible_data_list Database of possible data, optional, for speed
#' @export
#' @return A dataframe
#' @examples
#'
#' fit <- fitted_model()
#'
#' # Example using parameters and  minimal arguments, assumes search for one case
#' # But nothing learned about parameters
#' diagnosis <- diagnose_strategies(
#'   analysis_model = make_model("X->Y"),
#'   queries = "Y[X=1]==1",
#'   use_parameters = TRUE,
#'   fit = fit)

#' # Example with minimal arguments, assumes search for one case, updating on parameters
#' diagnosis <- diagnose_strategies(
#'   analysis_model = make_model("X->Y"),
#'   queries = "Y[X=1]==1",
#'   sims = 10)
#'
#'# Example comparing two data strategies with two queries
#'  rm(list = ls())
#'	if(!exists("fit")) fit  <- fitted_model()
#'
#'	analysis_model <-
#'    make_model("X->M->Y")  %>%
#'    set_restrictions(c("Y[M=1]<Y[M=0]", "M[X=1]<M[X=0]")) %>%
#'    set_parameter_matrix()
#'
#' given   <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1)) %>%
#'              collapse_data(analysis_model, remove_family = TRUE)
#'
#' queries <- list(ATE = "Y[X=1]-Y[X=0]", PC = "Y[X=1]-Y[X=0]")
#'
#' subsets <- list(TRUE, "Y==1 & X==1") # Subsets for queries
#'
#' data_strategies <- list(
#'   strategy1 = list(N=1, withins = TRUE, vars = "M", conditions = list("Y==1 & X==1")),
#'   strategy2 = list(N=1, withins = TRUE, vars = "M", conditions = list("Y==X")))
#'
#' diagnosis <- diagnose_strategies(
#'   analysis_model = analysis_model,
#'   data_strategies = data_strategies,
#'   given = given,
#'   queries = queries,
#'   subsets = subsets,
#'   sims = 4000,
#'   fit = fit)
#' diagnosis
#'
#' # Wide or deep illustration
#'
#' rm(list = ls())
#'
#' model <- make_model("K-> X -> Y <- K")
#' data_strategies = list(
#' 		N4L0 =  list(N=4, withins = FALSE, vars = list(c("X", "Y")), conditions = TRUE),
#' 		N2L2 =  list(N=2, withins = FALSE, vars = list(c("X", "K", "Y")), conditions = TRUE),
#' 		N3L1 =  list(N=list(1,2), withins = FALSE, vars = list(c("X", "K", "Y"), c("X", "Y")), conditions = TRUE))
#'
#' 	possible_data_list = lapply(data_strategies, function(ds)
#' 		with(ds, make_possible_data(model = model, given = NULL,
#' 		N = N, withins = withins, conditions = conditions, vars = vars)))
#'
#' lapply(possible_data_list, length)
#'
#' diagnosis <- diagnose_strategies(
#'   analysis_model = model,
#'   data_strategies = data_strategies,
#'   queries = "Y[X=1] - Y[X=0]",
#'   fit = fit,
#'   possible_data_list = possible_data_list)
#' diagnosis
#'
#' # Simple illustration of updating on a probability given a uniform prior.
#' # MSE and Expected posterior variance  should be both 1/18
#' diagnosis <- diagnose_strategies(
#'   analysis_model = make_model("X -> Y"),
#'   data_strategies = list(
#' 		take_one =  list(N=1, withins = FALSE, vars = list(c("X")), conditions = TRUE)),
#'   queries = "X==1",
#'   fit = fit,
#'   sims = 6000)
#' diagnosis
#'
#'

diagnose_strategies <- function(reference_model = NULL,
																analysis_model,
																given = NULL,
																queries,
																subsets = TRUE, # Subsets for queries
																expand_grid = FALSE, # For queries
																data_strategies = list(strategy1 = list(N=1, withins = TRUE, vars = NULL, conditions = list(TRUE))),
																sims = 1000,
																iter = NULL,
																refresh = 1000,
																use_parameters = FALSE,
																estimands_database = NULL,
																estimates_database = NULL,
																possible_data_list = NULL,
																fit = NULL,
																add_MSE = TRUE # Add MSE
																){
	# Housekeeping
	if(is.null(fit) & !(use_parameters)) fit  <- gbiqq::fitted_model()
  if(is.null(iter)) iter = max(sims, 4000)

  if(use_parameters) sims <- 1  # If parameters are used there is no simulation

	# Add a strategy with no data at the top of the list
	data_strategies[["Prior"]] <- list(N=0, withins = FALSE, vars = NULL, conditions = list(TRUE))
	# data_strategies <- data_strategies[c(length(data_strategies), 1:(length(data_strategies) - 1))]

	# 1. REFERENCE MODEL
	# If not provided, reference model should be the analysis model, updated
	if(is.null(reference_model)) {

		if(is.null(given) | use_parameters) {
		  reference_model <- analysis_model
		} else {
			data <- expand_data(given, analysis_model)
		  reference_model <- gbiqq(analysis_model, data, stan_model = fit, iter = iter, refresh = 0)
	}}

	# 2. REFERENCE PARAMETERS DISTRIBUTION
	#########################################
	if(use_parameters) {

		using <- "parameters"
		param_dist <- data.frame(t(get_parameters(reference_model)))
	} else {

	if(is.null(given) | is.null(reference_model$posterior_distribution)) {

			using <- "priors"

			if(is.null(reference_model$prior_distribution)) reference_model <-  set_prior_distribution(reference_model, n_draws = 8000)

			if(nrow(reference_model$prior_distribution) < sims) {
				stop("model contains priors distribution with too few parameter draws")}

			params_to_use <- sample(1:(nrow(reference_model$prior_distribution)), sims)
			param_dist    <- reference_model$prior_distribution[params_to_use,]

		} else {

		using <- "posteriors"
		params_to_use <- sample(1:(nrow(reference_model$posterior_distribution)), sims)
		param_dist <- reference_model$posterior_distribution[params_to_use,]

		}
	}

	# 3 ESTIMANDS DATABASE: If not provided, this gets generated using reference model and queries
	##################################################################################

	if(is.null(estimands_database)) {
	# Distribution of causal types (just one draw if parameters = TRUE)
		type_distribution <- get_type_prob_multiple(reference_model, using = using)

	# posterior_draws* [sims*]  length(queries) matrix of estimands
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
				                     		 make_possible_data(model = reference_model,
				                     		 									 given = given,
				                     		 									 N = N,
				                     		 									 withins = withins,
				                     		 									 conditions = conditions,
				                     		 									 vars = vars)))
		}

  # 5 POSSIBLE DATA DISTRIBUTION: FOR EACH STRATEGY nrow(possible_data) * sims matrix of data probabilities
	##################################################################################################
  ## FLAG: THIS FUNCTION IS THE SLOWEST STEP: HOW TO SPEED UP?

	data_probabilities_list <-
		lapply(possible_data_list[-length(possible_data_list)], function(possible_data) {

			# Arguments generated  prior to applying apply for speed
			A_w <-  get_data_families(reference_model, drop_impossible = TRUE, drop_none = TRUE, mapping_only = TRUE)[possible_data$event, ]
			strategy     <- possible_data$strategy
			strategy_set <- unique(possible_data$strategy)

	p <- apply(param_dist, 1, function(pars) {
			make_data_probabilities(
				model = reference_model,
				pars = pars,
				possible_data = possible_data,
				A_w  = A_w,
				strategy = strategy,
				strategy_set = strategy_set)
			}
			)

	# p should be of dimensionality simulations x data_possibilities
	if(use_parameters) return(matrix(p, nrow = 1))
	return(t(p))

		})

	# For prior existing data is certain
	data_probabilities_list[["prior"]] <- matrix(1, nrow = sims)


	# 6 ESTIMATES DISTRIBUTION: ESTIMAND FOR EACH QUERY FROM UPDATED DATA
	##################################################################################################
	if(is.null(estimates_database)) {
		estimates_database <- lapply(possible_data_list, function(possible_data){

		        make_estimates_database(analysis_model,
		        												given = given,
		        												possible_data = dplyr::select(possible_data, -strategy),
														        queries = queries,
		        												subsets = subsets,
		        												expand_grid = expand_grid,
		        												iter = iter,
		        												use_parameters = use_parameters,
		        												refresh = refresh)
		})
	}

  # 7 MAGIC: Calculate diagnostics for each data strategy
	##################################################################################################
	estimands <- apply(data.frame(estimands_database), 2, mean) #estimand calculation from full distribution
	if(!use_parameters) mands  <- data.frame(estimands_database[params_to_use,]) # estimand draws from sampled parameters

diagnosis <-

	lapply(1:length(estimates_database), function(j) { # applied over each strategy, inlcuding prior

		mate <- estimates_database[[j]]
		# Error if each datatype observed  n_query * n_data_possibilities
		estimate <- sapply(1:length(mate), function(i)  mate[[i]]$mean)
		if(is.null(dim(estimate))) estimate <- matrix(estimate, nrow = 1) #eg use_parameters

		# post_var should be a n_data * nq matrix
		post_var <- sapply(1:length(mate), function(i) (mate[[i]]$sd)^2)
		if(length(queries)==1) post_var <- matrix(post_var, nrow =1)
		prob     <- data_probabilities_list[[j]]

		# Return
		df <- data.frame(
			strategy  = names(data_strategies)[j],
			estimates_database[[1]][[1]][,c("Query", "Subset")],    # This picks up query and subset labels
			estimand  = estimands,

			# prob is dim(sims * data_possibilities) estimate is dim(queries * data_possibilities)
			# if use_parametes 1*2, t(q * 2)
			estimates = apply(prob%*%t(estimate), 2, mean) # Double averaging: over parameters and data type draw
      )

		if(use_parameters)   df$post_var= (abs(df$estimate))*(1 - (abs(df$estimate))) # p(1-p) for probability of type

		if(!use_parameters)  {

				if(add_MSE){
				# Go through for each query
				df$MSE    <-

				sapply(1:ncol(mands), function(q) {

	          est_given_data <- estimate[q,]
						squared_error <- sapply(est_given_data, function(est) (est-mands[,q])^2)
													 	mean(apply(prob * squared_error, 1, sum))
													}
													)
				}


#			df$post_var = post_var%*%apply(prob,2, mean)
			Expected_post_vars <- prob%*%t(post_var)
			df$post_var = apply(Expected_post_vars, 2, mean)
			df$post_var_sd = apply(Expected_post_vars, 2, sd)
		}

		df
	})


# Clean up and export: Put priors in front and bind
diagnosis <- diagnosis[c(length(diagnosis), 1:(length(diagnosis) - 1))]
diagnosis <- do.call("rbind", diagnosis)
rownames(diagnosis) <- NULL

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

### NOTE TO SELF
## EPXETED POSTERIOR VARIANCE WHEN PARAMETERS = TRUE HAS TO BE ON THE *CASE* LEVEL VARIANCE NOT HTE PARAMETER VARIANCE, WHICH WILL BE 0!


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


