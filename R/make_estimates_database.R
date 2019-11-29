#' Generates a database of results using gbiqq over possible data
#'
#' This function can run many models and can take a long time depending on the size of possible data.
#'
#' @param model A causal model as created by \code{make_model}
#' @param given A compact data frame with observe data
#' @param possible_data A data frame with an events column and possible data columns (if a strategy columns is included it is ignored)
#' @param queries list of statements for causal queries
#' @param subsets list of statements for subsets for queries
#' @param expand_grid logical, If TRUE combinations of queries and subsets are expanded
#' @param iter Number of iterations for stan estimation, defaults to 4000
#' @export
#' @return A list with query output dataframes for each data strategy
#' @examples
#' model <- make_model("X->M->Y")  %>%
#'    set_restrictions(c("Y[M=1]<Y[M=0]", "M[X=1]<M[X=0]")) %>%
#'    set_parameter_matrix()
#'
#' given <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1)) %>%
#'         collapse_data(model, remove_family = TRUE)
#'
#'
#' possible_data <- make_possible_data(model, given, vars = "M", condition = "X==1 & Y==1")
#'
#' estimates_database <- make_estimates_database(
#'       model,
#'       given = given,
#'       possible_data = possible_data,
#'       queries = "Y[X=1]>Y[X=0]")
#'
#' estimates_database <- make_estimates_database(
#'       model,
#'       given = given,
#'       possible_data = possible_data,
#'       queries = c(ATE = "Y[X=1]-Y[X=0]", PC = "Y[X=1]>Y[X=0]"),
#'       subsets = c(TRUE, "Y==1 & X==1"))

make_estimates_database <- function(model,
																		given,
																		possible_data = NULL,
																		queries = "Y[X=1]>Y[X=0]",
																		subsets = TRUE,
																		expand_grid = FALSE,
																		iter = 4000,
																		refresh = 0,
																		use_parameters = FALSE,
																		...) {


	if(is.null(possible_data)) possible_data <- make_possible_data(model, given, ...)
	if("strategy" %in% names(possible_data)) possible_data <- dplyr::select(possible_data, -strategy)

	if(use_parameters){
	return(
		lapply(2:ncol(possible_data), function(j) {

		data.frame(
			query_model(model,
									queries = queries,
									using = "parameters",
									subset = subsets,
									expand_grid =
									expand_grid),
			data_pattern = j -1
		)})
		)
	  }


	# If parameters not used then model is updated using gbiqq

	if(!exists("fit")) fit  <- fitted_model()


	## Update model for each possible data type and query updated model
	## Note: 2 here only because of particular shape of possible data

	lapply(2:ncol(possible_data), function(j) {

		data_events <- possible_data[, c(1, j)]

		data <- expand_data(data_events, model)

		updated <- gbiqq::gbiqq(model = model, data = data, stan_model = fit, iter = iter, refresh = refresh)

		data.frame(
			query_model(updated, queries = queries, using = "posteriors", subsets = subsets, expand_grid = expand_grid),
			data_pattern = j -1
			)

	})

}
