#' Generates a probability distribution over possible data outcomes
#'
#' NOTE: This needs to be checked for whether it is taking account of strategy probabilities properly
#'
#' @param model A causal model as created by \code{make_model}
#' @param given A data frame with observations
#' @param subset data strategy
#' @export
#' @return A dataset
#' @examples
#'
#' library(dplyr)
#' model <- make_model("X->M->Y") %>%
#'          set_parameters(type = "flat")
#' possible_data <- make_possible_data(model, N= 2, vars = list(model$variables), within = FALSE)
#' make_data_probabilities(model, pars = model$parameters, possible_data)
#'
#' given <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1)) %>%
#'   collapse_data(model, remove_family = TRUE)
#' possible_data <- make_possible_data(model, given = given, condition = "X==1 & Y==1", vars = "M", within = TRUE )
#' make_data_probabilities(model, pars = model$parameters, possible_data)
#'
make_data_probabilities <- function(model, pars,  possible_data, A_w = NULL, strategy = NULL, strategy_set = NULL) {

	# Ambiguity matrix for data types
	if(is.null(A_w)) A_w <- (get_likelihood_helpers(model)$A_w)[possible_data$event, ]

	w_full = A_w %*% (draw_event_prob(model, parameters = pars, using = "parameters"))

	# Flag: not very elegant. Better maybe to carry strategy along with possible data
  strategy <- possible_data$strategy
  if(is.null(strategy)){
	strategy <- merge(possible_data[,1:2], collapse_data(expand_data(possible_data[, 1:2], model), model), by = "event")$strategy}

	if(is.null(strategy_set)) strategy_set <- unique(strategy)

	# Probability of outcomes within each strategy set
  # Need to be sure about ordering of data
	x <- apply(dplyr::select(possible_data, - event), 2, function(d)
		sapply(strategy_set, function(j) dmultinom(d[strategy==j],
																						prob = w_full[strategy==j])
		))
	if(!is.null(nrow(x))) x <- apply(x, 2, prod)

	# Normalization
	x/sum(x)
}




#' Generates an "average" probability distribution over possible data outcomes
#'
#' NOTE: This needs to be checked for whether it is taking account of strategy probabilities properly
#'
#' @param model A causal model as created by \code{make_model}
#' @param data_event A data frame with event plus one or multiple data columns
#' @export
#' @return A dataset
#' @export
#' @examples
#' library(dplyr)
#' model <- make_model("X->M->Y")
#' possible_data <- make_possible_data(model, vars = list(c("X", "Y")), N= 1, within = FALSE)
#' average_data_probabilities(model, possible_data, using = "priors")
#'
#'
average_data_probabilities <- function(model, data_event, using, sims = 500) {

	event <- data_event$event

	x <- replicate(sims, {
	pars <-	draw_parameters(model, using = using)
  probs <- make_data_probabilities(
		model,
		pars = pars,
		possible_data = data_event
		)
		}
	)
	apply(x, 1, mean)
 }

