#' Generates a probability distribution over possible data outcomes
#'
#' NOTE: This needs to be checked for whether it is taking account of strategy probabilities properly
#'
#' @param model A causal model as created by \code{make_model}
#' @param pars A parameter vector
#' @param possible_data Possible events data
#' @param A_w Ambiguity matrix for data types, optional
#' @param strategy vector providing data strategy set for event, optional
#' @param strategy_set vector containing possible strategies, optional
#' @export
#' @return A dataset
#' @examples
#'
#' model <- make_model("X->M->Y")
#' possible_data <- make_possible_data(model, N= 2, vars = list(model$node), within = FALSE)
#' make_data_probabilities(model, pars = get_parameters(model), possible_data)
#'
#' given <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1)) %>%
#'   collapse_data(model, remove_family = TRUE)
#' possible_data <- make_possible_data(model, given = given, condition = "X==1 & Y==1", vars = "M", within = TRUE )
#' make_data_probabilities(model, pars = get_parameters(model), possible_data)
#'
make_data_probabilities <- function(model, pars,  possible_data, A_w = NULL, strategy = NULL, strategy_set = NULL) {

	# Check data consistency
	if(!all(names(possible_data)[1:2] == c("event", "strategy"))) stop("possible_data should lead with event and strategy columns")

	# If only one possible data, return 1
	if(ncol(possible_data)==3) return(1)

	# Ambiguity matrix for data types
	if(is.null(A_w)) A_w <- 	get_data_families(model, drop_impossible = TRUE, drop_none = TRUE, mapping_only = TRUE)[possible_data$event, ]

	w_full = A_w %*% (get_event_prob(model, parameters = as.numeric(pars)))

  if(is.null(strategy)){strategy <- possible_data$strategy}
	if(is.null(strategy_set)) strategy_set <- unique(strategy)

	# Probability of outcomes within each strategy set
  # Need to be sure about ordering of data
	x <- apply(dplyr::select(possible_data, - c(event, strategy)), 2, function(d)
		sapply(strategy_set, function(j) dmultinom(d[strategy==j], prob = w_full[strategy==j])
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
#' @param possible_data An events dataframe
#' @param using String, indicating `priors` or `posteriors`
#' @param sims Integer, number of simulations
#' @export
#' @return A dataset
#' @export
#' @examples
#' library(dplyr)
#' model <- make_model("X->M->Y")
#' possible_data <- make_possible_data(model, vars = list(c("X", "Y")), N= 1, within = FALSE)
#' average_data_probabilities(model, possible_data, using = "priors", sims = 10)
#'
#'
average_data_probabilities <- function(model, possible_data, using, sims = 500) {

  # Predefined for speed
	A_w <- (get_likelihood_helpers(model)$A_w)[possible_data$event, ]
	strategy     <- possible_data$strategy
	strategy_set <- unique(possible_data$strategy)

	x <- replicate(sims, {
	  make_data_probabilities(
			model,
			pars = draw_parameters(model, using = using),
			possible_data = possible_data,
			A_w  = A_w,
			strategy = strategy,
			strategy_set = strategy_set
			)
			}
		)
	apply(x, 1, mean)
 }
