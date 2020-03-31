#' Generates a probability distribution over possible data outcomes
#'
#' NOTE: This needs to be checked for whether it is taking account of strategy probabilities properly
#' @inheritParams CQtools_internal_inherit_params
#' @param pars A parameter vector
#' @param possible_data Possible events data
#' @param A_w Ambiguity matrix for data types, optional
#' @param strategy_set vector containing possible strategies, optional
#' @param normalize logical if TRUE probabilites are normalized to sum to 1
#' @export
#' @return A dataset
#' @examples
#'
#' model <- make_model("X->M->Y")
#' possible_data <- make_possible_data(model, N= 2, vars = list(model$node), within = FALSE)
#' make_data_probabilities(model, pars = get_parameters(model), possible_data)
#'
#' given <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1)) %>%
#'   collapse_data(model, drop_family = TRUE)
#' possible_data <- make_possible_data(model, observed = given, condition = "X==1 & Y==1", vars = "M", within = TRUE )
#' make_data_probabilities(model, pars = get_parameters(model), possible_data)
#'
make_data_probabilities <- function(model, pars,  possible_data, A_w = NULL, strategy = NULL, strategy_set = NULL, normalize = FALSE) {

	# Check data consistency
	if(!all(names(possible_data)[1:2] == c("event", "strategy"))) stop("possible_data should lead with event and strategy columns")

	# If only one possible data, return 1
	if(normalize & ncol(possible_data)==3) return(1)

	# Ambiguity matrix for data types
	if(is.null(A_w)) A_w <- 	get_data_families(model, drop_impossible = TRUE, drop_all_NA = FALSE, mapping_only = TRUE)[possible_data$event, ]

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
	if(normalize) x <- x/sum(x)

	x

}



