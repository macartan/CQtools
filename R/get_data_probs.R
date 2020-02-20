#' Get data probabilities
#'
#' Takes in a matrix of possible (single case) observations and returns the probability of each.
#' FLAG: Some redundancy with make_data_probabilities
#'
#' @param model A  model
#' @param data Data in long format
#' @param parameters A numeric vector. Values of parameters may be specified. By default, parameters is drawn from priors.
#' @export
#' @examples
#' model <- make_model("X->Y")
#' data <- simulate_data(model, n = 4)
#' get_data_probs(model, data)
get_data_probs <- function(model, data, parameters = NULL){

	if(is.null(parameters)) parameters <- get_parameters(model)

	events  <- collapse_data(data = data, model = model)$event
	A_w     <- get_data_families(model, mapping_only = TRUE, drop_all_NA = FALSE)
	probs   <- A_w %*% get_event_prob(model, parameters = parameters)
	np      <- rownames(probs)
	unlist(sapply(events, function(j) probs[np==j]))
}


