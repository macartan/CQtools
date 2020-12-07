#' Generates a probability distribution over possible data outcomes
#'
#' NOTE: This needs to be checked for whether it is taking account of strategy probabilities properly
#' @inheritParams CQtools_internal_inherit_params
#' @param pars A parameter vector.
#' @param possible_data Possible events data
#' @param A_w Ambiguity matrix for data types, optional
#' @param strategy_set vector containing possible strategies, optional
#' @param normalize logical if TRUE probabilities are normalized to sum to 1
#' @export
#' @return A dataset
#' @examples
#'
#' model <- make_model("X->M->Y")
#' possible_data <- make_possible_data(model, N= 2, vars = list(model$node), within = FALSE)
#' make_data_probabilities(model, pars = get_parameters(model), possible_data)
#'
#' given <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1)) %>%
#'   collapse_data(model)
#' possible_data <- make_possible_data(model, observed = select(given, - strategy), condition = "X==1 & Y==1", vars = "M", within = TRUE )
#' make_data_probabilities(model, pars = get_parameters(model), possible_data)
#' make_data_probabilities(model, pars = get_parameters(model), given)
#'
make_data_probabilities <- function(
    model, 
    pars,  
    possible_data, 
    A_w = NULL, 
    strategy = NULL, 
    strategy_set = NULL, 
    normalize = FALSE,
    type_prob = NULL,
    event_probs = NULL,
    P = NULL
    ) {

	# Check data consistency
	if(!all(names(possible_data)[1:2] == c("event", "strategy"))) 
	    stop("possible_data should lead with event and strategy columns")

	# If only one possible data, return 1
	if(normalize & ncol(possible_data %>% select(-event, -strategy))==1) 
	    return(1)

	
	if (is.null(event_probs)){
    
    	 # Type probabilities
    	 if (is.null(type_prob))
    	     type_prob <- get_type_prob(model = model, P = P, parameters = as.numeric(pars))
    	 
    	 event_probs <- get_event_prob(model, parameters = as.numeric(pars), type_prob=type_prob)

    	 }
    
    # Ambiguity matrix for data types
    if(is.null(A_w)) 
        A_w <- CausalQueries:::get_data_families(model, drop_impossible = TRUE, drop_all_NA = FALSE, mapping_only = TRUE)[possible_data$event, ]
    
    # Event probabilities
	w_full <- A_w %*% event_probs

    if(is.null(strategy))
        strategy <- possible_data$strategy
	
	if(is.null(strategy_set)) 
	    strategy_set <- unique(strategy)

  # Probability of outcomes within each strategy set
  # Need to be sure about ordering of data
	x <- 
	    dplyr::select(possible_data, - c(event, strategy)) %>%
	    apply(2, function(d)
		sapply(strategy_set, function(j) dmultinom(d[strategy==j], prob = w_full[strategy==j])
		))
	
  # Take product of probabilities across strategies	
	if(length(strategy_set)>1) x <- apply(x, 2, prod)

	# Normalization
	if(normalize) x <- x/sum(x)

	x

}



