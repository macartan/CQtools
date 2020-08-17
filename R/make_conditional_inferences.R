
#' Conditional inferences
#'
#' Calculate estimands conditional on observed data (currently, for single-case process tracing) together with data realization probabilities
#' Realization probabilities are the probability of the observed data given data is sought on observed node
#'
#' @inheritParams CQtools_internal_inherit_params
#'
#' @export
#' @examples
#' model <- make_model("X->Y")
#' model <- set_parameters(model, type = "flat")
#' conditional_inferences(model, query = "Y[X=1]>Y[X=0]")
#'
#' # Example of posteriors given monotonic X -> M -> Y model
#' model <- make_model("X-> M -> Y")  %>%
#'   set_restrictions(labels = list(M = "10", Y = "10"))
#' conditional_inferences(model, query = "Y[X=1]>Y[X=0]", given = "Y==1")
#'
#' # Running example
#' model <- make_model("S -> C -> Y <- R <- X; X -> C -> R") %>%
#'    set_restrictions(labels =
#'    list(C = "1110", R = "0001", Y = "0001"), keep = TRUE)
#' conditional_inferences(model, query = list(COE = "(Y[S=0] > Y[S=1])"),
#' given = "Y==1 & S==0")

conditional_inferences <- function(model, query, parameters=NULL,  given = NULL){

	if(is.null(parameters)) parameters <- get_parameters(model)

	vars <- model$nodes

	# Possible data
	vals <-  all_data_types(model, given = given) %>% select(-event)

	# Conditions
	conds <- t(apply(vals, 1, function(j) paste(vars, j, sep = "==")))
	conds[is.na(vals)] <- NA
	given <- apply(conds, 1, function(j) paste(j[!is.na(j)], collapse = " & ")) %>%
		 as.list
	given[given==""] <- paste0(model$nodes[1], ">-1") # Guaranteed true

	impossible <- lapply(given, function(s) all(!(CausalQueries:::map_query_to_causal_type(model, s)$types))) %>% unlist

	if(all(impossible)) return(data.frame(vars, posterior = NA, prob = NA))

	vals   <- vals[!impossible, ]
	given  <- given[!impossible]


	# Calculate estimands
	estimands <- query_model(
		model   = model,
		parameters  = parameters,
		using = "parameters",
		queries = query,
		given = given)$mean

	# Cac=lculate data probabilities
	probs <- unlist(get_data_probs(model, data = vals))
	probs <- probs[rownames(vals)]


	out <- data.frame(cbind(vals, estimands, probs))

	names(out) <- c(vars, "posterior", "prob")

	data.frame(out)
}
