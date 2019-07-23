#' Make  data for multi-step strategy
#'
#' Creates a database of possible data from a data strategy.
#' Users can gather additional data on variables specified via \code{vars} for any possible cases in the model ("any"). Or they can
#' gather data in all cases within a given dataset ("within"). Or they can specify  the subset of cases for which within-case data should be collected (e.g. "Y == 1").
#'
#' @param model A causal model object generated by \code{make_model}.
#' @param given A data.frame with observations.
#' @param N An integer. Number of variables to seek.
#' @param within Logical. Whether to seek variables within existing data.
#' @param condition  A list of character strings indicating for which cases data should be gathered. Options are: (i) to gather additional data on variables specified via \code{vars} for any possible cases in the model ("any"), (ii) to gather data in all cases within a given dataset ("within"), or (iii) to specify the subset of cases for which within-case data should be collected (e.g. "Y == 1").
#' @param vars A character vector. Variables to be sought or NA. If NA \code{make_possible_data} gathers data on all variables containing NA for the specified data strategy.
#' @export
#' @return A dataset
#' @examples
#' library(dplyr)
#' model <- make_model("X->M->Y")  %>%
#'    set_restrictions(causal_type_restrict = "Y[M=1]<Y[M=0] | M[X=1]<M[X=0] ") %>%
#'    set_parameter_matrix()
#' df <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1))
#' given <- summarize_data(model, df)[, -2]
#'
#' # Look for data on M for all possible cases in the given data
#' make_possible_data(model, N = 2)
#' make_possible_data(model, given, vars = list("M"), within = TRUE, N = 2)
#'
#' # Not possible:
#' make_possible_data(model, given, vars = "M", within = TRUE, N = 7)
#'
#' # Within conditions
#' make_possible_data(model, given, within = TRUE, N = 2, condition = "X==1 & Y==1", vars = "M")
#' make_possible_data(model, given, within = TRUE, N = 3, condition = "Y==1", vars = "M")
#' make_possible_data(model, given, within = TRUE, condition = "X == 1 | Y == 1", vars = "M")
#'
#' # Look for data on K but not M
#' model <- make_model("X->M->Y <-K")   %>%
#'    set_parameter_matrix()
#' df <- data.frame(X = c(0,0,1,1,1), K = NA, M = NA, Y = c(0,0,0,1,1))
#' given <- summarize_data(model, df)[, -2]
#' make_possible_data(model, given, within = TRUE, N = 1, vars = "K")
#'
#' # Look for data on M when X = 1 and Y = 0
#' make_possible_data(model,
#'                    given,
#'                    condition =  "X == 1 & Y == 0",
#'                    vars ="M")
#'
#' model <- make_model("X->M->Y")   %>%
#'    set_parameter_matrix()
#' make_possible_data(model,
#'                    given = NULL,
#'                    N = list(3,1),
#'                    within = FALSE,
#'                    condition =  list(TRUE, "X == 1 & Y == 0"),
#'                    vars = list(c("X", "Y"), "M"))

make_possible_data <- function(model,
															 given = NULL,
															 N = list(1),
															 within = FALSE,
															 condition = list(TRUE),
															 vars = list(NULL)) {

	if(!is.null(given)) if(!identical(names(given), c("event", "count"))){
		stop("'given' df should have two columns: event and count")}

	if(!identical(length(condition), length(vars), length(N)) )
		stop("N, cases and vars must have the same length")

	g_df <- gbiqq:::make_possible_data_single(model,
																						given = given,
																						within = within,
																						N = N[[1]],
																						condition = condition[[1]],
																						vars = vars[[1]] )

	if(length(N) == 1){
		attr(g_df, "possible_data_args") <- list(N = N,within = within, condition = condition, vars = vars)
		return(g_df)
	}


	for (i in 2:length(N)){
		possible <- data.frame(select(g_df, event), count = 1)
		possible <- simulate_data(model, data_events = possible)
		possible<-  possible[with(possible, eval(parse(text = condition[[i]]))),]
		possible <- gbiqq:::collapse_data(possible, model)
		use_this <- g_df[g_df$event %in% possible$event[possible$count>0],] >= N[[i]]
		if(nrow(use_this) > 1) use_this <- apply(use_this, 2, any)
		use_data <-  g_df[,use_this]
		skip   <-  data.frame(dplyr::select(g_df, event, strategy), g_df[,!use_this ])
		if(ncol(skip)>2) names(skip)[3:ncol(skip)] <-  names(g_df)[!use_this]

		# Avoid running make_possible_single if there are possible data with enough space to alocate N
		if(all(!use_this[3:ncol(g_df)])){
			message(paste0("Not enough space to allocate N = ",N[[i]], " in step ",i))
			return(g_df)
		}

		out <- lapply(3:ncol(use_data), function(s) {
			use_data <-  use_data[,c(1,s)]
			names(use_data)  <- c("event", "count")
			data_single <- make_possible_data_single(model,
																							 given = 	use_data,
																							 within = TRUE,
																							 N     = N[[i]],
																							 condition = condition[[i]],
																							 vars  = vars[[i]])

			colnames(data_single)[3:ncol(data_single)] <- paste0(s-2, "-", colnames(data_single)[3:ncol(data_single)])
			data_single
		})
		# out <- lapply(1:length(out), function(n_s){ out1 <- out[[n_s]]
		# 	                                          colnames(out1)[3:ncol(out1)] <- paste0(n_s, "-", colnames(out1)[3:ncol( out1)])
		# 	                                          out1 })
		# rbind not working yet since output is of different length
		# x <- do.call("rbind", out)
		# x <- t(t(x)[!duplicated(t(x)),])
		# given <- cbind(given[,1:2], x)
		out   <- Reduce(function(x, y) merge(x, y,  by = c("event", "strategy"), all = TRUE), 	out)
		out   <- merge( skip, out,  by = c("event", "strategy"), all = TRUE)
		g_df  <- dplyr:::mutate_if(out, is.numeric, ~replace(., is.na(.), 0))

	}
	# out2 <- lapply(1:length(out2), function(n_s){
	# 	out1 <- out2[[n_s]]
	# colnames(out1)[3:ncol(out1)] <- paste0(n_s, "-", colnames(out1)[3:ncol( out1)])
	# out1 })
	#given <- Reduce(function(x, y) merge(x, y,  by = c("event", "strategy"), all = TRUE), 	out2, right = TRUE)
	#given <- given[,!duplicated(t(given))]
	#given




	return(g_df)
}




#' Make possible data for a single strategy step
#'
#' Creates a database of possible data from a data strategy.
#' Users can gather additional data on variables specified via \code{vars} for any possible cases in the model ("any"). Or they can
#' gather data in all cases within a given dataset ("within"). Or they can specify  the subset of cases for which within-case data should be collected (e.g. "Y == 1").
#'
#' @param model A causal model as created by \code{make_model}
#' @param given A data frame in compact form with first column indicating event type and second column indicating number of events of that type.
#' @param N Number of variables to seek
#' @param within logical Whether to seek variables within existing data
#' @param condition  A list of character strings indicating for which cases data should be gathered. Options are: (i) to gather additional data on variables specified via \code{vars} for any possible cases in the model ("any"), (ii) to gather data in all cases within a given dataset ("within"), or (iii) to specify the subset of cases for which within-case data should be collected (e.g. "Y == 1").
#' @param vars Variables to be sought or NA. If NA \code{make_possible_data} gathers data on all variables containing NA for the specified data strategy.
#' @export
#' @return A dataset
#' @examples
#' library(dplyr)
#' model <- make_model("X->M->Y")  %>%
#'    set_restrictions(causal_type_restrict = "Y[M=1]<Y[M=0] | M[X=1]<M[X=0] ") %>%
#'    set_parameter_matrix()
#' df <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1))
#' given <- summarize_data(model, df)[, -2]
#'
#' # Look for data on M for all possible cases in the given data
#' make_possible_data_single(model, N = 2)
#' make_possible_data_single(model, given = given, vars = "M", within = TRUE, N = 2)
#' make_possible_data_single(model, given = given,
#'                           within = TRUE, vars = "M",
#'                           N = 2,
#'                           condition = "X==1 & Y==1")
#'
#' model <- make_model("X -> M -> Y <- K")  %>%
#'    set_restrictions(causal_type_restrict = "(Y[M=1, K= .]<Y[M=0, K= .]) | M[X=1]<M[X=0] ") %>%
#'    set_parameter_matrix()
#' given <- data.frame(X = c(0,0,0,1,1,1), K = NA,  M = NA, Y = c(0,0,1,0,1,1)) %>%
#'          collapse_data(model = model)
#' make_possible_data_single(model, given = given,
#'                           within = TRUE,
#'                           N = 2,
#'                           condition = "X==1",
#'                           vars = "M")
#'
#' make_possible_data_single(model, given = given,
#'                           within = TRUE,
#'                           N = 2,
#'                           condition = "X==1 & Y==1",
#'                           vars = c("M", "K"))
#'
#'model <- make_model("X->M->Y")  %>%
#' set_restrictions(causal_type_restrict = "Y[M=1]<Y[M=0] | M[X=1]<M[X=0]") %>%
#' set_parameter_matrix()
#' given <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1)) %>%
#'          collapse_data(model)
#' make_possible_data_single(model,
#'                           given = given,
#'                           within = TRUE,
#'                           vars = "M",
#'                           N = 1,
#'                           condition = "X==1 & Y==1")
#' make_possible_data_single(model,
#'                           given = given,
#'                           within = TRUE,
#'                           vars = "M",
#'                           N = 1,
#'                           condition = "X==1")

make_possible_data_single <- function(model,
																			given = NULL,
																			N = 1,
																			within = FALSE,
																			condition = TRUE,
																			vars = NULL) {

	if(is.null(vars) & within) stop("Please specify vars to be examined")


	if(within & is.null(given)) stop("If 'within' is specified 'given' must be provided")

	if(!within){possible_data <-  gbiqq:::all_possible(model, N, vars)}

	if(within){

		if(is.null(given)) stop("given not provided, but 'within' requested")

		all_event_types <- dplyr::select(summarize_data(model, all_data_types(model)), event, strategy)

		possible <- get_max_possible_data(model)
		possible <- possible[with(possible, eval(parse(text = condition))),]
		possible <- gbiqq:::collapse_data(possible, model)
		A_w      <- get_likelihood_helpers(model)$A_w

		# What is the set of types in which we can seek new data
		acceptable_bucket <- (A_w %*% possible[,"count"])>0
		acceptable_bucket <-rownames(acceptable_bucket)[acceptable_bucket]

		all_buckets          <- given
		all_buckets$capacity <- all_buckets$count
		all_buckets$capacity[!(given$event %in% acceptable_bucket)] <-0

		if(sum(all_buckets$capacity) < N) {message("Not enough space to allocate N"); return(given)}
		strategies <- as.matrix(partitions::blockparts(all_buckets$capacity, N))
		colnames(strategies) <- 1:ncol(strategies)
		all_buckets <- cbind(all_buckets, strategies)

		# This function goes through a bucket strategy and generates all possible datasets that could be produced by the strategy
		get_results_from_strategy <- function(strategy){
			buckets          <- all_buckets[all_buckets$capacity>0 ,]
			buckets          <- buckets[buckets[,strategy]>0,]
			data_list        <- lapply(1:nrow(buckets), function(j)   fill_bucket(model, buckets, vars, row = j, column = strategy))
			b_names          <- sapply(1:nrow(buckets), function(b_i) (all_buckets$event %in% buckets$event[b_i])*(ncol(data_list[[b_i]])-2))
			addresses        <- do.call(rbind, lapply(apply(b_names, 2, perm),function(b)data.frame(b)+1))
			addresses        <- apply(addresses, 1, paste, collapse = ".")
			strategy_results <- Reduce(function(x, y) merge(x, y,  by = "event", all = TRUE), data_list)
			out              <- merge(given, strategy_results,  by = "event", all = TRUE )
			out              <- dplyr::mutate_at(out, vars(-c("event", "count")),  list(~ dplyr::coalesce(., count)))
			colnames(out)[3:ncol(out)] <- paste0(strategy-3, "-",addresses)
			# Hack --avoid column names duplicates
			dups <- colnames(out)[duplicated(colnames(out))]
			l_dups <- length(dups)
			if(l_dups > 0) colnames(out)[duplicated(colnames(out))] <- paste0(dups, "_",seq_len(l_dups))
			out[,-2]
		}

		# Run over all strategies
		all_strategies <- sapply(4:ncol(all_buckets), function(s) get_results_from_strategy(s), simplify = FALSE)
		all_strategies <-	Reduce(function(x, y) merge(x, y,  by = "event", all = TRUE), all_strategies)
		all_strategies <- dplyr:::mutate_if(all_strategies, is.numeric, ~replace(., is.na(.), 0))
		possible_data  <- merge(select(given, event), all_strategies,  by = "event", all = TRUE)

		# Add strategies
		possible_data <- merge(all_event_types, possible_data, by = "event")
		possible_data <- dplyr:::mutate_if(possible_data, is.numeric, ~replace(., is.na(.), 0))


	}
	return(possible_data)
}