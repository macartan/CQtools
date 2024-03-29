#' Make data for multi-step strategy
#'
#' Creates a database of possible data from a data strategy.
#' Users can gather additional data on node specified via \code{vars} for any possible cases in the model ("any"). Or they can
#' gather data in all cases within an observed dataset ("within"). Or they can specify  the subset of cases for which within-case data should be collected (e.g. "Y == 1").
#'
#' @inheritParams CQtools_internal_inherit_params
#' @param N An integer. Number of node to seek.
#' @param withins A list of logicals. Whether to seek node within existing data. Defaults to TRUE.
#' @param conditions  A list of character strings indicating for which cases data should be gathered. Options are: (i) to gather additional data on node specified via \code{vars} for any possible cases in the model ("any"), (ii) to gather data in all cases within an observed dataset ("within"), or (iii) to specify the subset of cases for which within-case data should be collected (e.g. "Y == 1").
#' @param vars A character vector. Variables to be sought or NA. If NA \code{make_possible_data} gathers data on all node containing NA for the specified data strategy.
#' @param prefix for columns of output; useful if multiple dataframes are later merged
#' @param unique = TRUE  If same data is gathered via different routes it still only gets represented once
#' @export
#' @return A dataset with columns: event, strategy, plus possibly multiple cases profiles
#' @examples
#' model <- make_model("X->M->Y")  %>%
#'    set_restrictions(c("Y[M=1]<Y[M=0]"), "(M[X=1]<M[X=0] ") %>%
#'    set_parameter_matrix()
#' df <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1))
#' observed <- collapse_data(df, model)[, -2]
#'
#' # Complex *not within*
#' make_possible_data(model, N=list(1,2), withins = FALSE, vars = list(c("X", "Y"), c("X")), conditions = TRUE)
#'
#' # Complex, sequential: not within, then within
#' make_possible_data(model, N=list(1,1), observed = observed,
#'              withins = c(FALSE, TRUE),
#'              vars = list(c("X", "Y"), c("M")),
#'              conditions = TRUE)
#'
#' # Look for data on M for all possible cases in the observed data
#' make_possible_data(model, N = 0, within = FALSE)
#' make_possible_data(model, N = 2, within = FALSE)
#' make_possible_data(model, observed = observed, vars = "M", N = 2, conditions = c("X==Y"))
#' make_possible_data(model, observed = observed, vars = "M", N = list(1,1), conditions = list("X==Y", "X==Y"))
#'
#' # Not possible:
#' make_possible_data(model, observed, vars = "M", within = TRUE, N = 7)
#'
#' # Partly possible: only one step completed
#' model2 <- make_model("A -> B -> C -> D")
#' observed2 <- data.frame(A = c(0,0,0,1,1,1), B = NA, C = NA, D = c(0,0,1,0,1,1)) %>%
#'   collapse_data(model2, drop_family = TRUE)
#' make_possible_data(model2, observed2, vars = list("B", "C"), within = TRUE, N = list(1,1), conditions = list("A==D", "A==D & B==1"))
#'
#' # Within conditions
#' make_possible_data(model, observed, within = TRUE, N = 2, conditions = "X==1 & Y==1", vars = "M")
#'
#' # Look for data on K but not M
#'
#'
#' # From book
#'
#' model <- make_model("X->M->Y")  %>%
#'  set_restrictions(c("(Y[M=1]<Y[M=0])", "(M[X=1]<M[X=0])"))
#'
#' 	observed <-  collapse_data(data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1)),
#' 													model, drop_family = TRUE)
#' 	make_possible_data(model,
#'  	observed = observed,
#'  	vars = list("M", "M", "M", "M"),
#'  	withins = TRUE,
#'  	N = list(1,1,1,1),
#'  	conditions = list("X==0 & Y==0", "X==1 & Y==0", "X==0 & Y==1", "X==1 & Y==1"))

make_possible_data <- function(
    model,
    observed = NULL,
    N = list(1),
	withins = TRUE,
    conditions = list(TRUE),
    vars = NULL,
    prefix = NULL,
    unique = TRUE) {

	if(is.null(vars)) vars <- list(model$node)

	if(is.null(observed) & withins[1]) {message("No datan observed; 'withins' changed to FALSE"); withins[1] <- FALSE}

	if(!is.null(observed)) if(!identical(names(observed), c("event", "count"))){
		stop("'observed' df should have two columns: event and count")}

	if(is.null(observed)) observed <- CausalQueries:::minimal_event_data(model)

	if("strategy" %in% names(observed)) observed <- dplyr::select(observed, - strategy)

	if(length(vars)==1 & (length(N)>1)) vars <- rep(vars, length(N))
	if(length(withins)==1 & (length(N)>1)) withins <- rep(withins, length(N))
	if(length(conditions)==1 & (length(N)>1)) conditions <- rep(conditions, length(N))
	if(!identical(length(conditions), length(N)) )
		stop("N, and conditions  must have the same length")



		if(!identical(length(vars), length(N)) )
		stop("Vars should be of length 1 or else have the same length as conditions  and N")

	g_df <- CQtools:::make_possible_data_single(
		model,
		observed = observed,
		withins = withins[[1]],
		N = N[[1]],
		conditions = conditions[[1]],
		vars = vars[[1]] )


	if(length(N) == 1){
		attr(g_df, "possible_data_args") <- list(N = N,withins = withins, conditions = conditions, vars = vars)
		names(g_df)[1:2] <- c("event", "count")
		g_df <- (CQtools:::check_event_data(g_df, model))
		colnames(g_df)[-c(1:2)] <- 1:(ncol(g_df)-2)

		return(g_df)
	}


	for (i in 2:length(N)){

		out <- lapply(2:ncol(g_df), function(s) {
			use_df <-  g_df[,c(1,s)]
			names(use_df)  <- c("event", "count")
			data_single <- make_possible_data_single(model,
																							 observed = 	use_df,
																							 withins = withins[[i]],
																							 N     = N[[i]],
																							 conditions = conditions[[i]],
																							 vars  = vars[[i]])

			colnames(data_single)[2:ncol(data_single)] <- paste0(s-1, "-", colnames(data_single)[2:ncol(data_single)])
			data_single
		})

		out   <- Reduce(function(x, y) merge(x, y,  by = c("event"), all = TRUE), 	out)
		g_df  <- dplyr:::mutate_if(out, is.numeric, ~replace(., is.na(.), 0))

	}

	if(!is.null(prefix)) names(g_df)[-1] <- paste0(prefix, "_", names(g_df)[-1])

	# Flag
	names(g_df)[1:2] <- c("event", "count")
	g_df <- check_event_data(g_df, model)

	g_df[,!duplicated(t(g_df))]

	colnames(g_df)[-c(1:2)] <- 1:(ncol(g_df)-2)
	
	if(unique) g_df <- g_df[, !duplicated(t(g_df))]

	g_df
}




#' Make possible data for a single strategy step
#'
#' Creates a database of possible data from a data strategy.
#' Users can gather additional data on node specified via \code{vars} for any possible cases in the model ("any"). Or they can
#' gather data in all cases within an observed dataset ("withins"). Or they can specify  the subset of cases for which withins-case data should be collected (e.g. "Y == 1").
#' @keywords internal
#' @inheritParams CQtools_internal_inherit_params
#' @param N Number of node to seek
#' @param withins logical Whether to seek node within existing data
#' @param conditions  A list of character strings indicating for which cases data should be gathered. Options are: (i) to gather additional data on node specified via \code{vars} for any possible cases in the model ("any"), (ii) to gather data in all cases within an observed dataset ("within"), or (iii) to specify the subset of cases for which within-case data should be collected (e.g. "Y == 1").
#' @param vars Variables to be sought or NA. If NA \code{make_possible_data} gathers data on all node containing NA for the specified data strategy.
#' @export
#' @return A dataset
#' @examples
#' model <- make_model("X->M->Y")  %>%
#'    set_restrictions(c("Y[M=1]<Y[M=0]", "M[X=1]<M[X=0]")) %>%
#'    set_parameter_matrix()
#' df <- data.frame(X = c(0,0,0,1,1,1), M = NA, Y = c(0,0,1,0,1,1))
#' observed <- collapse_data(df, model, drop_family = TRUE)
#'
#' make_possible_data_single(model, observed = observed, vars = "M", withins = TRUE, N = 2)

make_possible_data_single <- function(model,
																			observed = NULL,
																			N = 1,
																			withins = FALSE,
																			conditions = TRUE,
																			vars = NULL) {

	if(is.null(vars) & withins) stop("Please specify vars to be examined")

	if(withins & is.null(observed)) stop("If 'withins' is specified 'observed' must be provided")

	if(is.null(observed)) observed <- CausalQueries:::minimal_event_data(model)[,-2]

	if(N==0 & withins) return(observed)

	# If not withins simply select possible data
	if(!withins) return(
		CQtools:::complex_combine(list(
			observed,
		  new = CQtools:::all_possible(model, N, vars) %>% dplyr::select(-count))))


	# Otherwise its more complicated: select from *within* available data
	if(is.null(observed)) stop("observed not provided, but 'withins' requested")

	# all_event_types <- collapse_data(all_data_types(model), model, drop_family = TRUE)

	possible <- all_data_types(model) %>%
		filter(eval(parse(text = conditions))) %>%
		mutate(event = as.character(event))

	# This part to allow searching in cases where there is space to seek listed vars
	possible <- possible[apply(possible[vars], 1, function(j) all(is.na(j))),]

	all_buckets <- left_join(dplyr::select(possible, event), observed, by = "event")[c("event", "count")]
	all_buckets$count[is.na(all_buckets$count)] <-0
	all_buckets <- mutate(all_buckets, capacity = count)

	if(sum(all_buckets$capacity) < N) {message("Not enough units to allocate N.
																						 Perhaps you are seeking data within cases in which data is already observed?"); return(observed)}

	strategies <- as.matrix(partitions::blockparts(all_buckets$capacity, N))
		colnames(strategies) <- 1:ncol(strategies)
		all_buckets <- cbind(all_buckets, strategies)

		# This function goes through a bucket strategy and generates all possible datasets
		# that could be produced by the strategy
		##################################################################################
		get_results_from_strategy <- function(strategy){
			buckets          <- all_buckets[all_buckets$capacity>0 ,]
			buckets          <- buckets[buckets[,strategy]>0, c(1:3, strategy)]
			data_list        <- lapply(1:nrow(buckets), function(j)   fill_bucket(model, buckets, vars, row = j, column = 4))

      # If cases have been drawn from within set, remove these now
			if(withins) data_list[["remove_bucket"]] <- data.frame(event = buckets$event, x = -(buckets$capacity))

			data_list[["observed"]] <- observed

			out <- CQtools:::complex_combine(data_list)

			names(out)[-1] <-paste0(strategy-3, ".",  1:(ncol(out)-1))
			out

		}

		# Run over all strategies
		all_strategies <- sapply(4:ncol(all_buckets), function(s) get_results_from_strategy(s),
														 simplify = FALSE)
		out <- Reduce(function(x, y) merge(x, y,  by = "event", all = TRUE), all_strategies) %>%
			mutate_if(is.numeric, ~replace(., is.na(.), 0))

		out$event <- as.character(out$event)

  	out

  }


#' Complex combine
#'
#' Used to combine permutations of rows of dataframes in a list
#' @param data_list list of dataframes. All dataframes should contain event column but have unique event elements
#' @examples
#' data_list <- list(
#' data.frame(event = c("a", "b"), w = 1:2, x = 3:4),
#' data.frame(event = c("c", "d", "e"), y = 5:7, z = 8:10, q = 11:13))
#' complex_combine(data_list)
#' data_list <- list(
#' data.frame(event = c("a", "b"), w = 1:2, x = 3:4),
#' data.frame(event = c("c", "d", "b"), y = 5:7, z = 8:10, q = 11:13))
#' complex_combine(data_list)


complex_combine <- function(data_list) {

	locations <- CausalQueries:::perm(unlist(lapply(data_list, ncol)) - 2)

	dfs <- lapply(1:nrow(locations), function(i) {
		parts <-  lapply(1:length(data_list), function(j) {
			df <- data_list[j][[1]][, c(1, locations[i,j][[1]]+2)]
			names(df) <- c("event", "cases")
			df
		})
		out <- do.call("rbind", parts)
		out <- aggregate(cases ~ event, data = out, sum)
		names(out) <- c("event", i)
		out
	})
		Reduce(function(x, y) merge(x, y,  by = "event", all = TRUE), dfs)
 }

