#' helper for ways to allocate N units into n data types: tidies partition::composition output
#'
#' @param N Number of observations to be distributed
#' @param n Number of possible values observations could take
#' @examples
#' allocations(4,2)
allocations <- function(N, n) {
	x <- partitions::compositions(N,n)
	x <- data.frame(as.matrix(x))
	colnames(x) <- 1:ncol(x)
	x
}

#' Helper: order event data and add strategy family
#' @param df dataframe with event data
#' @param model a model made by make_model
#' @examples
#' model <- make_model("X -> M -> Y")
#' df <- data.frame(M = c(1, NA), X = c(1,1), Y = c(NA,NA)) %>%
#'  collapse_data(model) %>%
#'  filter(count == 1)
#' df
#' gbiqqtools:::check_event_data(df, model)
#'
check_event_data <- function(df, model) {
	if(!(names(df)[[1]] == "event")) stop("event_data must include an initial `event` column")
	if("strategy" %in% names(df)) df <- dplyr::select(df, - strategy)
	structure <- collapse_data(expand_data(df[, 1:2], model), model)[, 1:2]
	out <- dplyr::left_join(structure, df, by = "event")
	out[is.na(out)] <- 0
	out
	}



#' helper to fill buckets dataframe
#' @param buckets dataframe with columns event, count and capacity vars plus strategy allocation var
#' @param vars vars to be observed
#' @export
#' @examples
#' model <- make_model("X->M->Y")
#' buckets = data.frame(event = "X0Y0", count = 3, capacity = 3, strategy = 2)
#' # Find different data that might result from looking at "M" in 2 out of 3 X0Y0 data types
#' fill_bucket(model, buckets, vars = "M")
fill_bucket <- function(model, buckets, vars, row = 1, column = 4){

	if(!(all(vars %in% model$node))) stop("Vars not in model$node")

	# Figure out set of possible finer units
	df <- expand_data(data_events = data.frame(
												event = buckets$event[row], count = 1), model)
	possible_findings <- perm(rep(1, length(vars)))
	df <- df %>% slice(rep(1:n(), each = nrow(possible_findings)))
	df[vars] <- possible_findings
	df <- collapse_data(df, model, remove_family = TRUE)
	# Assign n across new possible finer events
	new_events <- cbind(event = df[df$count ==1, "event"],
											gbiqqtools:::allocations(buckets[row, column], sum(df$count)))

	# tidy up
	remaining  <- data.frame(event = buckets[row, 1], matrix(buckets$count[row] - buckets[row, column], ncol = ncol(new_events)-1, nrow = 1))
	names(remaining) <- names(new_events)
	rbind(new_events,remaining)
}


#' Helper for getting all data on specified node with N observed
#'
#' @param model model made with gbiqq::make_model
#' @param N Integer, number of observed cases
#' @param vars String vector listing node observed
#' @param condition Statement indicating condition satisfied by observed data
#' @examples
#' model <- make_model("X->M->Y")
#' gbiqqtools:::all_possible(model, N=2, vars = c("X", "M"))
#' gbiqqtools:::all_possible(model, N=2, vars = c("X", "Y"), condition = "Y==0")
all_possible <- function(model, N, vars = NULL, condition = TRUE){

	if(is.null(vars)) vars <- model$node

	df <- get_max_possible_data(model) %>% filter(eval(parse(text = condition)))

	if(!all(is.na(vars))) df[, !names(df) %in% vars] <- NA

	df  <- collapse_data(df, model, remove_family = TRUE)

	possible_data <- gbiqqtools:::allocations(N, n = sum(df$count>0))

	out <- matrix(0, nrow(df), ncol(possible_data))
	out[df$count > 0,] <- as.matrix(possible_data)
	cbind(df, out)
}


#' Encode data
#'
#' Takes data in long format, including NA values or blanks and returns vector with each row encoded as a data type.
#'
#' @param model A  model
#' @param data Data in long format
#' @export
#' @examples
#' model <- make_model("X -> Y")
#' data <- simulate_data(model, n = 4)
#' data[1,1] <- ""
#' data[3,2] <- NA
#'
#' encode_data(model, data)
encode_data <- function(model, data){
	data[data ==""] <- NA
	vars <- model$node
	apply(data, MARGIN = 1, FUN = function(row){
		paste0(vars[!(is.na(row))],row[!(is.na(row))], collapse = "")})
}

