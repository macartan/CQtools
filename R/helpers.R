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
#' @inheritParams CQtools_internal_inherit_params
#' @param df dataframe with event data
#' @examples
#' model <- make_model("X -> M -> Y")
#' df <- data.frame(M = c(1, NA), X = c(1,1), Y = c(NA,NA)) %>%
#'  collapse_data(model) %>%
#'  filter(count == 1)
#' df
#' CQtools:::check_event_data(df, model)
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
#' @inheritParams CQtools_internal_inherit_params
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
	df <- collapse_data(df, model, drop_family = TRUE)
	# Assign n across new possible finer events
	new_events <- cbind(event = df[df$count ==1, "event"],
											CQtools:::allocations(buckets[row, column], sum(df$count)))

	# tidy up
	remaining  <- data.frame(event = buckets[row, 1], matrix(buckets$count[row] - buckets[row, column], ncol = ncol(new_events)-1, nrow = 1))
	names(remaining) <- names(new_events)
	rbind(new_events,remaining)
}


#' Helper for getting all data on specified node with N observed
#'
#' @inheritParams CQtools_internal_inherit_params
#' @param N Integer, number of observed cases
#' @param vars String vector listing node observed
#' @param condition Statement indicating condition satisfied by observed data
#' @examples
#' model <- make_model("X->M->Y")
#' CQtools:::all_possible(model, N=2, vars = c("X", "M"))
#' CQtools:::all_possible(model, N=2, vars = c("X", "Y"), condition = "Y==0")
all_possible <- function(model, N, vars = NULL, condition = TRUE, possible_data = TRUE, complete_data = TRUE){

	if(is.null(vars)) vars <- model$node

	#df <- all_data_types(model, complete_data = complete_data)%>%
	##df <- get_max_possible_data(model)
  #	filter(eval(parse(text = condition))) %>% select(-event)

	df <- all_data_types(model, possible_data = possible_data, complete_data = complete_data, given = condition)

	if(!all(is.na(vars))) df[, !names(df) %in% vars] <- NA

	df  <- collapse_data(df, model, drop_family = TRUE)

	possible_data <- CQtools:::allocations(N, n = sum(df$count>0))

	out <- matrix(0, nrow(df), ncol(possible_data))
	out[df$count > 0,] <- as.matrix(possible_data)
	cbind(df, out)
}


#' Encode data
#'
#' Takes data in long format, including NA values or blanks and returns vector with each row encoded as a data type.
#'
#' @inheritParams CQtools_internal_inherit_params
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
	out <- apply(data, MARGIN = 1, FUN = function(row){
		paste0(vars[!(is.na(row))],row[!(is.na(row))], collapse = "")})
	out[out == ""] <- "None"
	out
}





#' Get joint distribution of nodal types
#'
#' Identifies possible conditional probabilities of nodal types. May be used to identify patterns of non independence.
#'
#' @inheritParams CQtools_internal_inherit_params
#' @param generic_parameters Logical. Whether to require selection of a generic parameter. Defaults to TRUE.
#' @keywords internal
#' @examples
#' model <- make_model('X -> Y') %>%
#'   set_confound(list('X <-> Y'))
#'
#' get_nodal_joint_probability(model)
#'

get_nodal_joint_probability <- function(model, parameters = NULL, generic_parameters = TRUE) {

    parameters_df <- model$parameters_df

    # Figure parameters to use
    if (!is.null(parameters))
        parameters <- CausalQueries:::clean_param_vector(model, parameters)
    if (is.null(parameters)) {
        if (!generic_parameters)
            parameters <- get_parameters(model)
        if (generic_parameters)
            parameters <- clean_param_vector(model, runif(nrow(parameters_df)))
    }

    nodal_type <- parameters_df$nodal_type
    P <- data.frame(get_parameter_matrix(model))
    type_prob <- get_type_prob(model, parameters = parameters)

    joint <- sapply(unique(nodal_type), function(par1) {
        sapply(unique(nodal_type), function(par2) prob_par1_given_par2(par1, par2, nodal_type, P, type_prob))
    }) %>% data.frame(stringsAsFactors = FALSE)

    # Add in node and nodal_type
    select(parameters_df, node, nodal_type) %>% distinct %>% mutate(nodal_type = factor(nodal_type)) %>%
        right_join(cbind(nodal_type = (rownames(joint)), joint), by = "nodal_type")
}

#' helper to get conditional probability of nodal types
#'
#' @param par1 parameter 1
#' @param par2 parameter 2
#' @param nodal_type vector of nodal types. See \code{\link[CausalQueries]{get_nodal_types}}
#' @param type_prob vector of type probabilities
#'
prob_par1_given_par2 <- function(par1, par2, nodal_type, P, type_prob) {
    par1_in_type <- dplyr::filter(P, nodal_type %in% par1) %>% apply(2, function(j) any(j == 1))
    par2_in_type <- dplyr::filter(P, nodal_type %in% par2) %>% apply(2, function(j) any(j == 1))
    sum(type_prob[par1_in_type & par2_in_type])/sum(type_prob[par2_in_type])
}



#'combine two lists by names
#' # FLAG : CHECK EXAMPLE BELOW
#' @param list1 a list
#' @param list2 a list typically different than list1
#' @examples
#' list1 = list(A = 1:3, B = 4)
#' list2 = list(A = 1, C = 1:2)
#' CQtools:::combine_lists(list1, list2)
#'
combine_lists <- function(list1, list2) {

	matches <- names(list1) %in% names(list2)
	matching_names <- names(list1)[matches]
	matches2 <- names(list2) %in% names(list1)

	if (any(matches)) {
		combined_list <- sapply(matching_names, function(nam) {
			out <- c(list1[[nam]], list2[[nam]])
			out[!duplicated(names(out))]
		}, simplify = FALSE)

		combined_list <- c(combined_list, list1[!matches])
		combined_list <- c(combined_list, list2[!matches2])
	} else {
		combined_list <- c(list1, list2)
	}

	combined_list
}



# helper from hadley
# https://stackoverflow.com/questions/4752275/test-for-equality-among-all-elements-of-a-single-vector
zero_range <- function(x, tol = .Machine$double.eps^0.5) {
	if (length(x) == 1)
		return(TRUE)
	x <- range(x)/mean(x)
	isTRUE(all.equal(x[1], x[2], tolerance = tol))
}

