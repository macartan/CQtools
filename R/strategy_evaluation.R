#' Strategy evaluation
#'
#' @param strategies  s x 3 matrix of strstegies. Each strategy is a triple that indicates (a) first node sought (b) action if first node is 0 (c) action if first node = 1
#' @param query A causal query, for instance "(Y[S=1] != Y[S=0])"
#' @param given A statement describing known data; for instance "Y==0 & X==1",
#' @param prices A vector of prices of length length(model$nodes)
#'
#' This function calculates the expected posterior variance from each strategy and the expected number of clues sought
#'
#' @export
#' @examples
#' model <-
#'   make_model("S -> C -> Y <- R <- X; X -> C -> R") %>%
#'      set_restrictions(
#'      labels = list(C = c("C1110", "C1111"),
#'                    R = c("R0001", "R0000"),
#'                    Y = c("Y0001")), keep = TRUE)
#'
#' # Evaluation of single strategy
#' strategy_evaluation(model,
#' 		strategies  = c("S", NA, NA),  given = "Y==0",
#' 		query = "(Y[S=1] != Y[S=0])")
#'
#' strategy_evaluation(model,
#' 		strategies  = c("S", "Y", "Y"),  given = "Y==0",
#' 		query = "(Y[S=1] != Y[S=0])")
#'
#' strategy_evaluation(model,
#' 		strategies  = c("S", NA, "X"),  given = "Y==0",
#' 		query = "(Y[S=1] != Y[S=0])")
#'
#'\dontrun{
#' strategy_evaluation(model,
#' 		strategies  = c("X", NA, NA),  given = "Y==0 & R==1",
#' 		query = "(Y[S=1] != Y[S=0])")
#'}
#'
#' # Evaluation of many strategies (a little slow)
#'
#' strategies <- two_step_strategies(model)[1:5,]
#'
#' result <-  strategy_evaluation(model,
#'   										strategies  = strategies,
#'   										query = "(Y[S=1] != Y[S=0])",
#'   										given = "Y==0 & X==1",
#'   										prices = c(1, 1.5, 1.2, .5, .7))
#' result

strategy_evaluation <- function(model,
																strategies,
																query,
																given = NULL,
																prices = NULL) {

	prior_mean <- query_model(model = model, queries = query, subsets = given)$mean
	prior_variance <- prior_mean*(1-prior_mean)
  if(is.na(prior_mean)) stop("Prior not defined. Check for impossible conditions.")

	prior = c(expected_n = 0, expected_variance = prior_variance, expected_cost = 0)

	vars <- model$nodes
	if(is.null(prices)) prices <- rep(1, length(vars))

	if(!is.null(given)) {
		given_vars <- stringr::str_extract_all(given, boundary("word"))[[1]]
		given_vars <- given_vars[(given_vars %in% vars)]
	} else {given_vars <- NULL}

	if(is.vector(strategies)) if(length(strategies) !=3){stop("A single strategy should be of length 3")
		} else {strategies <- matrix(strategies, 1, 3)}

	x <- apply(strategies, 1, function(j) {
	 gbiqqtools:::strategy_evaluation_single(model,
	   										strategy  = j,
	   										query = query,
	   										given = given,
	   										prices = prices,
	   										vars = vars,
	   										given_vars = given_vars)}) %>%
		  t() %>%
			data.frame(stringsAsFactors = FALSE)

	# Export
	labels <- c("priors", apply(strategies, 1, paste, collapse = "-"))

	cbind(labels,
				rbind(prior, x))
	}


#' Internal function for evaluating a single strategy

strategy_evaluation_single <- function(model,
																strategy,
																query,
																given = NULL,
																prices = NULL,
																vars = model$nodes,
																given_vars = NULL  #Names of vars in "given"
																) {

	# Strategy 1

	s1 <- strategy[[1]]

	if(is.na(s1)) stop("First strategy should not be NA")

	na_vars <- vars[!(vars %in% unlist(c(s1, given_vars)))]

	givens  <- paste(given, paste0("is.na(", na_vars, ")", collapse = " & "), sep = " & ")

	# Conditional inferences
	ci  <- conditional_inferences(model, query = query, given = givens)

	# Step 1 conclusions depending on findings on strategy 1
	p0 <- filter(ci, ci[s1]==0)[1, c("posterior", "prob")]
	p1 <- filter(ci, ci[s1]==1)[1, c("posterior", "prob")]
	p0[is.na(p0)] <- 0  # NAs arise if no possible events
	p1[is.na(p1)] <- 0  # NAs arise if no possible events

	# Probability of finding 0/1 in step 1 (normalized)
	probs  <- c(p0$prob, p1$prob)/(p0$prob + p1$prob)

	# Gather posteriors
	posts <- c(p0$posterior, p1$posterior)

	evs <-
		sapply(0:1, function(j){
			s <- strategy[[2+j]]
			ifelse(is.na(s),
						 (posts[1+j])*(1-posts[1+j]),
						 ifelse(probs[1+j] == 0, -1,  # - 1 is placeholder for posteriors on zero probability event
						 expected_learning(model,
						 									query = query,
						 									strategy = paste(s),
						 									given = paste(given, "&", strategy[[1]], "==", j),
						 									parameters = NULL)$E_post_var)
						 )
		})

	# Flag: Fill in NAs if probs = 0
	evs[probs == 0]  <- 0
	strategy[2:3][probs == 0]  <- "NONE"

	# Expected Number of clues sought: 1 plus 1 if additional steps sought
	expected_n        <- 1 + probs%*%as.vector((!is.na(strategy[2:3])))

	# Expected variance given conditional strategies
	expected_variance <- probs%*%evs

	# Costs
	costs <- c(sum(prices[vars %in% factor(strategy[-3])]),
						 sum(prices[vars %in% factor(strategy[-2])]))
	expected_costs <- probs%*%costs

	c(expected_n =        expected_n,
		expected_variance = expected_variance,
		expected_costs =    expected_costs)
}



#' Make set of two step strategies
#'
#' @param model A causal model made by make_model
#' @param vars Variables to
#' @export
#' @examples
#' two_step_strategies(make_model("X->M->Y))
two_step_strategies <- function(model, vars = model$nodes){

	x <- sapply(1:length(vars), function(j)
		expand.grid(one = vars[j],
								two_0 = c(NA, vars[-j]),
								two_1 = c(NA, vars[-j]),
								stringsAsFactors = FALSE),
		simplify = FALSE)

	do.call("rbind", x) %>%
      	data.frame(stringsAsFactors = FALSE)
}
