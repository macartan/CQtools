#' Get leave one out data likelihood
#'
#' Assesses data likelihood using leave out procedure
#' 
#' @inheritParams CQtools_internal_inherit_params
#' @param data Data in long format
#' @param sims Number of simulations
#' @param ... arguments passed to update_model
#' @export
#' @examples
#' model <- make_model("X->Y")
#' data <- simulate_data(model, n = 4)
#' get_loo_likelihood(data, model)
#' 
get_loo_likelihood <- function(data, model, sims=100, ...){
    
    short_data <- collapse_data(data, model) 
    rows       <- 1:nrow(short_data)
    
    case_likelihood <-
        lapply(rows,
               function(j) {
                   if (short_data$count[j] == 0) return(1)
                   loo_data <-
                       mutate(short_data, count = ifelse(count > 0 & rows == j, count - 1, count))
                   loo_left <-
                       mutate(short_data, count = ifelse(count > 0 & rows == j, 1, 0)) %>%
                       expand_data(model)
                   loo_post <-
                       update_model(model, loo_data, data_type = "compact", ...)
                   k <-
                       min(sims, nrow(loo_post$posterior_distribution))
                   loo_prob <-
                       apply(loo_post$posterior_distribution[1:k, ], 1, function(pars)
                           get_data_probs(model, loo_left, parameters = pars)[short_data$event[j]])
                   mean(loo_prob)
               })
    
    list(data = short_data, model = model, case_likelihood = case_likelihood, 
         loo_likelihood = prod((case_likelihood %>% unlist)^(short_data$count)))
}
