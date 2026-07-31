#' Fisher combining function
#'
#' This function takes an array of p-values and returns a combined p-value using fisher's combining function:
#' \eqn{-2 \sum_i \log(p_i)}
#'
#' @param pvalues Array of p-values
#' @return Combined p-value using fisher's method
#' @export
#' @examples
#' fisher(pvalues = c(.05, .1, .5))
#'
fisher <- function(pvalues){
  # compute fisher combination
  value <- -2 * log(prod(pvalues))

  return(value)
}


#' Tippett combining function
#'
#' This function takes an array of p-values and returns a combined p-value using Tippett's combining function:
#' \eqn{\max_i \{1-p_i\}}
#'
#' @param pvalues Array of p-values
#' @return Combined p-value using Tippett's method
#' @export
#' @examples
#' tippett(pvalues = c(.05, .1, .5))
#'
tippett <- function(pvalues){
  # compute tippett
  value <- max(1-pvalues)

  return(value)
}


#' Liptak combining function
#'
#' This function takes an array of p-values and returns a combined p-value using Liptak's combining function:
#' \eqn{\sum_i \Phi^{-1}(1-p_i)} where \eqn{\Phi} is the CDF of the Normal distribution
#'
#' @importFrom stats qnorm
#' @param pvalues Array of p-values
#' @return Combined p-value using Liptak's method
#' @export
#' @examples
#' liptak(pvalues = c(.05, .1, .5))
#'
liptak <- function(pvalues){

  value <- sum(qnorm(1-pvalues))

  return(value)
}


#' Run NPC
#'
#' This function takes a data frame and group and outcome column names as input
#' and returns the nonparametric combination of tests (NPC) omnibus p-value
#'
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_cols The names of the columns in df that corresponds to the outcome variable
#' @param alternative A string. Options are "greater", "less", or "two-sided"
#' @param shift Value of shift to apply in one- or two-sample problem. Can be vector-valued
#' @param strata_col The name of the column in df that corresponds to the strata
#' @param test_stat Test statistic function
#' @param perm_func Function to permute group, default is permute_group which randomly permutes group assignment
#' @param combn Combining function method to use, takes values 'fisher', 'tippett', or 'liptak', or a user defined function
#' @param reps Number of iterations to use when calculating permutation p-value
#' @param perm_set Matrix of permutations to use instead of reps iterations of perm_func
#' @param complete_enum Boolean, whether to calculate P-value under complete enumeration of permutations
#' @param seed An integer seed value
#' @return The omnibus p-value
#' @export
#' @examples
#' data <- data.frame(group = c(rep(1, 4), rep(2, 4)),
#' out1 = c(0, 1, 0, 0, 1, 1, 1, 0),
#' out2 = rep(1, 8))
#' output <- npc(df = data, group_col = "group", alternative = "greater",
#'               outcome_cols = c("out1", "out2"), perm_func = permute_group,
#'               combn = "tippett", shift = c(0,0), reps = 10^4, seed = 42)
#'
npc <- function(df, group_col, outcome_cols, alternative, shift, 
                strata_col = NULL,
                test_stat = "diff_in_means",
                perm_func = permute_group,
                combn = "tippett",
                reps = 10^4,
                perm_set = NULL,
                complete_enum = FALSE,
                seed = NULL){

  # Check the combination function
  if(is.character(combn) && combn == "fisher"){
    combn_func <- fisher
  } else if(is.character(combn) && combn == "tippett"){
    combn_func <- tippett
  } else if(is.character(combn) && combn == "liptak"){
    combn_func <- liptak
  } else {
    combn_func <- combn
  }
  
  # Get the number of reps if a permutation set was provided
  if(!is.null(perm_set)){
    reps = nrow(perm_set)
  }

  # Get matrix of p-values
  n_out <- length(outcome_cols)
  obs_p_value <- rep(NA, n_out)
  p_value_mat <- matrix(NA, nrow = reps, ncol = n_out)

  for(i in 1:n_out){
    
    # Access the associated shift parameter
    if(is.null(shift)){
      shift_i <- 0
    } else { 
      shift_i <- as.numeric(shift)[i] 
    }
    
      output <- permutation_test(df = df, group_col = group_col,
                                 outcome_col = outcome_cols[i], strata_col = strata_col,
                                 test_stat = test_stat, perm_func = perm_func,
                                 shift = shift_i, reps = reps, perm_set = perm_set,
                                 complete_enum = complete_enum, alternative = alternative,
                                 return_test_dist = T, return_perm_dist = T,
                                 seed = seed)

      obs_p_value[i] <- output$p_value

    if(i == 1){
      perm_set <- output$perm_indices_mat
      reps <- nrow(perm_set)
    }
      p_value_mat[,i] <- (reps - rank(output$test_stat_dist, ties.method = "min") + 1)/(reps+1)
  }
  # Get combined p-values
  combn_pvalues <- rep(NA, reps)
  obs_combn_pvalue <- combn_func(obs_p_value)

  for(j in 1:reps){
    combn_pvalues[j] <- combn_func(p_value_mat[j, ])
  }
  
  # Get omnibus p-values
  omnibus_p <- (sum(combn_pvalues >= obs_combn_pvalue)+1) / (reps+1)
  
  # Return both the omnibus p-value and the permutation set used
  return(list(omnibus_p = omnibus_p,
              perm_set = perm_set))
}

#' Adjust p-values for multiple testing
#'
#' This function takes an array of p-values and returns adjusted p-values using
#' user-inputted FWER or FDR correction method
#'
#' @param pvalues Array of p-values
#' @param method The FWER or FDR correction to use, either 'holm-bonferroni',
#' 'bonferroni', or 'benjamini-hochberg'
#' @return Adjusted p-values
#' @export
#' @examples
#' adjust_p_value(pvalues = c(.05, .1, .5), method='holm-bonferroni')
#'
adjust_p_value <- function(pvalues, method='holm-bonferroni'){
  # get number of p-values
  n <- length(pvalues)
  if(method == 'holm-bonferroni'){
    order <- rank(pvalues, ties.method = 'last')
    adj_pvalues <- pmin(pvalues * (n - order + 1), rep(1, n))
    prev_index <- which(order == 1)
    for (i in 1:n) {
      current_index <- which(order == i)
      adj_pvalues[current_index] <- max(adj_pvalues[prev_index], adj_pvalues[current_index])
      prev_index <- current_index
    }
  } else if (method == "bonferroni"){
    adj_pvalues <- pmin(pvalues*n, rep(1, n))
  } else if (method == "benjamini-hochberg"){
    order <- rank(pvalues, ties.method = 'last')
    adj_pvalues <- pmin(pvalues * (n / order), rep(1, n))
    prev_index <- which(order == n)
    for (i in n:1) {
      current_index <- which(order == i)
      adj_pvalues[current_index] <- min(adj_pvalues[prev_index], adj_pvalues[current_index])
      prev_index <- current_index
    }
  } else {
    stop("Method must be 'holm-bonferroni', 'bonferroni', or 'benjamini-hochberg'")
  }

  return(adj_pvalues)
}


#' Run NPC for an array of parameter values
#'
#' This function takes a data frame, a matrix of parameter values
#' and group and outcome column names as input
#' and returns a matrix of the nonparametric combination of tests (NPC) omnibus p-values,
#' with each p-value indexed by the corresponding parameter value.
#'
#' @param df A data frame
#' @param n The total number of observations
#' @param m The number of observations in the treatment group (Group 1)
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_cols The names of the columns in df that corresponds to the outcome variable
#' @param alternative A string. Options are "greater", "less", or "two-sided"
#' @param param_values The matrix containing the parameter values to test
#' @param test_stat Test statistic function
#' @param combn Combining function method to use, takes values 'fisher', 'tippett', or 'liptak', or a user defined function
#' @param reps Number of iterations to use when calculating permutation p-value
#' @param perm_set Matrix of permutations to use instead of reps iterations of perm_func
#' @param seed An integer seed value
#' @return A list containing global_p_values, which contains the omnibus p-values, and perm_set, the permutations used.
#' @export
#' @examples
#' TODO: write down an example!!
#' data <- NA
#' out1 = NA
#' out2 = NA
#' df_param <- NA
#' output <- npc_grid(df = data, n = 20, m = 10, 
#'                    group_col = 'group',
#'                    outcome_cols = c('out1', 'out2'),
#'                    alternative = 'greater',
#'                    param_values = df_param,
#'                    test_stat = 'diff_in_means',
#'                    combn = 'tippett',
#'                    reps = 10^4, 
#'                    perm_set = NULL,
#'                    seed = 374923084)
#'
npc_many <- function(df,
                     group_col, outcome_cols, alternative,
                     param_values,
                     test_stat, combn, 
                     reps = NULL, perm_set = NULL,
                     seed = NULL){
  
  # TODO: Add support for user choice of 
  # strata_col, complete_enum, and perm_func
  
  # Print an error message if both reps and perm_set are null
  # TODO: throw an error instead of just printing something
  if(is.null(reps) & is.null(perm_set)){
    print("The parameter 'reps' and the parameter 'perm_set' 
          are both null. Please provide the desired number of 
          permutations as 'reps' or provide the actual
          set of permutations you want to use as 'perm_set'.")
  }
  
  # Convert the parameter and permutation objects to arrays
  param_values <- as.matrix(param_values) # Converts to array.
  if (!is.null(perm_set)){
    perm_set <- as.matrix(perm_set) # Converts to array.
  }
  
  # Sets up storage for the global p-values and the permutation set.
  L <- nrow(param_values)
  global_p_values_vector <- array(data = NA, dim = L)
  
  # For each null parameter value
  for (l in seq_len(L)) {
    
    # Accesses the l^th null parameter value.
    null_param <- param_values[l,]
    
    # Calculates the global p-values!
    npc_output <- npc(
      df           = df,
      group_col    = group_col,
      outcome_cols = outcome_cols,
      test_stat    = test_stat,
      perm_func    = permute_group,
      combn        = combn,
      shift        = null_param,   
      perm_set     = perm_set,
      reps         = reps,
      seed         = seed,
      alternative  = alternative,
    )
    
    # Saves the omnibus p-value for the l^th parameter value.
    global_p_values_vector[l] <- npc_output$omnibus_p
    
    # If no permutation set was provided by the user,
    # saves the permutation set generated by npc() for the first parameter value
    # and re-uses it for subsequent parameter values.
    if(l == 1 & is.null(perm_set)){
      perm_set <- npc_output$perm_set
    }
  }
  
  # Sets up an array of the p-values, indexed by the parameter values.
  num_params <- ncol(param_values)
  global_p_values <- array(data = NA, 
                           dim = c(L, num_params + 1))
  
  colnames(global_p_values) <-  c(paste0('param', seq_len(num_params)),
                                  'global_p_value')
  
  global_p_values[,'global_p_value'] <- global_p_values_vector
  
  global_p_values[, paste0('param', seq_len(num_params))] <- param_values
  
  # Returns both the p-value array and the array of permutations used.
  # If the user provided a permutation set, the returned array should
  # be the same as the original.
  return(list(global_p_values = global_p_values,
              perm_set        = perm_set))
}
