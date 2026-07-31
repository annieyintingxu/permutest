plot_global_p <- function(data, group_col, outcome_cols,
                          param_values, global_p_values,
                          threshold, figure_name, folder_name,
                          pop_diff = NULL,
                          save = TRUE){
  
  # Prepare the global p-values for plotting
  df <- data.frame(global_p_values)
  r <- threshold
 
  print(
    ggplot() +
      geom_tile(data = subset(df, global_p_value >= r), aes(param1, param2, fill = 
                                                              global_p_value)) +
      scale_fill_viridis_c(name = paste("p \u2265", r),
                           guide = guide_colourbar(order = 2)) +
      new_scale_fill() +
      geom_tile(data = subset(df, global_p_value < r), aes(param1, param2, 
                                                           fill = global_p_value)) +
      scale_fill_gradient(name = paste("p <", r), low = "grey90", high = "grey40",
                          guide = guide_colourbar(order = 3)) +
      labs(title = paste(figure_name),
           x = "Parameter 1",
           y = "Parameter 2") 
  )
  
  if(save){
    file_name <- paste0(folder_name, "/", figure_name, ".png")
    ggsave(file_name, width = 6, height = 5, dpi = 300)
  }
}


subtract_null <- function(theta_0, group1, group2) {
  group1_subtract <- array(data = NA, dim = dim(group1))
  group1_subtract[,1] <- group1[,1] - as.numeric(theta_0[1])
  group1_subtract[,2] <- group1[,2] - as.numeric(theta_0[2])
  return(rbind(group1_subtract, group2))  # group 2 not shifted
}


peek <- function(param_value, group1, group2, n, m, pop_diff,
                 comb_funcs, groups_permuted,
                 B,
                 data_type = NULL, r = NULL, s = NULL){
  
  # Get the number of combining functions we want to use
  num_func <- length(comb_funcs)
  
  # Shift group 1
  subtracted_data <- subtract_null(theta_0 = param_value, group1 = group1, group2 = group2)
  
  # Initialize for index and data storage 
  group1_indices_permuted <- array(data = NA, dim = c(B, m))
  data_permuted_1 <- array(data = NA, dim = c(B, m, 2))
  data_permuted_2 <- array(data = NA, dim = c(B, n-m, 2))
  
  
  # Select the permuted data
  for(b in 1:B){ # For the b_th permutation
    
    # Some acrobatics to get from groups_permuted to indices_permuted
    group1_indices <- data.frame(cbind(1:n, groups_permuted[b,])) |> 
      dplyr::filter(X2 == 1) |>
      dplyr::select(X1) |> 
      unlist()
    
    # Save the indices
    group1_indices_permuted[b,] <- group1_indices
    
    # Subset the data
    data_permuted_1[b,,] <- subtracted_data[group1_indices,]
    data_permuted_2[b,,] <- subtracted_data[-group1_indices,]
  }
  
  ## TEST STATISTIC ##
  # Original difference
  diff_orig_1 <- mean(subtracted_data[1:m, 1]) - mean(subtracted_data[(m+1):n, 1]) 
  diff_orig_2 <- mean(subtracted_data[1:m, 2]) - mean(subtracted_data[(m+1):n, 2]) 
  
  # Difference for the permuted datasets. 
  # We do row-wise mean because each row is one dataset.
  means_group1_1 <- apply(X = data_permuted_1[,,1], MARGIN = 1, FUN = mean) 
  means_group1_2 <- apply(X = data_permuted_1[,,2], MARGIN = 1, FUN = mean)
  means_group2_1 <- apply(X = data_permuted_2[,,1], MARGIN = 1, FUN = mean) 
  means_group2_2 <- apply(X = data_permuted_2[,,2], MARGIN = 1, FUN = mean)
  
  # TODO: check this part...why are the differences centered around 0
  # with most of the differences between -5 and 5 for all of the nulls?
  diff_perm_1 <- means_group1_1 - means_group2_1
  diff_perm_2 <- means_group1_2 - means_group2_2
  
  # Append the original test statistic to the permuted data's test statistics
  # NOTE: the original data will be the last entry!
  diff_perm_and_orig_1 <- c(diff_perm_1, diff_orig_1)
  diff_perm_and_orig_2 <- c(diff_perm_2, diff_orig_2)
  
  # And add the original indices as the last row
  indices_permuted_and_orig <- rbind(group1_indices_permuted, 1:m)
  # colnames(indices_permuted_and_orig) <- c('data_index1', 'data_index2', 'data_index3')
  
  # for storage
  p_values <- array(data = NA, dim = c(1 + B, 2))
  combined_perm <- array(data = NA, dim = c(B+1, 1)) # Just Tippett for now.
  
  ## P-values ##
  for(b in 1:(B+1)){
    
    # Calculate the p-values for the bth permutation
    p_value_perm_1 <- sum(diff_perm_and_orig_1 >= 
                            diff_perm_and_orig_1[b]) / (1 + B)
    
    p_value_perm_2 <- sum(diff_perm_and_orig_2 >= 
                            diff_perm_and_orig_2[b]) / (1 + B)
    
    p_values[b,] <- c(p_value_perm_1, p_value_perm_2)
    
    if("tippett" %in% comb_funcs){
      
      # Tippett's combining function
      combined_perm[b,1] <- max(1-p_value_perm_1, 1-p_value_perm_2)
    } else{
      
      return("Tippett's was not in the comb_funcs argument. Check your spelling
          or adjust the code to return outputs from multiple combination
          functions.")
      
      # # Product combining function
      # combined_perm[b, 1] <- -2*sum(log(p_value_perm_1), log(p_value_perm_2))
      # 
      # # Liptak's combining function
      # combined_perm[b, 2] <- sum(qnorm(1 - p_value_perm_1), qnorm(1 - p_value_perm_2))
      # 
      # # Edgington's combining function
      # combined_perm[b, 4] <- sum(1-p_value_perm_1, 1-p_value_perm_2)
    }
  }
  
  # Peek: component differences in means, component p-values, combined p-values
  df_components <- cbind(c(1:B, 0), indices_permuted_and_orig, 
                         diff_perm_and_orig_1, diff_perm_and_orig_2, 
                         p_values, combined_perm)
  
  data_index_names <- paste0('data_index', c(1:m))
  
  colnames(df_components) <- c('draw_index', 
                               data_index_names,
                               'diff_means_1', 'diff_means_2',
                               'p_value_1', 'p_value_2',
                               'combined')
  
  ## Global p-values ##
  global_p_values <- rep(NA, num_func)
  
  for(j in 1:num_func){
    combined_orig <- combined_perm[B+1,j]
    
    # TODO: should this be >= or < if the individual p-values are calculated with >= ?
    global_p_values[j] <- sum(combined_perm[,j] >= combined_orig)/(1 + B)
  }
  
  return(list(global_p_values = global_p_values,
              df_components = df_components))
}
