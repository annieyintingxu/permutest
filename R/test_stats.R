#' Calculate difference in means
#'
#' This function takes a data frame, and group and outcome column names as input
#' and returns the difference in mean outcome between the two groups
#'
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_col The name of the column in df that corresponds to the outcome variable
#' @param treatment_value The value of group_col to be considered 'treatment'
#' @return The difference in mean outcome between the two groups
#' @examples
#' data <- data.frame(group = c(rep(1, 4), rep(2, 4)),
#'                    outcome = c(rep(3, 4), rep(5, 4)))
#'
#' diff_in_means(df = data,
#'               group_col = "group",
#'               outcome_col = "outcome",
#'               treatment_value = 1)
#'
#' @export
diff_in_means <- function(df, group_col, outcome_col, treatment_value=NULL){
  # get unique groups
  groups <- unique(df[[group_col]])
  # if more or less than 2 groups, throw error
  if(length(groups) != 2){
    stop("Error: dataset must contain exactly 2 unique groups to use this test statistic")
  }
  # check which group to consider 'treatment'
  if(is.null(treatment_value)){
    warning(paste0("No value for parameter treatment_value, using ", groups[1], " as treatment group"))
    treatment_group = groups[1]
    control_group = groups[2]
  } else if (treatment_value == groups[1]) {
    treatment_group = groups[1]
    control_group = groups[2]
  } else if(treatment_value == groups[2]) {
    treatment_group = groups[2]
    control_group = groups[1]
  } else {
    stop("treatment_value is not a valid value contained in the group column")
  }
  # calculate difference in means test statistic
  t <- mean(df[[outcome_col]][df[group_col] == treatment_group]) - mean(df[[outcome_col]][df[group_col] == control_group])
  # return test statistic value
  return(t)
}


#' Calculate difference in medians
#'
#' This function takes a data frame, and group and outcome column names as input
#' and returns the difference in median outcome between the two groups
#'
#' @importFrom stats median
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_col The name of the column in df that corresponds to the outcome variable
#' @param treatment_value The value of group_col to be considered 'treatment'
#' @return The difference in median outcome between the two groups
#' @examples
#' data <- data.frame(group = c(rep(1, 4), rep(2, 4)),
#'                    outcome = c(rep(3, 4), rep(5, 4)))
#'
#' diff_in_medians(df = data,
#'               group_col = "group",
#'               outcome_col = "outcome",
#'               treatment_value = 1)
#' @export
diff_in_medians <- function(df, group_col, outcome_col, treatment_value=NULL){
  # get unique groups
  groups <- unique(df[[group_col]])
  # if more or less than 2 groups, throw error
  if(length(groups) != 2){
    stop("Error: dataset must contain exactly 2 unique groups to use this test statistic")
  }
  # check which group to consider 'treatment'
  if(is.null(treatment_value)){
    warning(paste0("No value for parameter treatment_value, using ", groups[1], " as treatment group"))
    treatment_group = groups[1]
    control_group = groups[2]
  } else if (treatment_value == groups[1]) {
    treatment_group = groups[1]
    control_group = groups[2]
  } else if(treatment_value == groups[2]) {
    treatment_group = groups[2]
    control_group = groups[1]
  } else {
    stop("treatment_value is not a valid value contained in the group column")
  }
  # calculate difference in means test statistic
  t <- median(df[[outcome_col]][df[group_col] == treatment_group]) - median(df[[outcome_col]][df[group_col] == control_group])
  # return test statistic value
  return(t)
}


#' Calculate t-test statistic
#'
#' This function takes a data frame, and group and outcome column names as input
#' and returns the t test statistic
#'
#' @importFrom stats t.test
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_col The name of the column in df that corresponds to the outcome variable
#' @return The t test statistic
#' @export
ttest_stat <- function(df, group_col, outcome_col){
  # get unique groups
  groups <- unique(df[[group_col]])
  # if more or less than 2 groups, throw error
  if(length(groups) != 2){
    stop("Error: dataset must contain exactly 2 unique groups to use this test statistic")
  }
  # calculate the t-test statistic
  t <- t.test(df[[outcome_col]][df[group_col] == groups[1]], df[[outcome_col]][df[group_col] == groups[2]])$statistic
  t <- as.numeric(t)
  # return test statistic value
  return(t)
}

#' Calculate one-way anova test statistic
#'
#' This function takes a data frame, and group and outcome column names as input
#' and returns the one-way anova test statistic
#'
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_col The name of the column in df that corresponds to the outcome variable
#' @return The one-way anova test statistic:
#' \eqn{\sum_{g=1}^G n_g(\overline{X_g} - \overline{X})^2} where \eqn{g} indexes the groups
#' @export
one_way_anova_stat <- function(df, group_col, outcome_col){
  # get unique groups
  groups <- unique(df[[group_col]])
  # calculate overall mean
  overall_mean <- mean(df[[outcome_col]])
  # initialize t
  t <- 0
  for(g in groups){
    x <- df[[outcome_col]][df[group_col] == g]
    n <- length(x)
    t <- t + n * (mean(x) - overall_mean)^2
  }

  return(t)
}

#' Calculate the one-sample problem test statistic
#'
#' This function takes a data frame, and group and outcome column names as input
#' and returns the mean of the product of the outcome and group. This test statistic
#' is used for the one-sample problem.
#'
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_col The name of the column in df that corresponds to the outcome variable
#' @return The one-sample problem test statistic: the mean of the product of the outcome and group
#' @examples
#' data <- data.frame(group = c(rep(1, 4), rep(2, 4)),
#'                    outcome = c(rep(3, 4), rep(5, 4)))
#'
#' one_sample_mean(df = data,
#'               group_col = "group",
#'               outcome_col = "outcome")
#'
#' @export
one_sample_mean <- function(df, group_col, outcome_col){
  t <- mean(df[[outcome_col]] * df[[group_col]])

  return(t)
}


