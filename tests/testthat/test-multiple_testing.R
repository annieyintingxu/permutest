test_that("fisher works", {
  output <- fisher(seq(0.05, 0.9, length.out = 5))
  expect_equal(round(output, 2), 11.12)
})

test_that("liptak works", {
  output <- liptak(seq(0.05, 0.9, length.out = 5))
  expect_equal(round(output, 2), 0.57)
})

test_that("tippett works", {
  output <- tippett(seq(0.05, 0.9, length.out = 5))
  expect_equal(round(output, 2), 0.95)
})

# # Annie to-do: adjust this to test the updated npc function
# test_that("npc works", {
#   data <- data.frame(group = c(rep(1, 4), rep(2, 4)),
#                     out1 = c(0, 1, 0, 0, 1, 1, 1, 0),
#                     out2 = rep(1, 8))
# 
#   # test stat for out1 is smaller if X all 0s which about 0.015 chance so p-value should be about 0.985
#   output <- npc(df = data, group_col = "group", outcome_cols = c("out1", "out2"),
#                 alternative = 'greater', shift = c(0, 0),
#                 reps = 10^4, seed=42)
#   expect_equal(round(output, 3), 0.986)
# 
#   # all permutations result in same test stat, 
#   # and alternative is 'greater', so combined p-value should be 1
#   data <- data.frame(group = c(1, 1, 2, 2),
#                      out1 = rep(0, 4),
#                      out2 = rep(1, 4))
#   output <- npc(df = data, group_col = "group", outcome_cols = c("out1", "out2"), 
#                 alternative = 'greater',
#                 shift = c(2, 5), reps = 10^4)
#   expect_equal(round(output, 3), 1)
#   
#   # all permutations result in same test stat, so combined p-value should be small
#   data <- data.frame(group = c(1, 1, 2, 2),
#                      out1 = rep(0, 4),
#                      out2 = rep(1, 4))
#   output <- npc(df = data, group_col = "group", outcome_cols = c("out1", "out2"), 
#                 alternative = 'less',
#                 shift = c(1, 1), reps = 10^4)
#   expect_equal(round(output, 3), 1)
# 
#   # 4/24 permutations result in test stat of same size
#   data <- data.frame(group = c(1, 1, 2, 2),
#                      out1 = c(2, 2, 1, 1),
#                      out2 = c(2, 2, 1, 1))
#   output <- npc(df = data, group_col = "group", outcome_cols = c("out1", "out2"), 
#                 reps = 10^4, seed = 42)
#   expect_equal(round(output, 2), .17)
# })

# Annie to-do: adjust this to test the updated npc function
test_that("npc works", {
  data <- data.frame(group = c(rep(1, 4), rep(2, 4)),
                     out1 = c(0, 1, 0, 0, 1, 1, 1, 0),
                     out2 = rep(1, 8))
  
  # test stat for out1 is smaller if X all 0s which about 0.015 chance so p-value should be about 0.985
  output <- npc(df = data, group_col = "group", outcome_cols = c("out1", "out2"),
                alternative = 'greater', shift = c(0, 0),
                reps = 10^4, seed=42)
  expect_equal(round(output$omnibus_p, 3), 0.986)
  
  # all permutations should give larger test stat than the observed,
  # and alternative is 'greater', so combined p-value should be 1
  data <- data.frame(group = c(1, 1, 2, 2),
                     out1 = rep(0, 4),
                     out2 = rep(1, 4))
  output <- npc(df = data, group_col = "group", outcome_cols = c("out1", "out2"), 
                alternative = 'greater',
                shift = c(2, 5), reps = 10^4)
  expect_equal(round(output$omnibus_p, 3), 1)
  
  # all permutations should give larger test stat than the observed,
  # and alternative is 'less', so combined p-value should be 1/(reps + 1)
  # unless one or more of the randomly drawn permutations
  # happens to be the same as the observed permutation.
  # For seed = 92748, the p-value should be 0.1642836
  data <- data.frame(group = c(1, 1, 2, 2),
                     out1 = rep(0, 4),
                     out2 = rep(1, 4))
  output <- npc(df = data, group_col = "group", outcome_cols = c("out1", "out2"), 
                alternative = 'less',
                shift = c(1, 1), reps = 10^4,
                seed = 92748)
  expect_equal(round(output$omnibus_p, 3), 0.164)
  
  # 4/24 permutations result in test stat of same size
  data <- data.frame(group = c(1, 1, 2, 2),
                     out1 = c(2, 2, 1, 1),
                     out2 = c(2, 2, 1, 1))
  output <- npc(df = data, group_col = "group", outcome_cols = c("out1", "out2"), 
                alternative = 'greater', shift = c(0, 0),
                reps = 10^4, seed = 42)
  expect_equal(round(output$omnibus_p, 2), .17)
})

test_that("npc_many works", {
  ##### Test 1 for npc_many ######
  data <- data.frame(group = c(1, 1, 2, 2), 
                     out1 = c(1, 6, 50, 11),
                     out2 = c(-90, -30, -60, -15))
  
  param_array <- data.frame(param1 = c(-100, -26, 0),
                            param2 = c(-100, -26, 10))
  
  output1 <- npc_many(df = data,
                     group_col = 'group',
                     outcome_cols = c('out1', 'out2'),
                     alternative = 'greater',
                     param_values = param_array, 
                     test_stat = 'diff_in_means',
                     combn = 'tippett',
                     reps = 10^4,
                     perm_set = NULL,
                     seed = NULL)
  
  # Even though our dataset is very small, param=(-100, -100) is way 
  # to the left and down compared to 
  # the sample difference, so the p-value should be small.
  # The output should slightly differ every time because 
  # we randomly draw the permutations.
  # Previous runs suggest that the value should be between 0.15 and 0.18.
  output1_1 <- as.numeric(output1$global_p_values[1,'global_p_value'])
  expect_gt(output1_1, 0.150) 
  expect_lt(output1_1, 0.180) 
  
  # Param=(-26, -26) is quite close to the sample difference,
  # so the p-value should be not too small nor too large.
  # Previous runs suggest that it ranges between 0.5 and 0.9.
  output1_2 <- as.numeric(output1$global_p_values[2,'global_p_value'])
  expect_gt(output1_2, 0.5) 
  expect_lt(output1_2, 0.9) 
  
  # Param=(0, 10) is a fair bit to the right and up compared to the sample difference,
  # so the p-value should be close to 1. 
  # Previous runs suggest that it is about 1.000.
  output1_3 <- as.numeric(output1$global_p_values[3,'global_p_value'])
  expect_gt(output1_3, 0.999) 
  expect_lte(output1_3, 1) 
  
  
  ##### Test 2 for npc_many ######
  data <- data.frame(group = c(rep(1, 20), rep(2, 10)), 
                     out1 = rep(1, 30),
                     out2 = rep(200, 30))
  
  param_array <- data.frame(param1 = c(-1, 0, 1),
                            param2 = c(-1, 0, 2))
  
  output2 <- npc_many(df = data,
                     group_col = 'group',
                     outcome_cols = c('out1', 'out2'),
                     alternative = 'greater',
                     param_values = param_array, 
                     test_stat = 'diff_in_means',
                     combn = 'tippett',
                     reps = 10^4,
                     perm_set = NULL,
                     seed = NULL)
  
  # With param=(-1, -1), 
  # all of the permutations should have a test statistic
  # less than the observed, unless some of the random permutations
  # are the same exact permutation as the observed.
  # The p-value should be between 0 and 0.001.
  output2_1 <- as.numeric(output2$global_p_values[1,'global_p_value'])
  expect_gt(output2_1, 0) 
  expect_lt(output2_1, 10^-3) 
  
  # At param=(1, 2), 
  # all of the permutations should have test statistics
  # greater than the observed, unless some of the permutations
  # are the exact same permutation as the observed.
  # The p-value should be between 0.999 and 1.
  output2_2 <- as.numeric(output2$global_p_values[2,'global_p_value'])
  expect_gt(output2_2, 0.999) 
  expect_lte(output2_2, 1) 
  
  # Similarly, at param=(1, 2), 
  # all of the permutations should have test statistics
  # greater than the observed, unless some of the permutations
  # are the exact same permutation as the observed.
  # The p-value should be between 0.999 and 1.
  output2_3 <- as.numeric(output2$global_p_values[3,'global_p_value'])
  expect_gt(output2_3, 0.999) 
  expect_lte(output2_3, 1)
  
  # TODO: add at least one more test and pass in a seed
  # and use expect_equal() instead of inequalities
})

test_that("p-value adjustment works", {
  output <- adjust_p_value(pvalues = c(.01, .04, .03, .005), method = 'holm-bonferroni')
  expect_equal(output, c(.03, .06, .06, .02))

  output <- adjust_p_value(pvalues = c(.01, .04, .03, .005), method = 'bonferroni')
  expect_equal(output, 4*c(.01, .04, .03, .005))

  output <- adjust_p_value(pvalues = c(0.01, 0.001, 0.05, 0.20, 0.15, 0.15), method = 'benjamini-hochberg')
  expect_equal(output, c(0.030, 0.006, 0.100, 0.200, 0.180, 0.180))
})
