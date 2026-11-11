# npc_many(): defaults match npc(), and every grid point equals a direct npc()
# call with the same relabellings, across alternatives, test statistics and
# combining functions.

many_case <- function() {
  df <- data.frame(group = c(1, 1, 1, 2, 2, 2, 2),
                   y1 = c(3.1, 0.4, 2.2, 1.0, -0.7, 0.9, 2.6),
                   y2 = c(-1.2, 2.5, 0.3, 0.8, 1.9, -0.4, 0.1))
  perm_set <- t(utils::combn(7, 3, function(i) { g <- rep(2, 7); g[i] <- 1; g }))
  grid <- expand.grid(eta1 = c(-1, 0.5), eta2 = c(-0.5, 1))
  list(df = df, perm_set = perm_set, grid = grid)
}

test_that("npc_many: defaults are greater, diff_in_means and fisher", {
  cs <- many_case()
  run <- function(...) suppressWarnings(
    npc_many(cs$df, "group", c("y1", "y2"), param_values = cs$grid,
             perm_set = cs$perm_set, ...))$global_p_values[, "global_p_value"]
  expect_equal(run(), run(alternative = "greater", test_stat = "diff_in_means",
                          combn = "fisher"))
})

test_that("npc_many matches npc() at every grid point", {
  cs <- many_case()
  for (alt in c("greater", "less", "two-sided")) {
    for (ts in c("diff_in_means", "diff_in_medians")) {
      for (cb in c("fisher", "tippett", "liptak")) {
        many <- suppressWarnings(
          npc_many(cs$df, "group", c("y1", "y2"), alternative = alt,
                   param_values = cs$grid, test_stat = ts, combn = cb,
                   perm_set = cs$perm_set))$global_p_values[, "global_p_value"]
        direct <- sapply(seq_len(nrow(cs$grid)), function(l) suppressWarnings(
          npc(cs$df, "group", c("y1", "y2"), alternative = alt,
              shift = as.numeric(cs$grid[l, ]), test_stat = ts, combn = cb,
              perm_set = cs$perm_set))$omnibus_p)
        expect_equal(unname(many), direct, info = paste(alt, ts, cb))
      }
    }
  }
})
