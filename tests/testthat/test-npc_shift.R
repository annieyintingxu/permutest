# npc() should accept a single shift value for several outcomes by repeating it
# (with a message), as requested in review. A shift vector of the wrong length
# is an error.

shift_case <- function() {
  df <- data.frame(group = c(1, 1, 1, 2, 2, 2, 2),
                   y1 = c(3.1, 0.4, 2.2, 1.0, -0.7, 0.9, 2.6),
                   y2 = c(-1.2, 2.5, 0.3, 0.8, 1.9, -0.4, 0.1))
  perm_set <- t(utils::combn(7, 3, function(i) { g <- rep(2, 7); g[i] <- 1; g }))
  list(df = df, perm_set = perm_set)
}

run_npc <- function(shift, ...) {
  cs <- shift_case()
  suppressWarnings(npc(cs$df, "group", c("y1", "y2"), alternative = "greater",
                       shift = shift, combn = "tippett", perm_set = cs$perm_set, ...))
}

test_that("npc: the default shift works with two outcomes", {
  cs <- shift_case()
  out <- suppressMessages(suppressWarnings(
    npc(cs$df, "group", c("y1", "y2"), perm_set = cs$perm_set)))
  ref <- suppressWarnings(
    npc(cs$df, "group", c("y1", "y2"), shift = c(0, 0), perm_set = cs$perm_set))
  expect_equal(out$omnibus_p, ref$omnibus_p)
})

test_that("npc: a single shift is repeated for every outcome, with a message", {
  expect_message(run_npc(0.5), "repeated")
  expect_equal(suppressMessages(run_npc(0.5))$omnibus_p, run_npc(c(0.5, 0.5))$omnibus_p)
})

test_that("npc: a full-length shift vector gives no message", {
  expect_no_message(run_npc(c(0.5, -1)))
})

test_that("npc: a shift vector of the wrong length is an error", {
  expect_error(run_npc(c(0.5, -1, 2)), "length")
})
