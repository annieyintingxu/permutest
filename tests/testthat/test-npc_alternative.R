# Partial p-values for the relabellings must follow `alternative`, like the
# observed partial p-value does. Checked with symmetries that hold exactly for a
# fixed set of relabellings: negating the outcomes (and the shift) swaps
# "greater" and "less", and leaves "two-sided" unchanged.

# 3 treated, 4 control; every one of the 35 assignments used once.
# Unbalanced on purpose: with equal group sizes and complete enumeration the
# permutation distribution is exactly symmetric, which hides the ranking bug.
make_case <- function() {
  df <- data.frame(group = c(1, 1, 1, 2, 2, 2, 2),
                   y1 = c(3.1, 0.4, 2.2, 1.0, -0.7, 0.9, 2.6),
                   y2 = c(-1.2, 2.5, 0.3, 0.8, 1.9, -0.4, 0.1))
  perm_set <- t(utils::combn(7, 3, function(i) { g <- rep(2, 7); g[i] <- 1; g }))
  list(df = df, perm_set = perm_set)
}

npc_p <- function(df, alternative, shift, combn) {
  suppressWarnings(
    npc(df, "group", c("y1", "y2"), alternative = alternative, shift = shift,
        combn = combn, perm_set = make_case()$perm_set)$omnibus_p)
}

test_that("npc: 'less' on y equals 'greater' on -y", {
  cs <- make_case(); neg <- cs$df; neg$y1 <- -neg$y1; neg$y2 <- -neg$y2
  for (cb in c("tippett", "fisher", "liptak")) {
    for (s in list(c(0, 0), c(0.5, -1))) {
      expect_equal(npc_p(cs$df, "less", s, cb), npc_p(neg, "greater", -s, cb),
                   info = paste(cb, paste(s, collapse = ",")))
    }
  }
})

test_that("npc: 'two-sided' is unchanged by negating the outcomes", {
  cs <- make_case(); neg <- cs$df; neg$y1 <- -neg$y1; neg$y2 <- -neg$y2
  for (cb in c("tippett", "fisher", "liptak")) {
    for (s in list(c(0, 0), c(0.5, -1))) {
      expect_equal(npc_p(cs$df, "two-sided", s, cb), npc_p(neg, "two-sided", -s, cb),
                   info = paste(cb, paste(s, collapse = ",")))
    }
  }
})
