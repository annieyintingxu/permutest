# Functions documented with @export must actually be exported in NAMESPACE.
# (devtools::load_all() exposes unexported functions, so other tests would not
# notice a missing export; this reads the namespace's export list directly.)
test_that("npc_many is exported", {
  expect_true("npc_many" %in% getNamespaceExports("permutest"))
})
