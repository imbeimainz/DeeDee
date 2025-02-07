test_that("shiny", {
  # data(dd_list_original, package = "DeeDee")

  library("shiny")
  myapp <- deedee_app(dde = DeeDeeExperiment(se = se_macrophage_noassays, de_results = de_named_list))

  expect_is(myapp, "shiny.appobj")
})
