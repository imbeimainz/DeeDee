test_that("creating", {
  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  print(dde)

  expect_is(dde, "DeeDeeExperiment")

  dde_only_de <- DeeDeeExperiment(
    de_results = de_named_list
  )
  expect_is(dde, "DeeDeeExperiment")

  dde_nodd <- DeeDeeExperiment(
    se = se_macrophage_noassays,
  )
  expect_is(dde_nodd, "DeeDeeExperiment")


  expect_is(
    get_dea_df(dde, "ifng_vs_naive"), "DataFrame"
  )

  expect_error(
    get_dea_df(dde, "wrong_name")
  )

  dde_gone_wrong <- dde
  rowData(dde_gone_wrong)[["ifng_vs_naive_log2FoldChange"]] <- NULL
  expect_error(
    get_dea_df(dde_gone_wrong, "ifng_vs_naive")
  )

  expect_error(
    DeeDeeExperiment(
      rowData(se_macrophage_noassays),
      de_results = de_named_list
    )
  )

  expect_error(
    DeeDeeExperiment(
      assay(se_macrophage),
      de_results = de_named_list
    )
  )

  expect_error(
    DeeDeeExperiment()
  )
})


test_that("adding and removing", {
  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  new_del <- list(
    ifng2 = de_named_list$ifng_vs_naive,
    ifngsalmo2 = de_named_list$ifngsalmo_vs_naive
  )
  # add a new (set of) DE result(s)
  dde_new <- add_dea(dde, new_del)
  expect_is(dde_new, "DeeDeeExperiment")
  expect_equal(length(dea(dde)), 4)
  expect_equal(length(dea(dde_new)), 6)

  dde_removed <- remove_dea(dde, "ifngsalmo_vs_naive")
  expect_is(dde_removed, "DeeDeeExperiment")
  expect_equal(length(dea(dde_removed)), 3)
})


test_that("validity and so", {
  dde2 <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  expect_true(validObject(dde2))

  rowData(dde2)[["ifng_vs_naive_log2FoldChange"]] <- NULL
  expect_error(validObject(dde2))

  dde3 <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  # invalid replacements
  ## actually, can enable it is an empty list (i.e. no DE (yet) inserted)
  # dde3@dea <- list()
  # expect_error(validObject(dde3))

  dde3@dea <- list("foo", "bar")
  expect_error(validObject(dde3))
})
