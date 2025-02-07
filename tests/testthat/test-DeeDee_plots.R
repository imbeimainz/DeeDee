dde <- DeeDeeExperiment(
  se_macrophage_noassays,
  de_results = de_named_list
)

# needed for bad handling of namespace by ggvenn
library("dplyr")

test_that("plotting & plotting", {
  p_scatter <- deedee_scatter(dde)
  expect_is(p_scatter, "gg")


  p_heatmap <- deedee_heatmap(dde)
  expect_is(p_heatmap, "Heatmap")


  p_venn <- deedee_venn(dde)
  expect_is(p_venn, "gg")

  p_venn_up <- deedee_venn(dde, mode = "up")
  expect_is(p_venn_up, "gg")
  p_venn_down <- deedee_venn(dde, mode = "down")
  expect_is(p_venn_down, "gg")


  p_upset <- deedee_upset(dde)
  expect_is(p_upset, "gg")

  p_upset_both <- deedee_upset(dde, mode = "both")
  expect_is(p_upset_both, "gg")

  p_upset_up <- deedee_upset(dde, mode = "up")
  expect_is(p_upset_up, "gg")
  p_upset_down <- deedee_upset(dde, mode = "down")
  expect_is(p_upset_down, "gg")


  p_cat <- deedee_cat(dde)
  expect_is(p_cat, "gg")


  p_bars <- deedee_bars(dde)
  expect_is(p_bars, "gg")

  p_bars_wnumbers <- deedee_bars(dde, show_DEnumbers = TRUE)
  expect_is(p_bars_wnumbers, "gg")


  p_cat_down <- deedee_cat(dde, mode = "down")
  expect_is(p_cat_down, "gg")
  p_cat_both <- deedee_cat(dde, mode = "both")
  expect_is(p_cat_both, "gg")


  p_qq <- deedee_qq(dde)
  expect_is(p_qq, "gg")
  p_qq_m <- deedee_qqmult(dde)
  expect_is(p_qq_m, "gg")

})


