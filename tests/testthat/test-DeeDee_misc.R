
test_that("summary report", {
  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  temp_output_file <- "test_report.html"
  deedee_summary(
    dde = dde,
    template = system.file("extdata",
                           "summary_template.Rmd",
                           package = "DeeDee"),
    output_path = temp_output_file
  )

  # trigger error by trying to overwrite
  expect_error(
    deedee_summary(
      dde = dde,
      template = system.file("extdata",
                             "summary_template.Rmd",
                             package = "DeeDee"),
      output_path = temp_output_file
    )
  )

  # browseURL(temp_output_file)
  expect_true(file.exists(temp_output_file))

  file.remove(temp_output_file)
})


