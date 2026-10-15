plot_input <- function(plotmode) {
  list(
    plotmode = plotmode,
    scheme = "auto",
    geom = "line",
    alpha = 0.8,
    lwd = 0.5,
    psize = 2,
    xlab = "",
    ylab = "",
    xmin = NA_real_,
    xmax = NA_real_,
    ymin = NA_real_,
    ymax = NA_real_,
    w = 25,
    basesize = 12,
    legcols = 0,
    title = ""
  )
}

test_that("renderer builds time-series plots with summary envelopes", {
  times <- as.POSIXct("2020-01-01", tz = "UTC") + 0:2
  data <- data.frame(
    datetime = times,
    date = as.Date(times),
    value = c(1, 2, 3),
    variable = "temperature",
    unit = "C",
    series = "temperature",
    ckey = "temperature",
    .plot_group = "series|gap|summary",
    .is_summary = TRUE,
    .bin_start = times - 0.5,
    .bin_end = times + 0.5,
    .min = c(0, 1, 2),
    .max = c(2, 3, 4)
  )

  plot <- env_plot_build_plot(data, plot_input("ts"))

  expect_s3_class(plot, "ggplot")
  expect_equal(length(ggplot2::ggplot_build(plot)$data), 2L)
})

test_that("renderer builds depth and dual-isotope plots", {
  depth <- data.frame(
    value = c(1, 2, 3),
    depth_bottom = c(10, 20, 30),
    variable = "temperature",
    unit = "C",
    series = "logger-1",
    ckey = "logger-1",
    date = as.Date("2020-01-01"),
    dataset = "sample",
    logger = "logger-1"
  )
  depth_input <- plot_input("depth")
  depth_input$facet_by <- ""
  expect_s3_class(
    env_plot_build_plot(depth, depth_input),
    "ggplot"
  )

  dual <- data.frame(
    dataset = "sample",
    series = "core-1",
    ckey = "core-1",
    date = as.Date("2020-01-01"),
    variable = rep(c("d18O", "d2H"), each = 3),
    value = c(-10, -11, -12, -70, -71, -72)
  )
  dual_input <- plot_input("dual")
  dual_input$xvar <- "d18O"
  dual_input$yvar <- "d2H"
  dual_input$fitline <- "none"
  dual_input$gmwl <- FALSE
  dual_input$lmwl <- FALSE
  expect_s3_class(
    env_plot_build_plot(dual, dual_input),
    "ggplot"
  )
})
