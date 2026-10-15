test_that("CSV normalization preserves duplicate timestamps and reports gaps", {
  path <- tempfile(fileext = ".csv")
  writeLines(c(
    "datetime,value,logger,unit",
    "2020-01-01 00:00:00,1,L1,C",
    "2020-01-01 00:00:01,2,L1,C",
    "2020-01-01 00:00:01,4,L1,C",
    "2020-01-01 00:00:02,NA,L1,C",
    "2020-01-01 00:01:00,100,L1,C",
    "2020-01-01 00:01:01,3,L1,C"
  ), path)
  on.exit(unlink(path))

  result <- read_timeseries_file(path)
  expect_equal(nrow(result$data), 6L)
  expect_equal(result$diagnostics$duplicate_timestamp_rows, 1L)
  expect_equal(result$diagnostics$gap_count, 1L)
  expect_equal(sum(is.na(result$data$value)), 1L)
  expect_equal(result$data$.gap_id[5L], result$data$.gap_id[4L] + 1L)
  expect_true(all(diff(as.numeric(result$data$datetime)) >= 0))
})

test_that("CSV and RDS loaders identify supported schemas", {
  csv <- tempfile(fileext = ".csv")
  rds <- tempfile(fileext = ".rds")
  on.exit(unlink(c(csv, rds)))
  table <- data.frame(
    timestamp = c("2022-03-01 00:00:00", "2022-03-01 00:00:01"),
    value = c(1.5, 2.5),
    variable = c("temperature", "temperature"),
    unit = c("C", "C")
  )
  utils::write.csv(table, csv, row.names = FALSE)
  saveRDS(table, rds)

  expect_equal(read_timeseries_file(csv)$data$value, c(1.5, 2.5))
  expect_equal(read_timeseries_file(rds)$data$variable,
               c("temperature", "temperature"))
})

test_that("ambiguous schemas require an explicit mapping", {
  table <- data.frame(
    timestamp = c("2020-01-01", "2020-01-02"),
    date = c("2020-01-01", "2020-01-02"),
    temperature = c(1, 2),
    pressure = c(10, 11)
  )
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  utils::write.csv(table, path, row.names = FALSE)

  expect_error(read_timeseries_file(path), "Ambiguous datetime columns")
  mapped <- read_timeseries_file(path, list(
    datetime = "timestamp",
    value = "temperature"
  ))
  expect_equal(mapped$data$value, c(1, 2))
})

test_that("timestamps respect the selected timezone", {
  parsed <- parse_timeseries_datetime(
    c("2020-01-01 12:00:00", "2020-07-01 12:00:00"),
    "Europe/Vienna",
    "timezone.csv"
  )
  expect_equal(format(parsed, "%H:%M:%S", tz = "Europe/Vienna"),
               c("12:00:00", "12:00:00"))
  expect_equal(
    as.numeric(diff(parsed), units = "secs"),
    182L * 86400 - 3600
  )
})

test_that("aligned summary levels preserve sufficient statistics and extrema", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  times <- as.POSIXct("2020-01-01", tz = "UTC") + 0:599
  values <- sin(seq_along(times) / 10)
  values[121L] <- 50
  normalized <- normalize_timeseries(
    data.frame(
      datetime = format(times, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      value = values,
      logger = "L1"
    ),
    path,
    list(
      datetime = "datetime",
      value = "value",
      variable = NULL,
      unit = NULL,
      group_cols = "logger",
      timezone = "UTC"
    )
  )
  summaries <- build_resolution_summaries(normalized$data, max_levels = 4L)
  level_zero <- summaries$levels[["0"]]
  level_one <- summaries$levels[["1"]]
  expect_gt(nrow(level_zero), nrow(level_one))

  parent <- level_one[1L, ]
  children <- level_zero[
    level_zero$.series_id == parent$.series_id &
      level_zero$.gap_id == parent$.gap_id &
      level_zero$.bin_start >= parent$.bin_start &
      level_zero$.bin_end <= parent$.bin_end,
  ]
  expect_equal(parent$.count, sum(children$.count))
  expect_equal(parent$.sum, sum(children$.sum), tolerance = 1e-10)
  expect_equal(parent$.min, min(children$.min))
  expect_equal(parent$.max, max(children$.max))
  expect_true(parent$.min_time %in% children$.min_time)
  expect_true(parent$.max_time %in% children$.max_time)
  expect_equal(max(level_zero$.max), 50)
})

test_that("viewport selection uses raw data when feasible and coarser levels otherwise", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  times <- as.POSIXct("2020-01-01", tz = "UTC") + 0:999
  normalized <- normalize_timeseries(
    data.frame(datetime = times, value = seq_along(times)),
    path,
    list(
      datetime = "datetime",
      value = "value",
      variable = NULL,
      unit = NULL,
      group_cols = character(),
      timezone = "UTC"
    )
  )
  summaries <- build_resolution_summaries(normalized$data, max_levels = 8L)
  selection_raw <- select_plot_resolution(
    summaries, times[1L], times[5L], plot_width_px = 100
  )
  expect_identical(selection_raw$level, "raw")

  selection_summary <- select_plot_resolution(
    summaries, times[1L], times[length(times)], plot_width_px = 10
  )
  expect_true(is.numeric(selection_summary$level))
  expect_gt(selection_summary$level, 0)
  expect_error(select_plot_resolution(
    summaries, times[2L], times[1L], plot_width_px = 10
  ), "finite time range")
  expect_error(select_plot_resolution(
    summaries, times[1L], times[2L], plot_width_px = 0
  ), "positive number of pixels")
})

test_that("viewport query returns exact raw edge observations", {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path))
  times <- as.POSIXct("2020-01-01", tz = "UTC") + 0:999
  normalized <- normalize_timeseries(
    data.frame(datetime = times, value = seq_along(times)),
    path,
    list(
      datetime = "datetime",
      value = "value",
      variable = NULL,
      unit = NULL,
      group_cols = character(),
      timezone = "UTC"
    )
  )
  summaries <- build_resolution_summaries(normalized$data, max_levels = 8L)
  prepared <- list(
    path = path,
    mapping = list(),
    cache_path = tempfile(fileext = ".rds"),
    summary = summaries,
    raw = normalized$data,
    diagnostics = normalized$diagnostics
  )
  query <- query_timeseries_viewport(
    prepared,
    times[11L] + 0.25,
    times[990L] + 0.25,
    plot_width_px = 10
  )
  expect_gt(query$resolution$level, 0)
  expect_true(all(query$raw$datetime >= times[11L] + 0.25))
  expect_true(all(query$raw$datetime <= times[990L] + 0.25))
  expect_gt(nrow(query$summary), 0L)
})

test_that("persistent summaries are invalidated when source files change", {
  path <- tempfile(fileext = ".csv")
  cache_dir <- tempfile("summary-cache-")
  dir.create(cache_dir)
  on.exit(unlink(c(path, cache_dir), recursive = TRUE))
  writeLines(c(
    "datetime,value",
    "2020-01-01 00:00:00,1",
    "2020-01-01 00:00:01,2"
  ), path)
  first <- load_preprocessed_timeseries(path, cache_dir = cache_dir)
  first_cache <- first$cache_path
  Sys.sleep(1.1)
  writeLines(c(
    "datetime,value",
    "2020-01-01 00:00:00,1",
    "2020-01-01 00:00:01,2",
    "2020-01-01 00:00:02,3"
  ), path)
  second <- load_preprocessed_timeseries(path, cache_dir = cache_dir)
  expect_false(identical(first_cache, second$cache_path))
  expect_equal(sum(second$summary$levels[["0"]]$.count), 3)
})
