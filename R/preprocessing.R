PREPROCESSING_VERSION <- "1.0.0"

PREPROCESSING_ID_COLUMNS <- c(
  "logger", "position", "source", "trm", "tree", "core_id", "probe_set",
  "depth_bottom", "depth_top", "station", "sensor", "channel"
)

PREPROCESSING_MEMORY_CACHE <- new.env(parent = emptyenv())
PREPROCESSING_MEMORY_CACHE_ORDER <- character()

preprocessing_abort <- function(message) {
  stop(structure(
    list(message = message, call = NULL),
    class = c("preprocessing_error", "error", "condition")
  ))
}

discover_data_files <- function(directories = c("data_input", "data")) {
  existing <- directories[dir.exists(directories)]
  if (!length(existing)) {
    preprocessing_abort(paste0(
      "No input directories were found. Create one of: ",
      paste(directories, collapse = ", ")
    ))
  }

  paths <- unlist(lapply(existing, function(directory) {
    list.files(
      directory,
      pattern = "\\.(csv|rds)$",
      full.names = TRUE,
      ignore.case = TRUE,
      recursive = FALSE
    )
  }), use.names = FALSE)
  paths <- paths[!startsWith(basename(paths), ".")]
  if (!length(paths)) {
    return(data.frame(
      path = character(),
      dataset = character(),
      format = character(),
      size_bytes = numeric(),
      modified = as.POSIXct(character(), tz = "UTC"),
      stringsAsFactors = FALSE
    ))
  }

  info <- file.info(paths)
  keep <- !is.na(info$isdir) & !info$isdir & !is.na(info$size)
  paths <- paths[keep]
  info <- info[keep, , drop = FALSE]
  if (!length(paths)) {
    return(data.frame(
      path = character(),
      dataset = character(),
      format = character(),
      size_bytes = numeric(),
      modified = as.POSIXct(character(), tz = "UTC"),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    dataset = tools::file_path_sans_ext(basename(paths)),
    format = tolower(tools::file_ext(paths)),
    size_bytes = as.numeric(info$size),
    modified = as.POSIXct(info$mtime, tz = "UTC"),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

read_timeseries_table <- function(path, header_only = FALSE) {
  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "csv")) {
    if (header_only) {
      names <- names(data.table::fread(
        path,
        nrows = 0L,
        check.names = FALSE,
        data.table = FALSE,
        showProgress = FALSE,
        nThread = 1L
      ))
      return(names)
    }
    return(data.table::fread(
      path,
      check.names = FALSE,
      data.table = FALSE,
      showProgress = FALSE,
      nThread = 1L,
      na.strings = c("", "NA", "NaN")
    ))
  }
  if (!identical(extension, "rds")) {
    preprocessing_abort(paste0(
      "Unsupported file type: ", extension, ". Choose a CSV or RDS file."
    ))
  }

  object <- tryCatch(
    readRDS(path),
    error = function(error) {
      preprocessing_abort(paste0(
        "Could not read ", basename(path), ": ", conditionMessage(error)
      ))
    }
  )
  if (!is.data.frame(object)) {
    preprocessing_abort(paste0(
      basename(path), " contains a ", paste(class(object), collapse = "/"),
      " object; expected a data frame or tibble."
    ))
  }
  if (header_only) names(object) else as.data.frame(object, stringsAsFactors = FALSE)
}

infer_timeseries_schema <- function(data, path, mapping = list()) {
  columns <- names(data)
  if (is.null(columns) || !length(columns)) {
    preprocessing_abort(paste0(basename(path), " has no columns."))
  }
  normalized_names <- tolower(gsub("[^a-z0-9]+", "_", columns))

  pick <- function(field, candidates, preferred = character()) {
    explicit <- mapping[[field]]
    if (!is.null(explicit) && length(explicit) && nzchar(explicit)) {
      if (!explicit %in% columns) {
        preprocessing_abort(paste0(
          "Mapped ", field, " column '", explicit,
          "' does not exist in ", basename(path), "."
        ))
      }
      return(explicit)
    }

    matched <- columns[normalized_names %in% candidates]
    if (length(preferred)) {
      preferred_matches <- intersect(preferred, matched)
      if (length(preferred_matches) == 1L) return(preferred_matches)
    }
    if (length(matched) == 1L) return(matched)
    if (length(matched) > 1L) {
      preprocessing_abort(paste0(
        "Ambiguous ", field, " columns in ", basename(path), ": ",
        paste(matched, collapse = ", "),
        ". Specify the column explicitly in the schema mapping."
      ))
    }
    NULL
  }

  datetime <- pick(
    "datetime",
    c("datetime", "timestamp", "time", "date_time", "date", "sample_time")
  )
  if (is.null(datetime)) {
    preprocessing_abort(paste0(
      "No timestamp column was identified in ", basename(path),
      ". Map a column named datetime, timestamp, time, or date."
    ))
  }

  value <- pick(
    "value",
    c("value", "measurement", "reading", "observation", "result"),
    preferred = "value"
  )
  if (is.null(value)) {
    if (ncol(data) > 0L && !is.null(dim(data))) {
      numeric_columns <- columns[vapply(data, function(column) {
        is.numeric(column) || is.integer(column)
      }, logical(1))]
      numeric_columns <- setdiff(numeric_columns, datetime)
      if (length(numeric_columns) == 1L) value <- numeric_columns
      if (length(numeric_columns) > 1L) {
        preprocessing_abort(paste0(
          "Ambiguous measurement columns in ", basename(path), ": ",
          paste(numeric_columns, collapse = ", "),
          ". Specify the value column in the schema mapping."
        ))
      }
    }
  }
  if (is.null(value)) {
    preprocessing_abort(paste0(
      "No numeric measurement column was identified in ", basename(path),
      ". Map the value column explicitly."
    ))
  }

  variable <- pick(
    "variable",
    c("variable", "parameter", "measurement_type", "measurand"),
    preferred = "variable"
  )
  unit <- pick("unit", c("unit", "units"), preferred = "unit")
  group_columns <- mapping$group_cols
  if (is.null(group_columns)) {
    group_columns <- intersect(PREPROCESSING_ID_COLUMNS, columns)
  }
  unknown_groups <- setdiff(group_columns, columns)
  if (length(unknown_groups)) {
    preprocessing_abort(paste0(
      "Unknown grouping columns: ", paste(unknown_groups, collapse = ", ")
    ))
  }
  group_columns <- setdiff(group_columns, c(datetime, value, variable, unit))

  list(
    datetime = datetime,
    value = value,
    variable = variable,
    unit = unit,
    group_cols = group_columns,
    timezone = mapping$timezone %||% "UTC"
  )
}

`%||%` <- function(value, default) {
  if (is.null(value)) default else value
}

parse_timeseries_datetime <- function(value, timezone, path) {
  if (length(timezone) != 1L || is.na(timezone) ||
      !timezone %in% c("UTC", "GMT", OlsonNames())) {
    preprocessing_abort(paste0("Invalid timezone: ", timezone))
  }
  if (inherits(value, "POSIXct")) {
    parsed <- as.POSIXct(value, tz = timezone)
  } else if (inherits(value, "POSIXlt")) {
    parsed <- as.POSIXct(value, tz = timezone)
  } else if (inherits(value, "Date")) {
    parsed <- as.POSIXct(value, tz = timezone)
  } else if (is.numeric(value)) {
    parsed <- as.POSIXct(value, origin = "1970-01-01", tz = timezone)
  } else {
    text <- trimws(as.character(value))
    text[!nzchar(text)] <- NA_character_
    formats <- c(
      "%Y-%m-%dT%H:%M:%OS%z",
      "%Y-%m-%dT%H:%M:%OS",
      "%Y-%m-%d %H:%M:%OS%z",
      "%Y-%m-%d %H:%M:%OS",
      "%Y-%m-%d",
      "%d.%m.%Y %H:%M:%OS",
      "%d.%m.%Y",
      "%m/%d/%Y %H:%M:%OS",
      "%m/%d/%Y"
    )
    parsed <- as.POSIXct(rep(NA_real_, length(text)),
                         origin = "1970-01-01", tz = timezone)
    for (format in formats) {
      missing <- is.na(parsed) & !is.na(text)
      if (!any(missing)) break
      candidate <- as.POSIXct(
        text[missing],
        format = format,
        tz = timezone
      )
      parsed[missing] <- candidate
    }
  }
  if (!any(!is.na(parsed))) {
    preprocessing_abort(paste0(
      "No timestamps in ", basename(path),
      " could be parsed in timezone ", timezone, "."
    ))
  }
  parsed
}

make_timeseries_series_id <- function(data, columns) {
  if (!length(columns)) {
    return(rep("all", nrow(data)))
  }
  encoded <- lapply(columns, function(column) {
    value <- as.character(data[[column]])
    missing <- is.na(value)
    encoded_value <- paste0(nchar(value, type = "bytes"), ":", value)
    encoded_value[missing] <- "-1:"
    encoded_value
  })
  do.call(paste, c(encoded, sep = "|"))
}

normalize_timeseries <- function(data, path, schema, source_id = path,
                                 gap_factor = 5) {
  required <- c(schema$datetime, schema$value, schema$group_cols)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    preprocessing_abort(paste0(
      "Required columns are missing from ", basename(path), ": ",
      paste(missing, collapse = ", ")
    ))
  }
  timezone <- schema$timezone %||% "UTC"
  datetime <- parse_timeseries_datetime(data[[schema$datetime]], timezone, path)
  value_text <- data[[schema$value]]
  value <- suppressWarnings(as.numeric(as.character(value_text)))
  invalid_value <- !is.na(value_text) & is.na(value)
  infinite_value <- !is.na(value) & !is.finite(value)
  value[infinite_value] <- NA_real_
  value[invalid_value] <- NA_real_

  dataset <- tools::file_path_sans_ext(basename(path))
  variable <- if (!is.null(schema$variable) && nzchar(schema$variable)) {
    as.character(data[[schema$variable]])
  } else {
    rep(dataset, nrow(data))
  }
  unit <- if (!is.null(schema$unit) && nzchar(schema$unit)) {
    as.character(data[[schema$unit]])
  } else {
    rep(NA_character_, nrow(data))
  }
  out <- data.frame(
    dataset = rep(dataset, nrow(data)),
    source_file = rep(normalizePath(path, winslash = "/", mustWork = FALSE), nrow(data)),
    datetime = datetime,
    date = as.Date(datetime, tz = timezone),
    time_zone = rep(timezone, nrow(data)),
    variable = variable,
    value = value,
    unit = unit,
    stringsAsFactors = FALSE
  )
  for (column in schema$group_cols) {
    if (!column %in% names(out)) out[[column]] <- data[[column]]
  }

  series_columns <- c("source_file", "dataset", "variable", "unit", schema$group_cols)
  out$.series_id <- make_timeseries_series_id(out, series_columns)
  valid_time <- !is.na(out$datetime)
  dropped_time <- sum(!valid_time)
  out <- out[valid_time, , drop = FALSE]
  if (!nrow(out)) {
    preprocessing_abort(paste0(
      "No rows with valid timestamps remain in ", basename(path), "."
    ))
  }
  data.table::setorderv(out, c(".series_id", "datetime"), na.last = TRUE)
  row.names(out) <- NULL

  duplicate_key <- paste(out$.series_id, as.numeric(out$datetime), sep = "\r")
  duplicate_count <- sum(duplicated(duplicate_key))
  gaps <- detect_timeseries_gaps(out, gap_factor = gap_factor)
  out <- gaps$data
  valid_value_count <- sum(is.finite(out$value), na.rm = TRUE)
  if (!valid_value_count) {
    preprocessing_abort(paste0(
      "No finite numeric measurements remain in ", basename(path),
      ". Check the mapped value column."
    ))
  }

  diagnostics <- gaps$diagnostics
  diagnostics$input_rows <- nrow(data)
  diagnostics$valid_timestamp_rows <- nrow(out)
  diagnostics$dropped_timestamp_rows <- dropped_time
  diagnostics$valid_value_rows <- valid_value_count
  diagnostics$invalid_value_rows <- sum(invalid_value)
  diagnostics$non_finite_value_rows <- sum(infinite_value)
  diagnostics$duplicate_timestamp_rows <- duplicate_count
  attr(out, "diagnostics") <- diagnostics
  list(
    data = out,
    schema = schema,
    diagnostics = diagnostics,
    source_id = source_id
  )
}

detect_timeseries_gaps <- function(data, gap_factor = 5) {
  if (!is.finite(gap_factor) || gap_factor <= 1) {
    preprocessing_abort("The gap factor must be a finite number greater than 1.")
  }
  if (!all(c(".series_id", "datetime", "value") %in% names(data))) {
    preprocessing_abort(
      "Gap detection requires .series_id, datetime, and value columns."
    )
  }
  out <- data.table::as.data.table(data)
  data.table::setorderv(out, c(".series_id", "datetime"))
  out[, .time_seconds := as.numeric(datetime)]
  out[, .delta_seconds := .time_seconds - data.table::shift(.time_seconds),
      by = .series_id]
  cadence <- out[.delta_seconds > 0,
                 .(cadence_seconds = stats::median(.delta_seconds)),
                 by = .series_id]
  cadence[, .gap_threshold_seconds := cadence_seconds * gap_factor]
  out <- merge(out, cadence, by = ".series_id", all.x = TRUE, sort = FALSE)
  data.table::setorderv(out, c(".series_id", "datetime"))
  out[, .gap_before := !is.na(.delta_seconds) &
        !is.na(.gap_threshold_seconds) &
        .delta_seconds > .gap_threshold_seconds]
  out[, .segment_break := .gap_before |
        !is.finite(value) |
        data.table::shift(!is.finite(value), fill = FALSE),
      by = .series_id]
  out[, .gap_id := cumsum(.segment_break), by = .series_id]
  out[, .segment_break := NULL]
  out[, .time_seconds := NULL]
  diagnostics <- out[, .(
    row_count = .N,
    distinct_timestamp_count = data.table::uniqueN(datetime),
    first_timestamp = min(datetime),
    last_timestamp = max(datetime),
    cadence_seconds = if (all(is.na(cadence_seconds))) NA_real_
    else stats::median(cadence_seconds, na.rm = TRUE),
    gap_threshold_seconds = if (all(is.na(.gap_threshold_seconds))) NA_real_
    else stats::median(.gap_threshold_seconds, na.rm = TRUE),
    gap_count = sum(.gap_before, na.rm = TRUE)
  ), by = .series_id]
  out[, .delta_seconds := NULL]
  out[, .gap_threshold_seconds := NULL]
  list(data = as.data.frame(out), diagnostics = as.data.frame(diagnostics))
}

read_timeseries_file <- function(path, mapping = list(), gap_factor = 5) {
  if (!file.exists(path)) {
    preprocessing_abort(paste0("Input file does not exist: ", path))
  }
  table <- read_timeseries_table(path)
  if (!nrow(table)) {
    preprocessing_abort(paste0("Input file is empty: ", basename(path)))
  }
  schema <- infer_timeseries_schema(table, path, mapping)
  normalize_timeseries(
    table,
    path,
    schema,
    source_id = path,
    gap_factor = gap_factor
  )
}

source_cache_identity <- function(path, mapping, configuration = list()) {
  info <- file.info(path)
  if (is.na(info$size) || is.na(info$mtime)) {
    preprocessing_abort(paste0("Could not inspect source file: ", path))
  }
  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    size = as.numeric(info$size),
    modified = as.numeric(info$mtime),
    version = PREPROCESSING_VERSION,
    mapping = mapping,
    configuration = configuration
  )
}

preprocessing_cache_path <- function(path, mapping, cache_dir,
                                     configuration = list()) {
  identity <- source_cache_identity(path, mapping, configuration)
  payload <- serialize(identity, NULL, version = 2)
  temporary <- tempfile("preprocess-id-")
  on.exit(unlink(temporary), add = TRUE)
  writeBin(payload, temporary)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(cache_dir)) {
    preprocessing_abort(paste0("Cannot create cache directory: ", cache_dir))
  }
  key <- unname(tools::md5sum(temporary))
  list(path = file.path(cache_dir, paste0(key, ".rds")), identity = identity)
}

aggregate_summary_level <- function(data, width_by_series, level) {
  dt <- data.table::copy(data.table::as.data.table(data))
  dt[, .base_seconds := unname(width_by_series[.series_id])]
  dt[, .width_seconds := .base_seconds * (2 ^ level)]
  dt[, .bin_start_seconds := floor(as.numeric(datetime) / .width_seconds) *
        .width_seconds]
  dt[, .bin_start := as.POSIXct(
    .bin_start_seconds,
    origin = "1970-01-01",
    tz = "UTC"
  )]
  dt <- dt[is.finite(value) & !is.na(datetime)]
  if (!nrow(dt)) return(data.frame())

  grouped <- dt[, .(
    .count = .N,
    .sum = sum(value),
    .mean = mean(value),
    .min = min(value),
    .min_time = datetime[which.min(value)],
    .max = max(value),
    .max_time = datetime[which.max(value)]
  ), by = c(".series_id", ".gap_id", ".base_seconds", ".width_seconds",
            ".bin_start")]
  grouped[, .bin_end := .bin_start + .width_seconds]
  grouped[, .level := as.integer(level)]
  as.data.frame(grouped)
}

merge_summary_level <- function(previous, level) {
  if (!nrow(previous)) return(data.frame())
  dt <- data.table::copy(data.table::as.data.table(previous))
  dt[, .width_seconds := .base_seconds * (2 ^ level)]
  dt[, .bin_start := as.POSIXct(
    floor(as.numeric(.bin_start) / .width_seconds) * .width_seconds,
    origin = "1970-01-01",
    tz = "UTC"
  )]
  grouped <- dt[, .(
    .count = sum(.count),
    .sum = sum(.sum),
    .mean = sum(.sum) / sum(.count),
    .min = min(.min),
    .min_time = .min_time[which.min(.min)],
    .max = max(.max),
    .max_time = .max_time[which.max(.max)]
  ), by = c(".series_id", ".gap_id", ".base_seconds", ".width_seconds",
            ".bin_start")]
  grouped[, .bin_end := .bin_start + .width_seconds]
  grouped[, .level := as.integer(level)]
  as.data.frame(grouped)
}

build_resolution_summaries <- function(normalized, max_levels = 12L,
                                       observations_per_base_bin = 4,
                                       min_bin_seconds = 1,
                                       gap_factor = 5) {
  if (!all(c(".series_id", ".gap_id", "datetime", "value") %in%
           names(normalized))) {
    preprocessing_abort(
      "Summary generation requires .series_id, .gap_id, datetime, and value."
    )
  }
  if (!is.finite(max_levels) || max_levels < 1 ||
      !is.finite(observations_per_base_bin) || observations_per_base_bin < 1 ||
      !is.finite(min_bin_seconds) || min_bin_seconds <= 0) {
    preprocessing_abort("Invalid resolution-generation configuration.")
  }
  dt <- data.table::as.data.table(normalized)
  cadence <- if ("cadence_seconds" %in% names(dt)) {
    dt[, .(
      .cadence_seconds = if (any(is.finite(cadence_seconds))) {
        stats::median(cadence_seconds[is.finite(cadence_seconds)])
      } else {
        NA_real_
      }
    ), by = .series_id]
  } else {
    dt[, .(.cadence_seconds = NA_real_), by = .series_id]
  }
  cadence[, .base_seconds := pmax(
    min_bin_seconds,
    ifelse(is.finite(.cadence_seconds),
           .cadence_seconds * observations_per_base_bin,
           min_bin_seconds)
  )]
  width_by_series <- stats::setNames(
    cadence$.base_seconds,
    cadence$.series_id
  )
  level_zero <- aggregate_summary_level(dt, width_by_series, 0L)
  if (!nrow(level_zero)) {
    preprocessing_abort("No finite observations are available to summarize.")
  }

  levels <- vector("list", max_levels + 1L)
  levels[[1L]] <- level_zero
  for (level in seq_len(max_levels)) {
    next_level <- merge_summary_level(levels[[level]], level)
    if (!nrow(next_level)) break
    levels[[level + 1L]] <- next_level
    if (nrow(next_level) <= 1L) {
      levels <- levels[seq_len(level + 1L)]
      break
    }
  }
  levels <- Filter(Negate(is.null), levels)
  names(levels) <- as.character(seq_along(levels) - 1L)

  non_series_columns <- c(
    "datetime", "date", "value", "cadence_seconds", ".gap_before",
    ".gap_id", ".delta_seconds", ".time_seconds", ".gap_threshold_seconds"
  )
  series_columns <- setdiff(names(normalized), non_series_columns)
  series_columns <- c(".series_id", setdiff(series_columns, ".series_id"))
  series_map <- unique(normalized[series_columns])
  row.names(series_map) <- NULL
  time_ranges <- data.table::as.data.table(normalized)[, .(
    .start = min(datetime),
    .end = max(datetime)
  ), by = .series_id]

  list(
    version = PREPROCESSING_VERSION,
    generated_at = as.POSIXct(Sys.time(), tz = "UTC"),
    source_identity = NULL,
    configuration = list(
      max_levels = as.integer(max_levels),
      observations_per_base_bin = observations_per_base_bin,
      min_bin_seconds = min_bin_seconds,
      alignment = "UTC Unix epoch",
      gap_factor = gap_factor
    ),
    raw_object_bytes = as.numeric(object.size(normalized)),
    series_map = series_map,
    time_ranges = as.data.frame(time_ranges),
    diagnostics = attr(normalized, "diagnostics"),
    levels = levels
  )
}

save_rds_atomic <- function(object, path) {
  directory <- dirname(path)
  if (!dir.exists(directory) &&
      !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
    preprocessing_abort(paste0("Cannot create output directory: ", directory))
  }
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "-"),
    tmpdir = directory,
    fileext = ".tmp"
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = FALSE)
  if (file.exists(path)) {
    existing <- tryCatch(readRDS(path), error = function(error) NULL)
    if (!is.null(existing)) {
      unlink(temporary)
      return(invisible(path))
    }
    if (!file.remove(path)) {
      preprocessing_abort(paste0(
        "Cannot replace a corrupt cache file: ", path
      ))
    }
  }
  if (!file.rename(temporary, path)) {
    if (file.exists(path)) {
      winner <- tryCatch(readRDS(path), error = function(error) NULL)
      if (!is.null(winner)) return(invisible(path))
    }
    preprocessing_abort(paste0(
      "Could not atomically write preprocessing cache: ", path
    ))
  }
  invisible(path)
}

preprocessing_memory_limit <- function() {
  value <- getOption("dataExplorer.preprocessing_cache_bytes", 512 * 1024^2)
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value < 0) {
    preprocessing_abort(
      "Option dataExplorer.preprocessing_cache_bytes must be non-negative."
    )
  }
  as.numeric(value)
}

prune_preprocessing_cache <- function(cache_dir, current_path = NULL) {
  limit <- getOption(
    "dataExplorer.preprocessing_disk_cache_bytes",
    5 * 1024^3
  )
  if (!is.numeric(limit) || length(limit) != 1L ||
      !is.finite(limit) || limit < 0) {
    preprocessing_abort(
      "Option dataExplorer.preprocessing_disk_cache_bytes must be non-negative."
    )
  }
  files <- list.files(
    cache_dir,
    pattern = "\\.rds$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(files)) return(invisible(0L))
  info <- file.info(files)
  total <- sum(info$size, na.rm = TRUE)
  candidates <- order(info$mtime)
  removed <- 0L
  for (index in candidates) {
    if (total <= limit) break
    path <- files[index]
    if (!is.null(current_path) &&
        identical(normalizePath(path, winslash = "/", mustWork = FALSE),
                  normalizePath(current_path, winslash = "/", mustWork = FALSE))) {
      next
    }
    if (file.remove(path)) {
      total <- total - info$size[index]
      removed <- removed + 1L
    } else {
      warning(
        paste0("Could not remove expired summary cache file: ", path),
        call. = FALSE
      )
    }
  }
  invisible(removed)
}

preprocessing_memory_get <- function(key) {
  if (!exists(key, envir = PREPROCESSING_MEMORY_CACHE, inherits = FALSE)) {
    return(NULL)
  }
  PREPROCESSING_MEMORY_CACHE_ORDER <<- c(
    key,
    setdiff(PREPROCESSING_MEMORY_CACHE_ORDER, key)
  )
  get(key, envir = PREPROCESSING_MEMORY_CACHE, inherits = FALSE)
}

preprocessing_memory_set <- function(key, value) {
  limit <- preprocessing_memory_limit()
  value_size <- as.numeric(object.size(value))
  if (value_size > limit) {
    value$raw <- NULL
    value_size <- as.numeric(object.size(value))
    if (value_size > limit) {
      if (exists(key, envir = PREPROCESSING_MEMORY_CACHE, inherits = FALSE)) {
        rm(list = key, envir = PREPROCESSING_MEMORY_CACHE)
        PREPROCESSING_MEMORY_CACHE_ORDER <<- setdiff(
          PREPROCESSING_MEMORY_CACHE_ORDER,
          key
        )
      }
      return(value)
    }
  }
  assign(key, value, envir = PREPROCESSING_MEMORY_CACHE)
  PREPROCESSING_MEMORY_CACHE_ORDER <<- c(
    key,
    setdiff(PREPROCESSING_MEMORY_CACHE_ORDER, key)
  )
  sizes <- vapply(
    as.list(PREPROCESSING_MEMORY_CACHE),
    function(item) as.numeric(object.size(item)),
    numeric(1)
  )
  total <- sum(sizes)
  while (total > preprocessing_memory_limit() &&
         length(PREPROCESSING_MEMORY_CACHE_ORDER) > 1L) {
    oldest <- tail(PREPROCESSING_MEMORY_CACHE_ORDER, 1L)
    if (identical(oldest, key)) break
    rm(list = oldest, envir = PREPROCESSING_MEMORY_CACHE)
    PREPROCESSING_MEMORY_CACHE_ORDER <<- head(
      PREPROCESSING_MEMORY_CACHE_ORDER,
      -1L
    )
    total <- sum(vapply(
      as.list(PREPROCESSING_MEMORY_CACHE),
      function(item) as.numeric(object.size(item)),
      numeric(1)
    ))
  }
  if (total > limit) {
    value$raw <- NULL
    assign(key, value, envir = PREPROCESSING_MEMORY_CACHE)
    return(value)
  }
  value
}

load_preprocessed_timeseries <- function(path, mapping = list(),
                                         cache_dir = file.path(
                                           "data_processed", ".cache"
                                         ),
                                         max_levels = 12L,
                                         gap_factor = 5) {
  cache <- preprocessing_cache_path(
    path,
    mapping,
    cache_dir,
    configuration = list(max_levels = max_levels, gap_factor = gap_factor)
  )
  key <- cache$path
  cached <- preprocessing_memory_get(key)
  if (!is.null(cached)) return(cached)

  summary <- NULL
  if (file.exists(cache$path)) {
    summary <- tryCatch(readRDS(cache$path), error = function(error) {
      warning(
        paste0("Ignoring corrupt preprocessing cache: ", conditionMessage(error)),
        call. = FALSE
      )
      NULL
    })
    if (!is.null(summary) &&
        !identical(summary$source_identity, cache$identity)) {
      summary <- NULL
    }
  }

  normalized <- NULL
  if (is.null(summary)) {
    normalized <- read_timeseries_file(path, mapping, gap_factor = gap_factor)
    source_after_read <- source_cache_identity(
      path,
      mapping,
      configuration = list(max_levels = max_levels, gap_factor = gap_factor)
    )
    if (!identical(cache$identity, source_after_read)) {
      preprocessing_abort(paste0(
        basename(path), " changed while it was being read. Retry after the ",
        "file has finished being written."
      ))
    }
    summary <- build_resolution_summaries(
      normalized$data,
      max_levels = max_levels,
      gap_factor = gap_factor
    )
    summary$source_identity <- cache$identity
    save_rds_atomic(summary, cache$path)
    prune_preprocessing_cache(cache_dir, cache$path)
  }
  value <- list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    mapping = mapping,
    gap_factor = gap_factor,
    cache_path = cache$path,
    summary = summary,
    raw = if (is.null(normalized)) NULL else normalized$data,
    diagnostics = if (is.null(normalized)) summary$diagnostics
    else normalized$diagnostics
  )
  value <- preprocessing_memory_set(key, value)
  preprocessing_memory_get(key) %||% value
}

load_cached_raw_timeseries <- function(preprocessed) {
  if (!is.null(preprocessed$raw)) return(preprocessed$raw)
  normalized <- read_timeseries_file(
    preprocessed$path,
    preprocessed$mapping,
    gap_factor = preprocessed$gap_factor %||% 5
  )
  key <- preprocessed$cache_path
  current <- preprocessing_memory_get(key)
  current$raw <- normalized$data
  preprocessing_memory_set(key, current)
  normalized$data
}

select_plot_resolution <- function(summary, start, end, plot_width_px,
                                   target_per_pixel = 2,
                                   selected_series = NULL,
                                   raw_row_count = NULL) {
  start_seconds <- as.numeric(start)
  end_seconds <- as.numeric(end)
  if (length(start_seconds) != 1L || length(end_seconds) != 1L ||
      !is.finite(start_seconds) || !is.finite(end_seconds) ||
      start_seconds >= end_seconds) {
    preprocessing_abort("Viewport start and end must define a finite time range.")
  }
  if (length(plot_width_px) != 1L || !is.finite(plot_width_px) ||
      plot_width_px <= 0) {
    preprocessing_abort("Plot width must be a positive number of pixels.")
  }
  if (length(target_per_pixel) != 1L || !is.finite(target_per_pixel) ||
      target_per_pixel <= 0) {
    preprocessing_abort("Target representations per pixel must be positive.")
  }
  budget <- max(1L, floor(plot_width_px * target_per_pixel))
  levels <- summary$levels
  if (!length(levels)) {
    return(list(level = "empty", budget = budget, reason = "no summary levels"))
  }
  selected <- selected_series
  candidate_rows <- lapply(levels, function(level_data) {
    keep <- level_data$.bin_end > start_seconds &
      as.numeric(level_data$.bin_start) < end_seconds
    if (!is.null(selected)) keep <- keep & level_data$.series_id %in% selected
    level_data[keep, , drop = FALSE]
  })
  level_zero_count <- if (length(candidate_rows)) {
    sum(candidate_rows[[1L]]$.count)
  } else {
    0
  }
  if (!is.null(raw_row_count)) level_zero_count <- raw_row_count
  raw_available <- summary$raw_object_bytes <= preprocessing_memory_limit()
  if (level_zero_count <= budget && raw_available) {
    return(list(
      level = "raw",
      budget = budget,
      estimated_rows = level_zero_count,
      reason = "raw observations fit the viewport point and memory budgets"
    ))
  }

  representations <- vapply(candidate_rows, nrow, integer(1)) * 3L
  selected_level <- which(representations <= budget)[1L]
  if (is.na(selected_level)) selected_level <- length(levels)
  list(
    level = as.integer(names(levels)[selected_level]),
    budget = budget,
    estimated_rows = representations[selected_level],
    reason = if (representations[selected_level] <= budget) {
      "finest summary level within the point budget"
    } else if (level_zero_count <= budget && !raw_available) {
      "raw subset exceeds the configured process memory budget"
    } else {
      "coarsest available level; series count exceeds the point budget"
    }
  )
}

query_timeseries_viewport <- function(preprocessed, start, end,
                                      plot_width_px, target_per_pixel = 2,
                                      selected_series = NULL) {
  choice <- select_plot_resolution(
    preprocessed$summary,
    start,
    end,
    plot_width_px,
    target_per_pixel,
    selected_series
  )
  if (identical(choice$level, "empty")) {
    return(list(
      raw = data.frame(),
      summary = data.frame(),
      resolution = choice,
      diagnostics = list(returned_rows = 0L, partial_edge_bins = 0L)
    ))
  }
  start_seconds <- as.numeric(start)
  end_seconds <- as.numeric(end)
  if (identical(choice$level, "raw")) {
    raw <- load_cached_raw_timeseries(preprocessed)
    dt <- data.table::as.data.table(raw)
    data.table::setkeyv(dt, c(".series_id", "datetime"))
    ids <- selected_series %||% unique(dt$.series_id)
    ranges <- data.table::data.table(
      .series_id = ids,
      start = as.POSIXct(start_seconds, origin = "1970-01-01", tz = "UTC"),
      end = as.POSIXct(end_seconds, origin = "1970-01-01", tz = "UTC")
    )
    rows <- dt[ranges, on = .(.series_id, datetime >= start, datetime < end),
               nomatch = 0L, allow.cartesian = TRUE]
    return(list(
      raw = as.data.frame(rows),
      summary = data.frame(),
      resolution = choice,
      diagnostics = list(returned_rows = nrow(rows), partial_edge_bins = 0L)
    ))
  }

  level_number <- as.character(choice$level)
  summary_rows <- preprocessed$summary$levels[[level_number]]
  keep <- as.numeric(summary_rows$.bin_start) >= start_seconds &
    summary_rows$.bin_end <= end_seconds
  if (!is.null(selected_series)) {
    keep <- keep & summary_rows$.series_id %in% selected_series
  }
  complete <- summary_rows[keep, , drop = FALSE]
  raw <- data.frame()
  edge_count <- 0L
  if (is.null(preprocessed$raw) &&
      preprocessed$summary$raw_object_bytes > preprocessing_memory_limit()) {
    return(list(
      raw = raw,
      summary = complete,
      resolution = choice,
      diagnostics = list(
        returned_rows = nrow(complete),
        partial_edge_bins = 0L,
        edge_bins_omitted = TRUE
      )
    ))
  }
  if (is.null(preprocessed$raw)) {
    raw <- load_cached_raw_timeseries(preprocessed)
    preprocessed$raw <- raw
  } else {
    raw <- preprocessed$raw
  }
  if (nrow(raw)) {
    base_width <- stats::setNames(
      summary_rows$.base_seconds,
      summary_rows$.series_id
    )
    base_width <- base_width[!duplicated(names(base_width))]
    dt <- data.table::as.data.table(raw)
    data.table::setkeyv(dt, c(".series_id", "datetime"))
    ids <- selected_series %||% unique(dt$.series_id)
    ranges <- data.table::data.table(
      .series_id = ids,
      start = as.POSIXct(start_seconds, origin = "1970-01-01", tz = "UTC"),
      end = as.POSIXct(end_seconds, origin = "1970-01-01", tz = "UTC")
    )
    visible <- dt[ranges, on = .(.series_id, datetime >= start, datetime < end),
                  nomatch = 0L, allow.cartesian = TRUE]
    widths <- unname(base_width[visible$.series_id]) * (2 ^ choice$level)
    bins <- floor(as.numeric(visible$datetime) / widths) * widths
    fully_represented <- bins >= start_seconds & bins + widths <= end_seconds
    edge_rows <- visible[!fully_represented]
    if (nrow(edge_rows)) {
      edge_bins <- bins[!fully_represented]
      edge_rows[, .plot_group := paste(
        .series_id, .gap_id, edge_bins, "edge", sep = "|"
      )]
    }
    raw <- as.data.frame(edge_rows)
    edge_count <- nrow(raw)
  }
  list(
    raw = raw,
    summary = complete,
    resolution = choice,
    diagnostics = list(
      returned_rows = nrow(raw) + nrow(complete),
      partial_edge_bins = edge_count
    )
  )
}
