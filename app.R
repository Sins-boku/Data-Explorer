# =============================================================================
# app.R  -  Rosalia plotting app
#
# Reads the tidy .rds files in data_processed/ and builds every dropdown from
# the data itself, so a new dataset saved there appears without editing this
# file. A dataset only needs: a `variable` column, a `value` column, and
# either `datetime` or `date`.
#
# RUN:  open app.R, click "Run App"   (or  shiny::runApp("app.R"))
#
# SAVE YOUR DATA FIRST (re-run whenever it changes):
#   saveRDS(raw_soilcore, "data_processed/soilcore.rds")
#   saveRDS(met_d,        "data_processed/met_d.rds")
#   saveRDS(sm,           "data_processed/soil_moisture.rds")
#  a test comment
# =============================================================================

#library(shiny)
#library(tidyverse)
#library(plotly)

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(tibble)
library(scales)
library(plotly)

source(file.path("R", "preprocessing.R"))

INPUT_DIRS <- c("data_input", "data")
DATA_DIR <- "data_processed"
ID_CANDIDATES <- PREPROCESSING_ID_COLUMNS

# Soil colour ramp, keyed to the depth ITSELF: 10 cm = black, 20 cm = #3D1F00,
# ... 100 cm = grey. A depth between two entries (SoPhy's 15/30/45) is
# interpolated along the ramp instead of stealing the next colour.
DEPTH_PAL <- c("#000000", "#3D1F00", "#5C3317", "#7B3F00", "#996633",
               "#D2B48C", "#F5DEB3", "#FFFDD0", "#D9D9D9", "#808080")
DEPTH_STEP <- 10   # cm per entry of DEPTH_PAL

depth_cols <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  i <- pmin(pmax(x / DEPTH_STEP, 1), length(DEPTH_PAL))      # 10 cm -> 1st
  out <- rep(NA_character_, length(x))
  ok <- !is.na(i)
  if (any(ok)) {
    m <- grDevices::colorRamp(DEPTH_PAL)((i[ok] - 1) / (length(DEPTH_PAL) - 1))
    out[ok] <- grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
  }
  out
}

TRM_RAMP <- list(D = c("red3", "#FB6A4A"),    # drought   = reds
                 R = c("royalblue", "#6BAED6"))    # reference = blues

`%||%`   <- function(a, b) if (is.null(a)) b else a
`%|""|%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a

# -----------------------------------------------------------------------------
load_dataset <- function(path) {
  d <- readRDS(path)
  if (!all(c("variable", "value") %in% names(d))) return(NULL)
  if (!"datetime" %in% names(d)) {
    if (!"date" %in% names(d)) return(NULL)
    d$datetime <- as.POSIXct(paste(d$date, "12:00:00"), tz = "Europe/Vienna")
  }
  if (!"date" %in% names(d)) d$date <- as.Date(d$datetime)
  if (!"unit" %in% names(d)) d$unit <- NA_character_
  # Keep the known id columns AND anything else that looks like a label, so a
  # new reader's own columns (medium, probe, method, ...) can be used to split,
  # facet and colour without being added to ID_CANDIDATES.
  drop <- c("dataset", "datetime", "date", "variable", "value", "unit",
            "source_file", "timestep_min")
  extra <- setdiff(names(d), c(ID_CANDIDATES, drop))
  extra <- extra[!str_starts(extra, "flag_")]
  extra <- extra[map_lgl(extra, function(v) {
    x <- d[[v]]; n <- n_distinct(x, na.rm = TRUE)
    if (is.character(x) || is.factor(x) || is.logical(x)) return(n > 0 && n <= 100)
    is.numeric(x) && n > 0 && n <= 50 && all(x[!is.na(x)] %% 1 == 0)
  })]
  ids <- c(intersect(ID_CANDIDATES, names(d)), extra)
  d %>%
    mutate(dataset = tools::file_path_sans_ext(basename(path))) %>%
    select(dataset, datetime, date, variable, value, unit, all_of(ids)) %>%
    filter(!is.na(datetime))
}

# include_variable = FALSE is used for dual-isotope plots: there the two
# isotopes are the x and y of ONE point, so the series label must NOT contain
# the variable name or the same sample would end up in two different series
# and the two isotopes could never be paired.
make_series <- function(d, id_cols, include_variable = TRUE) {
  if (!length(id_cols)) {
    return(mutate(d, series = if (include_variable) variable else "all"))
  }
  d <- d %>%
    unite("id_label", all_of(id_cols), sep = " | ", remove = FALSE, na.rm = TRUE)
  if (!include_variable) {
    return(mutate(d, series = if_else(id_label == "", "all", id_label)))
  }
  d %>% mutate(series = if_else(id_label == "", variable,
                                paste0(variable, " (", id_label, ")")))
}

# colours keyed on the series present, following the chosen scheme
build_colours <- function(d, scheme, manual = NULL) {
  if (!"ckey" %in% names(d)) d$ckey <- d$series
  s <- sort(unique(d$ckey))
  
  if (scheme == "depth" && "depth_bottom" %in% names(d)) {
    key <- d %>% distinct(ckey, depth_bottom)
    out <- depth_cols(key$depth_bottom)
    names(out) <- key$ckey
    # rows with no depth (stem probes, precipitation, ...) are not part of the
    # soil ramp, so they get their own colours instead of all sharing one grey
    if (any(is.na(out)))
      out[is.na(out)] <- scales::hue_pal()(sum(is.na(out)))
    out <- out[s]
    names(out) <- s
    return(out)
  }
  
  if (scheme == "trm" && "trm" %in% names(d)) {
    key <- d %>% distinct(ckey, trm) %>% arrange(ckey)
    out <- setNames(rep("grey50", nrow(key)), key$ckey)
    for (tr in names(TRM_RAMP)) {
      idx <- which(key$trm == tr)
      if (!length(idx)) next
      out[idx] <- colorRampPalette(TRM_RAMP[[tr]])(length(idx))
    }
    return(out[s])
  }
  
  if (scheme == "date") {
    # colour is keyed on the sampling date itself, not the series, so the
    # legend lists dates (one entry per campaign) rather than series names
    lev <- sort(unique(d$ckey))
    return(setNames(rainbow(length(lev), end = 0.85), lev))
  }
  
  if (scheme == "manual" && !is.null(manual)) return(manual[s])
  
  setNames(scales::hue_pal()(length(s)), s)
}

# -----------------------------------------------------------------------------
source(file.path("Env-Plot", "R", "renderer.R"))

ui <- fluidPage(
  titlePanel("Rosalia data explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # options scroll on their own, the plot stays put
      style = "height: calc(100vh - 90px); overflow-y: auto; position: sticky; top: 10px;",
      selectInput("datasets", "Datasets", choices = NULL, multiple = TRUE),
      uiOutput("file_status"),
      uiOutput("schema_ui"),
      numericInput(
        "gap_factor",
        "Gap threshold (cadence multiples)",
        value = 5,
        min = 1.1,
        max = 100,
        step = 0.5
      ),
      
      radioButtons("plotmode", "Plot type",
                   c("Time series"   = "ts",
                     "Depth profile" = "depth",
                     "Dual isotope"  = "dual"), selected = "ts"),
      
      uiOutput("var_ui"),
      uiOutput("dual_ui"),
      uiOutput("id_ui"),
      uiOutput("id_keep_ui"),
      uiOutput("facet_by_ui"),
      uiOutput("depth_ui"),
      uiOutput("date_ui"),
      uiOutput("dates_pick_ui"),
      
      hr(),
      uiOutput("axis_mode_ui"),
      uiOutput("sec_ui"),
      
      hr(),
      checkboxInput("interactive", "Interactive time series (zoom, pan, hover)", TRUE),
      radioButtons("geom", "Draw as",
                   c("Line" = "line", "Points" = "point",
                     "Line + points" = "both", "Bars" = "col"),
                   selected = "line"),
      sliderInput("alpha", "Opacity", 0.1, 1, 0.9, step = 0.05),
      sliderInput("lwd", "Line width", 0.1, 2, 0.5, step = 0.1),
      sliderInput("psize", "Point size", 0.5, 6, 2, step = 0.5),
      
      hr(),
      selectInput("scheme", "Colour by",
                  c("Automatic" = "auto", "Treatment (D red / R blue)" = "trm",
                    "Depth (soil ramp)" = "depth", "Date (rainbow)" = "date",
                    "Pick by hand" = "manual"), selected = "auto"),
      uiOutput("colour_ui"),
      
      hr(),
      textInput("title", "Title", ""),
      textInput("xlab", "X label (blank = auto)", ""),
      textInput("ylab", "Y label / left axis", ""),
      uiOutput("ylab2_ui"),
      uiOutput("facet_lab_ui"),
      
      hr(),
      uiOutput("xlim_ui"),
      fluidRow(
        column(6, numericInput("ymin", "Y min", NA)),
        column(6, numericInput("ymax", "Y max", NA))
      ),
      helpText("Axis range: blank = automatic. Applies to every panel when facetted."),
      actionLink("lim_reset", "reset axis range"),
      
      hr(),
      numericInput("w", "Width (cm)", 25, 4, 60, 0.5),
      numericInput("h", "Height (cm)", 15, 4, 60, 0.5),
      numericInput("dpi", "DPI", 300, 72, 900, 10),
      numericInput("basesize", "Base font size (pt)", 12, 6, 24, 1),
      numericInput("legcols", "Legend columns (0 = auto)", 0, 0, 12, 1),
      downloadButton("dl_png", "PNG"),
      downloadButton("dl_pdf", "PDF")
    ),
    
    mainPanel(
      width = 9,
      style = "position: sticky; top: 10px;",
      uiOutput("plot_area"),
      helpText("Drag a box then double-click to zoom. Double-click on empty space to reset."),
      hr(),
      verbatimTextOutput("info")
    )
  )
)

# -----------------------------------------------------------------------------
server <- function(input, output, session) {
  file_discovery_error <- reactiveVal(NULL)
  
  files <- reactive({
    invalidateLater(30000, session)
    tryCatch(
      {
        result <- discover_data_files(INPUT_DIRS)
        file_discovery_error(NULL)
        result
      },
      error = function(error) {
        file_discovery_error(conditionMessage(error))
        data.frame(
          path = character(),
          dataset = character(),
          format = character(),
          size_bytes = numeric(),
          modified = as.POSIXct(character(), tz = "UTC")
        )
      }
    )
  })
  
  observe({
    source_files <- files()
    choices <- source_files$path
    labels <- paste0(source_files$dataset, " (", basename(dirname(choices)), ")")
    names(choices) <- make.unique(labels)
    updateSelectInput(
      session,
      "datasets",
      choices = choices,
      selected = intersect(isolate(input$datasets), unname(choices))
    )
  })

  output$file_status <- renderUI({
    if (!is.null(file_discovery_error())) {
      return(div(
        class = "alert alert-danger",
        file_discovery_error()
      ))
    }
    if (!nrow(files())) {
      return(helpText(
        "No CSV or RDS datasets found in data_input/ or data/."
      ))
    }
    NULL
  })

  output$schema_ui <- renderUI({
    paths <- input$datasets
    if (is.null(paths) || !length(paths)) return(NULL)
    tagList(lapply(seq_along(paths), function(index) {
      path <- paths[index]
      columns <- read_timeseries_table(path, header_only = TRUE)
      normalized <- tolower(gsub("[^a-z0-9]+", "_", columns))
      default_for <- function(aliases) {
        matches <- columns[normalized %in% aliases]
        if (length(matches) == 1L) matches else ""
      }
      group_choices <- intersect(PREPROCESSING_ID_COLUMNS, columns)
      id <- as.character(index)
      tagList(
        h5(paste("Column mapping:", basename(path))),
        selectInput(
          paste0("schema_datetime_", id),
          "Timestamp column",
          choices = c("Detect automatically" = "", stats::setNames(columns, columns)),
          selected = default_for(c(
            "datetime", "timestamp", "time", "date_time", "date", "sample_time"
          ))
        ),
        selectInput(
          paste0("schema_value_", id),
          "Measurement/value column",
          choices = c("Detect automatically" = "", stats::setNames(columns, columns)),
          selected = default_for(c(
            "value", "measurement", "reading", "observation", "result"
          ))
        ),
        selectInput(
          paste0("schema_variable_", id),
          "Variable column (optional)",
          choices = c("Use file name / infer" = "", stats::setNames(columns, columns)),
          selected = default_for(c(
            "variable", "parameter", "measurement_type", "measurand"
          ))
        ),
        selectInput(
          paste0("schema_unit_", id),
          "Unit column (optional)",
          choices = c("No unit / infer" = "", stats::setNames(columns, columns)),
          selected = default_for(c("unit", "units"))
        ),
        checkboxGroupInput(
          paste0("schema_groups_", id),
          "Series/group columns",
          choices = columns,
          selected = group_choices,
          inline = TRUE
        ),
        textInput(paste0("schema_timezone_", id), "Timezone", "UTC")
      )
    }))
  })

  schema_mappings <- reactive({
    paths <- input$datasets %||% character()
    lapply(seq_along(paths), function(index) {
      id <- as.character(index)
      list(
        datetime = input[[paste0("schema_datetime_", id)]],
        value = input[[paste0("schema_value_", id)]],
        variable = input[[paste0("schema_variable_", id)]],
        unit = input[[paste0("schema_unit_", id)]],
        group_cols = input[[paste0("schema_groups_", id)]],
        timezone = input[[paste0("schema_timezone_", id)]] %||% "UTC"
      )
    })
  })

  preprocessed <- reactive({
    req(input$datasets)
    paths <- input$datasets
    mappings <- schema_mappings()
    lapply(seq_along(paths), function(index) {
      tryCatch(
        load_preprocessed_timeseries(
          paths[index],
          mappings[[index]],
          cache_dir = file.path(DATA_DIR, ".cache"),
          gap_factor = input$gap_factor %||% 5
        ),
        error = function(error) {
          showNotification(conditionMessage(error), type = "error", duration = NULL)
          stop(error)
        }
      )
    })
  })
  
  # A depth profile is almost always read as "how did the profile change
  # between campaigns", so colouring by sampling date is the sensible default
  # there. Switching plot type flips it; the user can still override.
  observeEvent(input$plotmode, {
    updateSelectInput(session, "scheme",
                      selected = if (identical(input$plotmode, "depth")) "date" else "auto")
  }, ignoreInit = TRUE)
  observeEvent(input$datasets, {
    zoom(NULL)
  }, ignoreInit = TRUE)
  observeEvent(input$dates, {
    zoom(NULL)
  }, ignoreInit = TRUE)
  
  observeEvent(input$facet_by, {                                 
    if (isTRUE(input$facet_by %in% c("trm", "tree")))
      updateSelectInput(session, "scheme", selected = "trm")
  }, ignoreInit = TRUE)
  # Adding a dataset from another period used to leave the old, narrower date
  # range in place; the new variable then had zero rows and facet_wrap dropped
  # its panel without a word - "the second figure never shows up".
  dat <- reactive({
    objects <- preprocessed()
    maps <- lapply(objects, function(object) {
      series_map <- object$summary$series_map
      ranges <- object$summary$time_ranges
      metadata <- dplyr::left_join(series_map, ranges, by = ".series_id")
      metadata$.start_date <- as.Date(metadata$.start, tz = "UTC")
      metadata$.end_date <- as.Date(metadata$.end, tz = "UTC")
      for (timezone in unique(metadata$time_zone)) {
        rows <- which(metadata$time_zone == timezone)
        metadata$.start_date[rows] <- as.Date(
          metadata$.start[rows],
          tz = timezone
        )
        metadata$.end_date[rows] <- as.Date(
          metadata$.end[rows],
          tz = timezone
        )
      }
      metadata$date <- metadata$.start_date
      metadata
    })
    dplyr::bind_rows(maps)
  })

  data_date_range <- reactive({
    metadata <- dat()
    range(c(metadata$.start_date, metadata$.end_date), na.rm = TRUE)
  })

  observeEvent(dat(), {
    date_range <- data_date_range()
    if (all(is.finite(date_range))) {
      updateDateRangeInput(
        session,
        "dates",
        start = date_range[1],
        end = date_range[2],
        min = date_range[1],
        max = date_range[2]
      )
    }
  }, ignoreInit = TRUE)

  # The widget can be momentarily absent while renderUI is updating.
  date_window <- reactive({
    rng <- data_date_range()
    if (is.null(input$dates) || !all(is.finite(rng))) rng else as.Date(input$dates)
  })

  sampling_dates <- reactive({
    bounds <- date_window()
    dates <- unlist(lapply(preprocessed(), function(object) {
      timezone <- object$mapping$timezone %||% "UTC"
      as.Date(object$summary$levels[["0"]]$.bin_start, tz = timezone)
    }), use.names = FALSE)
    sort(unique(dates[dates >= bounds[1] & dates <= bounds[2]]))
  })
  
  # --- variable pickers ------------------------------------------------------
  output$var_ui <- renderUI({
    if (identical(input$plotmode, "dual")) return(NULL)
    v <- sort(unique(dat()$variable))
    selectInput("vars", "Variables", choices = v, multiple = TRUE,
                selected = head(v, 1))
  })
  
  output$dual_ui <- renderUI({
    if (!identical(input$plotmode, "dual")) return(NULL)
    v <- sort(unique(dat()$variable))
    tagList(
      selectInput("xvar", "X variable", choices = v,
                  selected = if ("d18O" %in% v) "d18O" else v[1]),
      selectInput("yvar", "Y variable", choices = v,
                  selected = if ("d2H" %in% v) "d2H" else tail(v, 1)),
      checkboxInput("gmwl", "Show GMWL (d2H = 8*d18O + 10)", TRUE),
      checkboxInput("lmwl", "Show LMWL", TRUE),
      fluidRow(
        column(6, numericInput("lmwl_slope", "LMWL slope", 7.6603, step = 0.001)),
        column(6, numericInput("lmwl_int", "LMWL intercept", 4.4649, step = 0.1))
      ),
      hr(),
      radioButtons("fitline", "Regression line",
                   c("None" = "none",
                     "One line for all points" = "all",
                     "One line per series" = "series"),
                   selected = "none"),
      checkboxInput("fit_eq", "Show slope / intercept / R2", TRUE)
    )
  })
  
  output$id_ui <- renderUI({
    ids <- intersect(ID_CANDIDATES, names(dat()))
    ids <- ids[map_lgl(ids, ~ n_distinct(dat()[[.x]], na.rm = TRUE) > 1)]
    selectInput("ids", "Split series by", choices = ids, multiple = TRUE,
                selected = intersect(c("logger", "trm"), ids))
  })
  
  # one tick box per level of every "Split series by" column, so single trees,
  # treatments or loggers can be dropped from the plot without re-filtering
  output$id_keep_ui <- renderUI({
    d <- dat(); ids <- intersect(input$ids %||% character(), names(d))
    if (!length(ids)) return(NULL)
    tagList(lapply(ids, function(v) {
      lv <- as.character(sort(unique(d[[v]])))
      if (identical(v, "tree")) lv <- c(intersect(c("6", "4", "7", "10"), lv),
                                        setdiff(lv, c("6", "4", "7", "10")))
      keep <- isolate(input[[paste0("keep_", v)]])
      checkboxGroupInput(paste0("keep_", v), paste("Show", v), inline = TRUE,
                         choices = lv,
                         selected = if (is.null(keep)) lv else intersect(lv, keep))
    }))
  })
  
  output$facet_by_ui <- renderUI({
    if (identical(input$plotmode, "ts")) return(NULL)
    ids <- intersect(ID_CANDIDATES, names(dat()))
    ids <- ids[map_lgl(ids, ~ n_distinct(dat()[[.x]], na.rm = TRUE) > 1)]
    selectInput("facet_by", "Facet panels by",
                choices = c("None" = "", "variable", ids), selected = "")
  })
  
  # --- show / hide individual depths ----------------------------------------
  output$depth_ui <- renderUI({
    d <- dat()
    if (!"depth_bottom" %in% names(d)) return(NULL)
    lev <- sort(unique(d$depth_bottom[!is.na(d$depth_bottom)]))
    if (!length(lev)) return(NULL)
    checkboxGroupInput("depths", "Depths to show (depth_bottom)",
                       choices = lev, selected = lev, inline = TRUE)
  })
  
  output$date_ui <- renderUI({
    r <- data_date_range()
    dateRangeInput("dates", "Date range", start = r[1], end = r[2],
                   min = r[1], max = r[2])
  })
  
  # Tick individual sampling dates on and off, the same way depths work.
  # Only offered when there is a manageable number of distinct dates, so it
  # does not try to list every timestamp of a 10-minute logger series.
  output$dates_pick_ui <- renderUI({
    req(input$dates)
    lev <- sampling_dates()
    if (!length(lev) || length(lev) > 60) return(NULL)
    lev_chr <- format(lev, "%Y-%m-%d")
    tagList(
      checkboxGroupInput("pick_dates", "Sampling dates to show",
                         choices = lev_chr, selected = lev_chr, inline = TRUE),
      actionLink("dates_all", "all"), " / ", actionLink("dates_none", "none")
    )
  })
  
  observeEvent(input$dates_all, {
    lev <- format(sampling_dates(), "%Y-%m-%d")
    updateCheckboxGroupInput(session, "pick_dates", selected = lev)
  })
  observeEvent(input$dates_none, {
    updateCheckboxGroupInput(session, "pick_dates", selected = character(0))
  })
  
  # x is a datetime in a time series, where the existing "Date range" already
  # sets the span, so numeric x limits are only offered for the other modes
  output$xlim_ui <- renderUI({
    if (identical(input$plotmode, "ts")) return(NULL)
    tagList(
      fluidRow(
        column(6, numericInput("xmin", "X min", NA)),
        column(6, numericInput("xmax", "X max", NA))
      ),
      # A facetted depth profile can need a different x range on each side of
      # the grid (e.g. the two reference trees vs the two drought trees), so
      # the right-hand panel column gets its own pair. Blank = same as left.
      if (identical(input$plotmode, "depth"))
        fluidRow(
          column(6, numericInput("xmin_r", "X min (right col.)", NA)),
          column(6, numericInput("xmax_r", "X max (right col.)", NA))
        )
    )
  })
  
  observeEvent(input$lim_reset, {
    for (i in c("xmin", "xmax", "xmin_r", "xmax_r", "ymin", "ymax"))
      updateNumericInput(session, i, value = NA)
  })
  
  output$axis_mode_ui <- renderUI({
    if (!identical(input$plotmode, "ts")) return(NULL)
    radioButtons("mode", "Y axes",
                 c("One axis" = "single", "Facet per variable" = "facet",
                   "Secondary axis" = "sec"), selected = "facet")
  })
  
  output$sec_ui <- renderUI({
    if (!identical(input$plotmode, "ts") || !identical(input$mode, "sec")) return(NULL)
    req(input$vars)
    selectInput("sec_var", "Variable on the RIGHT axis",
                choices = input$vars, selected = tail(input$vars, 1))
  })
  
  output$ylab2_ui <- renderUI({
    if (!identical(input$plotmode, "ts") || !identical(input$mode, "sec")) return(NULL)
    textInput("ylab2", "Right axis label", "")
  })
  
  # per-variable labels shown as facet strip text
  output$facet_lab_ui <- renderUI({
    if (!identical(input$plotmode, "ts") || !identical(input$mode, "facet")) return(NULL)
    req(input$vars)
    tagList(
      radioButtons("strippos", "Facet label position",
                   c("Above each panel" = "top",
                     "Left, as per-panel y-axis label" = "left"),
                   selected = "left"),
      helpText("Choose 'Left' to give each variable its own y-axis label (e.g. permille on top, % below):"),
      lapply(input$vars, function(v)
        textInput(paste0("flab_", make.names(v)), v, value = v))
    )
  })
  
  selected_metadata <- reactive({
    metadata <- dat()
    bounds <- date_window()
    d <- metadata %>%
      filter(.end_date >= bounds[1], .start_date <= bounds[2])

    if ("tree" %in% names(d)) {
      tree_levels <- c(6, 4, 7, 10)
      d <- d %>% mutate(tree = factor(tree, levels = c(
        tree_levels, setdiff(sort(unique(tree)), tree_levels)
      )))
    }
    for (column in intersect(input$ids %||% character(), names(d))) {
      keep <- input[[paste0("keep_", column)]]
      if (!is.null(keep)) {
        d <- d[as.character(d[[column]]) %in% keep, , drop = FALSE]
      }
    }
    d <- droplevels(d)

    if ("depth_bottom" %in% names(d) && !is.null(input$depths)) {
      d <- d %>% filter(
        is.na(depth_bottom) |
          as.character(depth_bottom) %in% as.character(input$depths)
      )
    }
    if (identical(input$plotmode, "dual")) {
      req(input$xvar, input$yvar)
      d <- d %>% filter(variable %in% c(input$xvar, input$yvar))
    } else {
      req(input$vars)
      d <- d %>% filter(variable %in% input$vars)
    }
    d
  })

  plot_viewport <- reactive({
    bounds <- date_window()
    current_zoom <- zoom()
    width_id <- if (isTRUE(input$interactive) &&
                    identical(input$plotmode, "ts")) "iplot" else "plot"
    plot_width <- session$clientData[[paste0("output_", width_id, "_width")]]
    if (is.null(plot_width) || !is.finite(plot_width) || plot_width <= 0) {
      plot_width <- 1000
    }
    list(bounds = bounds, zoom = current_zoom, width = plot_width)
  })

  plot_dat <- reactive({
    metadata <- selected_metadata()
    validate(need(nrow(metadata) > 0, "No data for that selection."))
    objects <- preprocessed()
    selected_ids <- unique(metadata$.series_id)
    viewport <- plot_viewport()
    if (identical(input$plotmode, "ts")) {
      pieces <- lapply(objects, function(object) {
        object_ids <- intersect(
          selected_ids,
          object$summary$series_map$.series_id
        )
        if (!length(object_ids)) return(NULL)
        timezone <- object$mapping$timezone %||% "UTC"
        start <- as.POSIXct(
          paste(viewport$bounds[1], "00:00:00"),
          tz = timezone
        )
        end <- as.POSIXct(
          paste(viewport$bounds[2] + 1, "00:00:00"),
          tz = timezone
        )
        if (!is.null(viewport$zoom)) {
          zoom_x <- as.POSIXct(
            viewport$zoom$x,
            origin = "1970-01-01",
            tz = timezone
          )
          if (length(zoom_x) == 2L && all(is.finite(zoom_x))) {
            start <- max(start, zoom_x[1])
            end <- min(end, zoom_x[2])
          }
        }
        if (start >= end) return(NULL)
        query <- query_timeseries_viewport(
          object,
          start,
          end,
          plot_width_px = viewport$width,
          selected_series = object_ids
        )
        raw <- query$raw
        if (nrow(raw)) {
          raw$.is_summary <- FALSE
          if (!".plot_group" %in% names(raw)) {
            raw$.plot_group <- paste(
              raw$.series_id, raw$.gap_id, "raw", sep = "|"
            )
          }
        }
        summaries <- query$summary
        if (nrow(summaries)) {
          summaries <- dplyr::left_join(
            summaries,
            object$summary$series_map,
            by = ".series_id"
          )
          summaries$datetime <- summaries$.bin_start +
            summaries$.width_seconds / 2
          summaries$date <- as.Date(
            summaries$datetime,
            tz = timezone
          )
          summaries$value <- summaries$.mean
          summaries$.is_summary <- TRUE
          summaries$.plot_group <- paste(
            summaries$.series_id,
            summaries$.gap_id,
            "summary",
            sep = "|"
          )
        }
        bind_rows(raw, summaries)
      })
      d <- bind_rows(pieces)
    } else {
      d <- bind_rows(lapply(objects, function(object) {
        object_ids <- intersect(
          selected_ids,
          object$summary$series_map$.series_id
        )
        if (!length(object_ids)) return(NULL)
        raw <- load_cached_raw_timeseries(object)
        timezone <- object$mapping$timezone %||% "UTC"
        start <- as.POSIXct(
          paste(viewport$bounds[1], "00:00:00"),
          tz = timezone
        )
        end <- as.POSIXct(
          paste(viewport$bounds[2] + 1, "00:00:00"),
          tz = timezone
        )
        raw %>%
          filter(
            .series_id %in% object_ids,
            datetime >= start,
            datetime < end
          )
      }))
      d$.is_summary <- FALSE
      d$.plot_group <- paste(d$.series_id, d$.gap_id, "raw", sep = "|")
    }

    d <- d %>% filter(!is.na(value))
    if (!is.null(input$pick_dates)) {
      d <- d %>% filter(format(date, "%Y-%m-%d") %in% input$pick_dates)
    }
    if (!"series" %in% names(d)) {
      d <- make_series(
        d,
        input$ids %||% character(),
        include_variable = !identical(input$plotmode, "dual")
      )
    }
    if (identical(input$scheme, "date")) {
      d <- d %>% mutate(ckey = format(date, "%Y-%m-%d"))
    } else if (identical(input$scheme, "trm") &&
               identical(input$plotmode, "depth")) {
      d <- d %>% mutate(ckey = paste0(series, " | ", format(date, "%Y-%m-%d")))
    } else {
      d <- d %>% mutate(ckey = series)
    }
    d
  })
  
  output$colour_ui <- renderUI({
    if (!identical(input$scheme, "manual")) return(NULL)
    s <- sort(unique(plot_dat()$series))
    validate(need(length(s) > 0, "No data for that selection."))
    if (length(s) > 15)
      return(helpText(paste(length(s), "series - too many to pick by hand.")))
    pal <- scales::hue_pal()(length(s))
    tagList(lapply(seq_along(s), function(i)
      textInput(paste0("col_", i), s[i], value = pal[i])))
  })
  
  manual_cols <- reactive({
    s <- sort(unique(plot_dat()$series))
    pal <- scales::hue_pal()(length(s))
    out <- map_chr(seq_along(s), function(i) {
      v <- input[[paste0("col_", i)]]
      if (is.null(v) || !nzchar(v)) pal[i] else v
    })
    ok <- map_lgl(out, ~ tryCatch({ grDevices::col2rgb(.x); TRUE },
                                  error = function(e) FALSE))
    out[!ok] <- pal[!ok]
    setNames(out, s)
  })
  
  zoom <- reactiveVal(NULL)
  observeEvent(input$dbl, {
    b <- input$brush
    zoom(if (!is.null(b)) list(x = c(b$xmin, b$xmax), y = c(b$ymin, b$ymax)) else NULL)
  })

  relayout <- reactive(
    plotly::event_data(
      "plotly_relayout",
      source = "timeseries",
      priority = "event"
    )
  )
  observeEvent(relayout(), {
    if (!isTRUE(input$interactive) || !identical(input$plotmode, "ts")) return()
    event <- relayout()
    if (is.null(event)) return()
    values <- as.list(event)
    if (isTRUE(values[["xaxis.autorange"]])) {
      if (!is.null(zoom())) zoom(NULL)
      return()
    }
    lower <- values[["xaxis.range[0]"]]
    upper <- values[["xaxis.range[1]"]]
    if (is.null(lower) || is.null(upper)) {
      range_values <- values[["xaxis.range"]]
      if (length(range_values) == 2L) {
        lower <- range_values[[1L]]
        upper <- range_values[[2L]]
      }
    }
    if (is.null(lower) || is.null(upper)) return()
    to_epoch <- function(value) {
      if (is.numeric(value)) {
        seconds <- if (abs(value) > 1e11) value / 1000 else value
        return(seconds)
      }
      as.numeric(as.POSIXct(as.character(value), tz = "UTC"))
    }
    range_x <- c(to_epoch(lower), to_epoch(upper))
    if (any(!is.finite(range_x)) || range_x[1] >= range_x[2]) return()
    current <- zoom()
    if (is.null(current) ||
        !isTRUE(all.equal(as.numeric(current$x), range_x, tolerance = 0.001))) {
      zoom(list(x = range_x, y = current$y %||% NULL))
    }
  })
  
  # --- the plot --------------------------------------------------------------
  build_plot <- reactive({
    data <- plot_dat()
    validate(need(nrow(data) > 0, "No data for that selection."))
    env_plot_build_plot(data, input, manual_cols(), zoom())
  })

  # Static output keeps the brush/double-click workflow; Plotly additionally
  # reports zoom and pan ranges so the server can query a finer representation.
  output$plot_area <- renderUI({
    if (isTRUE(input$interactive) && identical(input$plotmode, "ts")) {
      if (!requireNamespace("plotly", quietly = TRUE))
        return(helpText("Interactive mode needs plotly: install.packages('plotly')"))
      plotly::plotlyOutput("iplot", height = "640px")
    } else {
      plotOutput("plot", height = "640px", dblclick = "dbl",
                 brush = brushOpts("brush", resetOnNew = TRUE))
    }
  })
  
  output$plot <- renderPlot(build_plot())
  
  output$iplot <- plotly::renderPlotly({
    p <- plotly::ggplotly(
      build_plot(),
      tooltip = c("x", "y", "colour"),
      source = "timeseries"
    )
    p <- plotly::plotly_build(p)
    # the colour aesthetic is called ckey internally; the hover box should not
    # say so
    p$x$data <- lapply(p$x$data, function(tr) {
      if (!is.null(tr$text))
        tr$text <- gsub("ckey:", "series:", tr$text, fixed = TRUE)
      tr
    })
    plotly::layout(p, hovermode = "closest",
                   legend = list(orientation = "h", y = -0.15),
                   uirevision = "data-explorer")
  })
  
  output$info <- renderText({
    d <- plot_dat()
    paste0(nrow(d), " points | ", n_distinct(d$series), " series | ",
           paste(sort(unique(d$dataset)), collapse = ", "))
  })
  
  output$dl_png <- downloadHandler(
    filename = function() paste0("rosalia_", Sys.Date(), ".png"),
    content  = function(file) ggsave(file, build_plot(), width = input$w,
                                     height = input$h, units = "cm",
                                     dpi = input$dpi)
  )
  output$dl_pdf <- downloadHandler(
    filename = function() paste0("rosalia_", Sys.Date(), ".pdf"),
    content  = function(file) ggsave(file, build_plot(), width = input$w,
                                     height = input$h, units = "cm",
                                     device = cairo_pdf)
  )
}

shinyApp(ui, server)
