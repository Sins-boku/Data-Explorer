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

DATA_DIR <- "data_processed"

ID_CANDIDATES <- c("logger", "position", "source", "trm", "tree", "core_id",
                   "probe_set", "depth_bottom", "depth_top", "station")

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
ui <- fluidPage(
  titlePanel("Rosalia data explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      # options scroll on their own, the plot stays put
      style = "height: calc(100vh - 90px); overflow-y: auto; position: sticky; top: 10px;",
      selectInput("datasets", "Datasets", choices = NULL, multiple = TRUE),
      
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
      checkboxInput("interactive", "Interactive time series (hover for values)", FALSE),
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
  
  files <- reactive({
    invalidateLater(30000, session)
    list.files(DATA_DIR, pattern = "\\.rds$", full.names = TRUE)
  })
  
  observe({
    ch <- setNames(files(), tools::file_path_sans_ext(basename(files())))
    updateSelectInput(session, "datasets", choices = ch,
                      selected = isolate(input$datasets))
  })
  
  # A depth profile is almost always read as "how did the profile change
  # between campaigns", so colouring by sampling date is the sensible default
  # there. Switching plot type flips it; the user can still override.
  observeEvent(input$plotmode, {
    updateSelectInput(session, "scheme",
                      selected = if (identical(input$plotmode, "depth")) "date" else "auto")
  }, ignoreInit = TRUE)
  
  observeEvent(input$facet_by, {                                 
    if (isTRUE(input$facet_by %in% c("trm", "tree")))
      updateSelectInput(session, "scheme", selected = "trm")
  }, ignoreInit = TRUE)
  # Adding a dataset from another period used to leave the old, narrower date
  # range in place; the new variable then had zero rows and facet_wrap dropped
  # its panel without a word - "the second figure never shows up".
  observeEvent(input$datasets, {
    r <- range(dat()$date, na.rm = TRUE)
    if (all(is.finite(r)))
      updateDateRangeInput(session, "dates", start = r[1], end = r[2],
                           min = r[1], max = r[2])
  })
  
  # The widget can be momentarily absent (it is rendered by renderUI), so the
  # date window falls back to everything that is loaded.
  date_window <- reactive({
    rng <- suppressWarnings(range(dat()$date, na.rm = TRUE))
    if (is.null(input$dates) || !all(is.finite(rng))) rng else as.Date(input$dates)
  })
  
  dat <- reactive({
    req(input$datasets)
    out <- compact(map(input$datasets, possibly(load_dataset, NULL)))
    validate(need(length(out) > 0,
                  "None of those files are in long format (need variable + value)."))
    bind_rows(out)
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
    r <- range(dat()$date, na.rm = TRUE)
    dateRangeInput("dates", "Date range", start = r[1], end = r[2],
                   min = r[1], max = r[2])
  })
  
  # Tick individual sampling dates on and off, the same way depths work.
  # Only offered when there is a manageable number of distinct dates, so it
  # does not try to list every timestamp of a 10-minute logger series.
  output$dates_pick_ui <- renderUI({
    req(input$dates)
    dd <- dat() %>%
      filter(date >= date_window()[1], date <= date_window()[2])
    lev <- sort(unique(dd$date))
    if (!length(lev) || length(lev) > 60) return(NULL)
    lev_chr <- format(lev, "%Y-%m-%d")
    tagList(
      checkboxGroupInput("pick_dates", "Sampling dates to show",
                         choices = lev_chr, selected = lev_chr, inline = TRUE),
      actionLink("dates_all", "all"), " / ", actionLink("dates_none", "none")
    )
  })
  
  observeEvent(input$dates_all, {
    dd <- dat() %>% filter(date >= date_window()[1], date <= date_window()[2])
    lev <- format(sort(unique(dd$date)), "%Y-%m-%d")
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
  
  # --- filtered data ---------------------------------------------------------
  plot_dat <- reactive({
    d <- dat() %>%
      filter(date >= date_window()[1], date <= date_window()[2], !is.na(value))
    
    # panel order: 6 | 4 top, 7 | 10 bottom
    if ("tree" %in% names(d)) {
      lv <- c(6, 4, 7, 10)
      d <- d %>% mutate(tree = factor(tree, levels = c(
        lv, setdiff(sort(unique(tree)), lv))))
    }
    
    # a NULL box means the UI has not been drawn yet -> keep everything;
    # an empty box means the user unticked every level
    for (v in intersect(input$ids %||% character(), names(d))) {
      keep <- input[[paste0("keep_", v)]]
      if (!is.null(keep)) d <- d[as.character(d[[v]]) %in% keep, ]
    }
    d <- droplevels(d)
    
    if (!is.null(input$pick_dates)) {
      d <- d %>% filter(format(date, "%Y-%m-%d") %in% input$pick_dates)
    }
    
    if ("depth_bottom" %in% names(d) && !is.null(input$depths)) {
      d <- d %>% filter(is.na(depth_bottom) |
                          as.character(depth_bottom) %in% as.character(input$depths))
    }
    
    if (identical(input$plotmode, "dual")) {
      req(input$xvar, input$yvar)
      d <- d %>% filter(variable %in% c(input$xvar, input$yvar))
    } else {
      req(input$vars)
      d <- d %>% filter(variable %in% input$vars)
    }
    d <- make_series(d, input$ids %||% character(0),
                     include_variable = !identical(input$plotmode, "dual"))
    # `series` groups the lines; `ckey` drives the colour. Keeping them apart
    # is what lets one line per core be coloured by sampling date.
    if (identical(input$scheme, "date")) {
      d <- d %>% mutate(ckey = format(date, "%Y-%m-%d"))
    } else if (identical(input$scheme, "trm") &&
               identical(input$plotmode, "depth")) {
      # one legend entry per treatment AND sampling date, so each campaign
      # gets its own step of the D-red / R-blue ramp and is named in the key
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
  
  # --- the plot --------------------------------------------------------------
  build_plot <- reactive({
    d <- plot_dat()
    validate(need(nrow(d) > 0, "No data for that selection."))
    
    miss <- setdiff(input$vars %||% character(), unique(d$variable))
    if (length(miss) && !identical(input$plotmode, "dual"))
      showNotification(
        paste0("No data for ", paste(miss, collapse = ", "),
               " with the current filters (date range, depths, dates ticked)",
               " - that panel is missing."),
        type = "warning", duration = 10)
    
    # a blank box stays NA, which every limit argument below reads as "auto"
    nn <- function(x) if (is.null(x) || !is.finite(x)) NA_real_ else as.numeric(x)
    per_col_x <- FALSE   # TRUE once the x range is set per panel column
    
    cols <- build_colours(d, input$scheme,
                          if (identical(input$scheme, "manual")) manual_cols() else NULL)
    
    add_geoms <- function(p) {
      p + switch(input$geom,
                 line  = geom_line(alpha = input$alpha, linewidth = input$lwd),
                 point = geom_point(alpha = input$alpha, size = input$psize),
                 both  = list(geom_line(alpha = input$alpha, linewidth = input$lwd),
                              geom_point(alpha = input$alpha, size = input$psize)),
                 col   = geom_col(aes(fill = series), alpha = input$alpha,
                                  position = "identity", colour = NA))
    }
    
    # ---------------- dual isotope ----------------
    if (identical(input$plotmode, "dual")) {
      req(input$xvar, input$yvar)
      idc <- c("dataset", "series", "ckey", "date", intersect(ID_CANDIDATES, names(d)))
      w <- d %>%
        select(all_of(idc), variable, value) %>%
        distinct(across(all_of(c(idc, "variable"))), .keep_all = TRUE) %>%
        pivot_wider(names_from = variable, values_from = value)
      validate(need(all(c(input$xvar, input$yvar) %in% names(w)) &&
                      sum(!is.na(w[[input$xvar]]) & !is.na(w[[input$yvar]])) > 0,
                    "No samples have both variables. Try adding core_id under 'Split series by'."))
      
      p <- ggplot(w, aes(.data[[input$xvar]], .data[[input$yvar]], colour = ckey))
      if (isTRUE(input$gmwl))
        p <- p + geom_abline(slope = 8, intercept = 10,
                             linetype = "dashed", colour = "grey40")
      if (isTRUE(input$lmwl) && !is.null(input$lmwl_slope) && !is.null(input$lmwl_int))
        p <- p + geom_abline(slope = input$lmwl_slope, intercept = input$lmwl_int,
                             linetype = "dotdash", colour = "steelblue")
      p <- p + geom_point(alpha = input$alpha, size = input$psize)
      
      # --- ordinary least squares fit through the plotted points -------------
      fit_txt <- NULL
      if (!identical(input$fitline, "none")) {
        fml <- as.formula(paste0("`", input$yvar, "` ~ `", input$xvar, "`"))
        
        fit_one <- function(dd, lab) {
          dd <- dd[!is.na(dd[[input$xvar]]) & !is.na(dd[[input$yvar]]), ]
          if (nrow(dd) < 3) return(NULL)
          m  <- lm(fml, data = dd)
          cf <- coef(m)
          tibble(series = lab,
                 slope = cf[2], intercept = cf[1],
                 r2 = summary(m)$r.squared, n = nrow(dd))
        }
        
        fits <- if (identical(input$fitline, "all")) {
          fit_one(w, "all points")
        } else {
          bind_rows(lapply(split(w, w$series),
                           function(dd) fit_one(dd, dd$series[1])))
        }
        
        if (!is.null(fits) && nrow(fits) > 0) {
          if (identical(input$fitline, "all")) {
            p <- p + geom_abline(slope = fits$slope[1], intercept = fits$intercept[1],
                                 colour = "black", linewidth = 0.6)
          } else {
            p <- p + geom_smooth(method = "lm", se = FALSE,
                                 formula = y ~ x, linewidth = 0.6)
          }
          if (isTRUE(input$fit_eq)) {
            fit_txt <- paste(sprintf("%s: y = %.3fx %+.3f  (R2 = %.3f, n = %d)",
                                     fits$series, fits$slope, fits$intercept,
                                     fits$r2, fits$n), collapse = "\n")
          }
        }
      }
      
      p <- p + labs(x = input$xlab %|""|% input$xvar,
                    y = input$ylab %|""|% input$yvar,
                    subtitle = fit_txt)
      
      # ---------------- depth profile ----------------
    } else if (identical(input$plotmode, "depth")) {
      validate(need("depth_bottom" %in% names(d),
                    "The selected datasets have no depth_bottom column."))
      d <- filter(d, !is.na(depth_bottom))
      validate(need(nrow(d) > 0, "No rows with a depth."))
      
      # A depth profile line must never run across two different cores or two
      # different sampling dates, so the line grouping is built here from the
      # finest identity available, independently of "Split series by".
      grp_cols <- intersect(c("dataset", "core_id", "tree", "trm", "position",
                              "logger", "probe_set", "variable"), names(d))
      d <- d %>%
        mutate(.grp = paste(series, format(date, "%Y-%m-%d"),
                            !!!syms(grp_cols), sep = "_"))
      
      # ---- per-column x range ------------------------------------------
      # ggplot has one x range per panel only when scales are free, and it
      # takes it from the data, so the range is forced here: rows outside the
      # wanted range are dropped and an invisible geom_blank pins each panel
      # to the exact limits. Column 1 uses X min/max, the rest use the
      # "(right col.)" pair; a blank box on either side means "auto".
      fb0 <- input$facet_by %||% ""
      xl_l <- c(nn(input$xmin),   nn(input$xmax))
      xl_r <- c(nn(input$xmin_r), nn(input$xmax_r))
      if (nzchar(fb0) && !all(is.na(xl_r))) {
        per_col_x <- TRUE
        lev  <- if (is.factor(d[[fb0]])) levels(droplevels(d[[fb0]]))
        else as.character(sort(unique(d[[fb0]])))
        ncol_f <- if (identical(fb0, "tree")) ceiling(length(lev) / 2)
        else ceiling(sqrt(length(lev)))
        col_of <- ((seq_along(lev) - 1) %% ncol_f) + 1
        lo <- ifelse(col_of == 1, xl_l[1], xl_r[1])
        hi <- ifelse(col_of == 1, xl_l[2], xl_r[2])
        key <- match(as.character(d[[fb0]]), lev)
        d <- d[(is.na(lo[key]) | d$value >= lo[key]) &
                 (is.na(hi[key]) | d$value <= hi[key]), ]
        validate(need(nrow(d) > 0, "No data inside that x range."))
        pin <- tibble(lev = rep(lev, 2), value = c(lo, hi),
                      depth_bottom = d$depth_bottom[1]) %>%
          filter(!is.na(value))
        # the pinning layer must carry the facet column in exactly the same
        # type as the data, or a factor gets coerced to text and the panels
        # fall back into alphabetical order (10, 4, 6, 7)
        pin$lev <- if (is.factor(d[[fb0]]))
          factor(pin$lev, levels = levels(d[[fb0]])) else
            methods::as(pin$lev, class(d[[fb0]])[1])
        names(pin)[1] <- fb0
      }
      
      p <- ggplot(d, aes(value, depth_bottom, colour = ckey, group = .grp)) +
        scale_y_reverse(breaks = sort(unique(d$depth_bottom))) +
        labs(x = input$xlab %|""|% "value",
             y = input$ylab %|""|% "Depth (cm)")
      p <- p + switch(input$geom,
                      line  = geom_path(alpha = input$alpha, linewidth = input$lwd),
                      point = geom_point(alpha = input$alpha, size = input$psize),
                      both  = list(geom_path(alpha = input$alpha, linewidth = input$lwd),
                                   geom_point(alpha = input$alpha, size = input$psize)),
                      col   = geom_point(alpha = input$alpha, size = input$psize))
      if (per_col_x)
        p <- p + geom_blank(data = pin, aes(value, depth_bottom),
                            inherit.aes = FALSE)
      fb <- input$facet_by %||% ""
      if (nzchar(fb)) {
        p <- p + facet_wrap(as.formula(paste("~", fb)), scales = "free_x",
                            nrow = if (identical(fb, "tree")) 2 else NULL)
      } else if (n_distinct(d$variable) > 1) {
        p <- p + facet_wrap(~ variable, nrow = 1, scales = "free_x")
      }
      
      # ---------------- time series, secondary axis ----------------
    } else if (identical(input$mode, "sec") && !is.null(input$sec_var) &&
               input$sec_var %in% d$variable && n_distinct(d$variable) > 1) {
      
      left  <- filter(d, variable != input$sec_var)
      right <- filter(d, variable == input$sec_var)
      validate(need(nrow(left) > 0 && nrow(right) > 0,
                    "Need at least one variable on each axis."))
      lr <- range(left$value, na.rm = TRUE); rr <- range(right$value, na.rm = TRUE)
      sf <- diff(lr) / diff(rr); if (!is.finite(sf) || sf == 0) sf <- 1
      sh <- lr[1] - rr[1] * sf
      
      d2 <- bind_rows(left, mutate(right, value = value * sf + sh))
      p <- add_geoms(ggplot(d2, aes(datetime, value, colour = ckey, group = series))) +
        scale_y_continuous(
          name = input$ylab %|""|% paste(setdiff(input$vars, input$sec_var),
                                         collapse = ", "),
          sec.axis = sec_axis(~ (. - sh) / sf,
                              name = (input$ylab2 %||% "") %|""|% input$sec_var)) +
        labs(x = input$xlab %|""|% NULL)
      
      # ---------------- time series, one axis or facets ----------------
    } else {
      p <- add_geoms(ggplot(d, aes(datetime, value, colour = ckey, group = series))) +
        labs(x = input$xlab %|""|% NULL,
             y = input$ylab %|""|% {
               u <- unique(na.omit(d$unit)); if (length(u) == 1) u else "value"
             })
      if (identical(input$mode, "facet")) {
        labs_map <- setNames(
          map_chr(input$vars, function(v) {
            val <- input[[paste0("flab_", make.names(v))]]
            if (is.null(val) || !nzchar(val)) v else val
          }), input$vars)
        if (identical(input$strippos, "left")) {
          # Putting the strip on the left and outside the axis gives each panel
          # what looks like its own y-axis title - ggplot itself allows only one
          # y title per plot, so this is the standard way to label panels
          # separately (permille on one, % on another).
          p <- p + facet_wrap(~ variable, ncol = 1, scales = "free_y",
                              strip.position = "left",
                              labeller = as_labeller(labs_map)) +
            labs(y = NULL) +
            theme(strip.background = element_blank(),
                  strip.placement  = "outside")
        } else {
          p <- p + facet_wrap(~ variable, ncol = 1, scales = "free_y",
                              labeller = as_labeller(labs_map))
        }
      }
    }
    
    leg_name <- if (identical(input$scheme, "date")) "Sampling date"
    else if (identical(input$scheme, "trm") && identical(input$plotmode, "depth"))
      "Treatment | sampling date" else NULL
    p <- p + scale_colour_manual(values = cols, name = leg_name)
    if (identical(input$geom, "col") && identical(input$plotmode, "ts"))
      p <- p + scale_fill_manual(values = cols, name = NULL)
    
    z <- zoom()
    if (!is.null(z) && identical(input$plotmode, "ts")) {
      p <- p + coord_cartesian(
        xlim = as.POSIXct(z$x, origin = "1970-01-01",
                          tz = attr(d$datetime, "tzone") %||% "UTC"),
        ylim = z$y, expand = FALSE)
    } else if (!is.null(z)) {
      p <- p + coord_cartesian(xlim = z$x, ylim = z$y)
    }
    
    # Manual axis range. A blank box stays NA, which coord_cartesian reads as
    # "auto", so one side can be fixed and the other left free. coord_cartesian
    # only zooms the viewport, so limits apply to EVERY panel of a facet even
    # with scales = "free", and no data is dropped.
    xl <- if (per_col_x) c(NA_real_, NA_real_) else c(nn(input$xmin), nn(input$xmax))
    yl <- c(nn(input$ymin), nn(input$ymax))
    if (!all(is.na(c(xl, yl)))) {
      p <- p + coord_cartesian(xlim = if (all(is.na(xl))) NULL else xl,
                               ylim = if (all(is.na(yl))) NULL else yl)
    }
    
    bs <- input$basesize %||% 12
    
    # A bottom legend is laid out in ONE row by default, so with long keys it
    # runs off the side of the saved file. Wrap it over several columns,
    # estimated from the export width, the font size and the longest label.
    nkey <- n_distinct(d$ckey)
    lc <- input$legcols %||% 0
    if (!is.finite(lc) || lc < 1) {
      lab_cm <- max(nchar(as.character(unique(d$ckey)))) * 0.020 * bs + 1
      lc <- max(1, min(nkey, floor((input$w %||% 25) / lab_cm)))
    }
    # the legend title sits to the LEFT of the keys by default and eats the
    # width the keys need, so it is moved on top
    p <- p + guides(
      colour = guide_legend(ncol = lc, byrow = TRUE, title.position = "top"),
      fill   = guide_legend(ncol = lc, byrow = TRUE, title.position = "top"))
    
    p +
      labs(title = if (nzchar(input$title)) input$title else NULL) +
      theme_bw(base_size = bs) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = bs),
            axis.text   = element_text(size = bs),
            axis.title  = element_text(size = bs),
            strip.text  = element_text(size = bs))
  })
  
  # Static ggplot by default (brush zoom + exact export). With the box ticked a
  # time series goes to plotly instead: hovermode "closest" reports the nearest
  # POINT only, rather than a running readout like dygraphs.
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
    p <- plotly::ggplotly(build_plot(), tooltip = c("x", "y", "colour"))
    p <- plotly::plotly_build(p)
    # the colour aesthetic is called ckey internally; the hover box should not
    # say so
    p$x$data <- lapply(p$x$data, function(tr) {
      if (!is.null(tr$text))
        tr$text <- gsub("ckey:", "series:", tr$text, fixed = TRUE)
      tr
    })
    plotly::layout(p, hovermode = "closest",
                   legend = list(orientation = "h", y = -0.15))
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
