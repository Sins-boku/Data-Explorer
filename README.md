# Data Explorer

Data Explorer is an R Shiny application for inspecting environmental
measurements as time series, depth profiles, or dual-isotope plots. It accepts
CSV and RDS data files, normalizes their columns for plotting, and uses
gap-aware multi-resolution summaries to keep wide time-series views manageable.

## Requirements

- R 4.1 or later
- The R packages listed under `Imports` in [`DESCRIPTION`](DESCRIPTION)
- The `Env-Plot` Git submodule, which provides the shared renderer

## Setup and run

Clone the repository with its submodule:

```sh
git clone --recurse-submodules https://github.com/Sins-boku/Data-Explorer.git
cd Data-Explorer
```

If you already cloned without submodules, initialize them with:

```sh
git submodule update --init --recursive
```

Install the app's declared dependencies from an R session:

```r
install.packages(c(
  "data.table", "dplyr", "ggplot2", "plotly", "promises", "purrr",
  "readr", "scales", "shiny", "stringr", "tibble", "tidyr"
))
```

Start the app from the repository root:

```r
shiny::runApp()
```

Alternatively, open `app.R` in RStudio and select **Run App**. The app expects
the `Env-Plot` submodule directory to be present relative to the repository
root.

## Add data

Place `.csv` or `.rds` files directly in either `data_input/` or `data/`.
The app discovers files in these directories while running. RDS files must
contain a data frame or tibble. Raw source files are read but not modified.

For long/tabular data, the app needs:

- A timestamp column
- A numeric measurement/value column
- Optionally, a variable/parameter column and a unit column
- Any grouping columns needed to distinguish sensors, stations, depths, or
  other independent series

Common timestamp and measurement names are detected automatically. When
detection is ambiguous or a file uses other names, select the appropriate
columns in the per-file mapping controls. Choose the timezone used to interpret
timestamps; UTC is the default. Group columns determine which records belong to
separate series and are not combined during summarization.

If no files are found, create one of the input directories and add a supported
file. The app refreshes its file list periodically.

## Explore and export

The sidebar provides controls for:

- Selecting datasets, variables, series/group levels, dates, and depths
- Switching between time series, depth profiles, and dual-isotope plots
- Adjusting gap detection, axes, facets, geometry, colors, and plot appearance
- Zooming and panning interactive time-series plots; summary detail is selected
  according to the visible time span and plot width
- Downloading the displayed plot as PNG or PDF

The time-series renderer uses bin means for summary values and retains
minimum/maximum envelopes to make extrema visible. Missing or non-finite
measurements and detected long gaps separate line segments rather than being
treated as zero or connected across a gap. Summary statistics are
observation-weighted.

## Preprocessing and cache

Normalized data and summaries are cached under `data_processed/.cache/`.
The cache is generated automatically and can be removed to force preprocessing
again; source data is left unchanged. Cache entries are invalidated when the
source file identity or preprocessing configuration changes.

The process-level memory cache defaults to 512 MiB, and the disk summary cache
defaults to 5 GiB. These limits can be changed before starting the app:

```r
options(
  dataExplorer.preprocessing_cache_bytes = 512 * 1024^2,
  dataExplorer.preprocessing_disk_cache_bytes = 5 * 1024^3
)
shiny::runApp()
```

The memory cache is per R process; setting its limit to `0` prevents entries
from being retained there. Disk pruning keeps the summary for the file
currently being processed, even if that one file exceeds the configured disk
limit. The disk cache is shared by app sessions using the same working
directory. The current implementation does not schedule preprocessing on a
parallel worker.

## Tests

Run the preprocessing and renderer tests from the repository root:

```r
testthat::test_dir("tests/testthat", reporter = "summary")
```

## Tasks not completed / limitations

- The app has not been exercised end to end in a browser. In particular,
  interactive zoom/pan refresh behavior and all plot modes still need
  integration verification.
- Renderer tests were added, but the final combined test run was interrupted;
  only the preprocessing tests had a confirmed passing run. Do not treat the
  renderer tests as verified.
- There are no representative project datasets in the repository. The input
  mapping and plotting behavior therefore still need validation against the
  datasets users intend to analyze.
- Multi-session cache contention, simultaneous preprocessing, and stale
  overlapping viewport requests have not been integration-tested. Processing
  is synchronous rather than delegated to a background worker.
- Peak memory and end-to-end latency have not been measured on the implemented
  app path with representative production data. A synthetic million-row
  benchmark informed the initial storage and resolution choices, but does not
  replace those measurements.
