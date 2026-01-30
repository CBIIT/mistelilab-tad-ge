# Figure 2 (D-E), Measurement of boundary distances
Faisal Almansour
January 30, 2026

Analysis setup If you cloned or copied this folder, restore the exact
package versions with: `renv::restore()`

``` r
here::i_am("figure_2_d-e/260121_TAD-boundary_vs_non-TAD-region.qmd")
```

    here() starts at /Users/pegorarog/Documents/manuscript_repos/mistelilab-tad-ge

``` r
library(tidyverse)
```

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.4     ✔ readr     2.1.6
    ✔ forcats   1.0.1     ✔ stringr   1.6.0
    ✔ ggplot2   4.0.1     ✔ tibble    3.3.1
    ✔ lubridate 1.9.4     ✔ tidyr     1.3.2
    ✔ purrr     1.2.1     

    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(here)
library(SpatialTools)
```

    # This research was partially supported under NSF Grant ATM-0534173

``` r
library(fs)
library(reshape2)
```


    Attaching package: 'reshape2'

    The following object is masked from 'package:tidyr':

        smiths

``` r
library(data.table)
```


    Attaching package: 'data.table'

    The following objects are masked from 'package:reshape2':

        dcast, melt

    The following objects are masked from 'package:lubridate':

        hour, isoweek, isoyear, mday, minute, month, quarter, second, wday,
        week, yday, year

    The following objects are masked from 'package:dplyr':

        between, first, last

    The following object is masked from 'package:purrr':

        transpose

``` r
library(ggthemes)
library(ggpointdensity)
library(FSA)
```

    ## FSA v0.10.1. See citation('FSA') if used in publication.
    ## Run fishR() for related website and fishR('IFAR') for related book.

Adjust and standarize figures width and height

Plot text size settings, set the relative size for all the fonts in all
the plots. A value of 1 keeps is as it is. Values \> 1 make fonts
larger. Value \< 1 make fonts smaller. Play around with it.

``` r
new_size <- 1.5
```

Image acquisition and analysis conditions Specify x,y pixel size
(microns) and z-step (microns).

``` r
xy_res <- 0.108
z_step <- 0.1
```

Read the metadata

``` r
metadata_index <- tibble::tibble(
  experiment = c("HFF_EGFR", "HFF_MYC", "HBEC_EGFR", "HBEC_MYC"),
  metadata_path = here::here(c(
    "figure_2_d-e/layouts/hff_egfr_layout.txt",
    "figure_2_d-e/layouts/hff_myc_layout.txt",
    "figure_2_d-e/layouts/hbec_egfr_layout.txt",
    "figure_2_d-e/layouts/hbec_myc_layout.txt"
  ))
)

read_one_md <- function(exp_label, md_path) {
  readr::read_tsv(md_path, show_col_types = FALSE) %>%
    dplyr::mutate(experiment = exp_label)
}

md_tbl <- purrr::map2_dfr(
  metadata_index$experiment,
  metadata_index$metadata_path,
  read_one_md
) %>%
  dplyr::filter(!is.na(treat)) %>%
  dplyr::mutate(
    treat = factor(treat),
    treat_conc = factor(treat_conc, levels = sort(unique(treat_conc))),
    row = as.character(row),
    column = as.character(column),
    experiment_key = tolower(experiment)
  )

md_tbl_join <- md_tbl %>%
  dplyr::distinct(experiment_key, row, column, .keep_all = TRUE) %>%
  dplyr::select(experiment_key, row, column, treat, treat_conc)
```

### Plates layout

![](output/treat_layout-1.pdf)

![](output/treat_layout-2.pdf)

![](output/treat_layout-3.pdf)

![](output/treat_layout-4.pdf)

Read object-level data

``` r
data_index <- tibble::tibble(
  experiment = c("HFF_EGFR", "HFF_MYC", "HBEC_EGFR", "HBEC_MYC"),
  data_path = here::here(
    "figure_2_d-e/",
    "data",
    c("hff_egfr", "hff_myc", "hbec_egfr", "hbec_myc")
  )
)
stopifnot(all(fs::dir_exists(data_index$data_path)))
```

Set `glob` patterns for directory searches for cell level data and spot
data,resepectively. The spot level data is stored in a single .csv file
per well where all the spots in all the channels are stored. The
`channel` variable value indicates the channel in which the spot was
detected according to this mapping: - 2: Green (488) - 3: Red (561) - 4:
Far Red (640)

Read and process the cell-level data. Filter nuclei that are irregularly
shaped (`solidity`) and that are too small (`area`).

``` r
cell_tbl <- data_index %>%
  dplyr::mutate(
    cell_csvs = purrr::map(
      data_path,
      ~ fs::dir_ls(.x, recurse = TRUE, glob = "*well_nuclei_results/*.csv")
    )
  ) %>%
  tidyr::unnest(cell_csvs, names_repair = "unique") %>%
  dplyr::rename(file_name = cell_csvs) %>%
  dplyr::mutate(
    data = purrr::map(file_name, ~ data.table::fread(.x) |> tibble::as_tibble())
  ) %>%
  tidyr::unnest(data) %>%
  dplyr::rename(x = `centroid-0`, y = `centroid-1`) %>%
  dplyr::select(
    experiment,
    column,
    row,
    field_index,
    time_point,
    cell_index,
    x:solidity
  ) %>%
  dplyr::mutate(
    row = as.character(row),
    column = as.character(column),
    experiment_key = tolower(experiment)
  ) %>%
  dplyr::filter(solidity >= 0.95, area >= 10) %>%
  dplyr::arrange(experiment, column, row, field_index, time_point, cell_index)
```

Read the spot-level CSVs (compact, per-experiment folders) Convert the
x, y and z coordinates to microns (They were originally calculated in
pixels) and filter spots that do not belong to any nucleus. These have a
`cell_index` value of `0`.

``` r
spot_tbl <- data_index %>%
  dplyr::mutate(
    spot_csvs = purrr::map(
      data_path,
      ~ fs::dir_ls(.x, recurse = TRUE, glob = "*well_spots_locations/*.csv")
    )
  ) %>%
  tidyr::unnest(spot_csvs, names_repair = "unique") %>%
  dplyr::rename(file_name = spot_csvs) %>%
  dplyr::mutate(
    data = purrr::map(file_name, ~ data.table::fread(.x) |> tibble::as_tibble())
  ) %>%
  tidyr::unnest(data) %>%
  dplyr::rename(
    spot_index_raw = V1,
    x = x_location,
    y = y_location,
    z = z_location
  ) %>%

  dplyr::mutate(
    x_mic = x * xy_res,
    y_mic = y * xy_res,
    z_mic = z * z_step
  ) %>%

  dplyr::filter(cell_index != 0) %>%
  dplyr::mutate(
    row = as.character(row),
    column = as.character(column),
    experiment_key = tolower(experiment)
  ) %>%
  dplyr::select(
    experiment,
    file_name,
    column,
    row,
    field_index,
    time_point,
    cell_index,
    spot_index_raw,
    channel,
    mean_intensity,
    x,
    y,
    z,
    x_mic,
    y_mic,
    z_mic
  ) %>%
  dplyr::arrange(
    experiment,
    column,
    row,
    field_index,
    time_point,
    cell_index,
    channel
  ) %>%
  dplyr::group_by(
    experiment,
    row,
    column,
    field_index,
    time_point,
    cell_index,
    channel
  ) %>%
  dplyr::mutate(spot_index = dplyr::row_number()) %>%
  dplyr::ungroup()

spot_tbl <- spot_tbl %>%
  dplyr::mutate(
    experiment_key = tolower(experiment),
    row = as.character(row),
    column = as.character(column)
  )
```

Use the `spot_tbl` to calculate the number of spots per cell in each
color.

``` r
cell_n_spots_tbl <- spot_tbl %>%
  dplyr::group_by(
    experiment_key,
    row,
    column,
    field_index,
    time_point,
    cell_index,
    channel
  ) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = channel,
    names_prefix = "n_spots_channel_",
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::inner_join(md_tbl_join, by = c("experiment_key", "row", "column"))

cell_counts <- cell_tbl %>%
  dplyr::left_join(
    cell_n_spots_tbl,
    by = c(
      "experiment_key",
      "row",
      "column",
      "field_index",
      "time_point",
      "cell_index"
    )
  ) %>%
  dplyr::left_join(md_tbl_join, by = c("experiment_key", "row", "column"))
```

### Red Histograms, 5’ probe

    [[1]]

![](output/n-red-spots-well-1.pdf)


    [[2]]

![](output/n-red-spots-well-2.pdf)


    [[3]]

![](output/n-red-spots-well-3.pdf)


    [[4]]

![](output/n-red-spots-well-4.pdf)

### Green Histograms, 3’ probe

    [[1]]

![](output/n-green-spots-well-1.pdf)


    [[2]]

![](output/n-green-spots-well-2.pdf)


    [[3]]

![](output/n-green-spots-well-3.pdf)


    [[4]]

![](output/n-green-spots-well-4.pdf)

### FarRed Histograms, nascent RNA probe

    [[1]]

![](output/n-farred-spots-well-1.pdf)


    [[2]]

![](output/n-farred-spots-well-2.pdf)


    [[3]]

![](output/n-farred-spots-well-3.pdf)


    [[4]]

![](output/n-farred-spots-well-4.pdf)

Cell filtering based on spots number Keep only cells that have 2 spots
in the Green channel, 2 spots in the Red channel, and 2 or less spots in
the FarRed channel.

``` r
cell_n_spots_filt_tbl <- cell_n_spots_tbl %>%
  semi_join(
    cell_tbl,
    by = c(
      "row",
      "column",
      "experiment_key",
      "field_index",
      "time_point",
      "cell_index"
    )
  ) %>%
  filter(n_spots_channel_2 == 2, n_spots_channel_3 == 2, n_spots_channel_4 <= 2)

glimpse(cell_n_spots_filt_tbl)
```

    Rows: 19,591
    Columns: 11
    $ experiment_key    <chr> "hbec_egfr", "hbec_egfr", "hbec_egfr", "hbec_egfr", …
    $ row               <chr> "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5…
    $ column            <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11"…
    $ field_index       <dbl> 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3…
    $ time_point        <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ cell_index        <dbl> 16, 4, 6, 8, 9, 11, 12, 14, 20, 23, 30, 33, 36, 4, 6…
    $ n_spots_channel_2 <int> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
    $ n_spots_channel_3 <int> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
    $ n_spots_channel_4 <int> 1, 2, 1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
    $ treat             <fct> Inter-TAD, Inter-TAD, Inter-TAD, Inter-TAD, Inter-TA…
    $ treat_conc        <fct> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…

Filter out spots that do not belong to the correct population of cells
obtained in the previous chunk.

``` r
spot_filt_tbl <- spot_tbl %>%
  semi_join(
    cell_n_spots_filt_tbl,
    by = c(
      "column",
      "row",
      "experiment_key",
      "field_index",
      "time_point",
      "cell_index"
    )
  ) %>%
  ungroup() %>%
  as.data.table()

glimpse(spot_filt_tbl)
```

    Rows: 92,264
    Columns: 18
    $ experiment     <chr> "HBEC_EGFR", "HBEC_EGFR", "HBEC_EGFR", "HBEC_EGFR", "HB…
    $ file_name      <fs::path> "/Users/pegorarog/Documents/manuscript_repos/miste…
    $ column         <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11", "…
    $ row            <chr> "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", …
    $ field_index    <dbl> 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
    $ time_point     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ cell_index     <dbl> 16, 16, 16, 16, 16, 4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 8,…
    $ spot_index_raw <int> 34, 35, 24, 25, 19, 6, 7, 3, 4, 4, 5, 11, 16, 9, 13, 7,…
    $ channel        <dbl> 2, 2, 3, 3, 4, 2, 2, 3, 3, 4, 4, 2, 2, 3, 3, 4, 2, 2, 3…
    $ mean_intensity <dbl> 504.8000, 423.5600, 291.3000, 294.7333, 437.5122, 454.7…
    $ x              <dbl> 1917.7483, 1922.6778, 1921.6621, 1923.9839, 1924.6041, …
    $ y              <dbl> 528.0721, 563.4961, 529.0953, 560.0783, 559.7089, 1502.…
    $ z              <dbl> -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,…
    $ x_mic          <dbl> 207.11682, 207.64920, 207.53951, 207.79026, 207.85724, …
    $ y_mic          <dbl> 57.03179, 60.85758, 57.14229, 60.48845, 60.44857, 162.2…
    $ z_mic          <dbl> -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -…
    $ spot_index     <int> 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1…
    $ experiment_key <chr> "hbec_egfr", "hbec_egfr", "hbec_egfr", "hbec_egfr", "hb…

Calculate 2D spot distances for the Green-Red pairs.

``` r
gr_dist_2D <- spot_filt_tbl[
  channel %in% 2:3, # Only Green and Red spots
  reshape2::melt(
    dist2(
      as.matrix(.SD[channel == 2, list(x_mic, y_mic, z_mic)]),
      as.matrix(.SD[channel == 3, list(x_mic, y_mic, z_mic)])
    ),
    value.name = "gr_dist"
  ),
  by = list(row, column, experiment_key, field_index, time_point, cell_index)
]

setnames(gr_dist_2D, c("Var1", "Var2"), c("g_index", "r_index"))

glimpse(gr_dist_2D)
```

Calculate 2D spot distances for the FarRed-Red pairs.

``` r
fr_dist_2D <- spot_filt_tbl[
  channel %in% 3:4, # Only Red and Far Red spots
  reshape2::melt(
    dist2(
      as.matrix(.SD[channel == 4, list(x_mic, y_mic, z_mic), ]),
      as.matrix(.SD[channel == 3, list(x_mic, y_mic, z_mic)])
    ),
    value.name = "fr_dist"
  ),
  by = list(row, column, experiment_key, field_index, time_point, cell_index)
]

setnames(fr_dist_2D, c("Var1", "Var2"), c("f_index", "r_index"))

glimpse(fr_dist_2D)
```

    Rows: 27,800
    Columns: 9
    $ row            <chr> "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", …
    $ column         <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11", "…
    $ experiment_key <chr> "hbec_egfr", "hbec_egfr", "hbec_egfr", "hbec_egfr", "hb…
    $ field_index    <dbl> 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
    $ time_point     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ cell_index     <dbl> 16, 16, 4, 4, 4, 4, 6, 6, 8, 8, 9, 9, 9, 9, 11, 11, 11,…
    $ f_index        <int> 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ r_index        <int> 1, 2, 1, 1, 2, 2, 1, 2, 1, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1…
    $ fr_dist        <dbl> 3.32150924, 0.07795967, 0.19346924, 2.68443824, 2.71897…

Calculate 2D spot distances for the FarRed-Green pairs.

``` r
gf_dist_2D <- spot_filt_tbl[
  channel %in% c(2, 4), # Only Green and Far Red spots
  reshape2::melt(
    dist2(
      as.matrix(.SD[channel == 2, list(x_mic, y_mic, z_mic), ]),
      as.matrix(.SD[channel == 4, list(x_mic, y_mic, z_mic)])
    ),
    value.name = "gf_dist"
  ),
  by = list(row, column, experiment_key, field_index, time_point, cell_index)
]

setnames(gf_dist_2D, c("Var1", "Var2"), c("g_index", "f_index"))

glimpse(gf_dist_2D)
```

    Rows: 27,800
    Columns: 9
    $ row            <chr> "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", …
    $ column         <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11", "…
    $ experiment_key <chr> "hbec_egfr", "hbec_egfr", "hbec_egfr", "hbec_egfr", "hb…
    $ field_index    <dbl> 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
    $ time_point     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ cell_index     <dbl> 16, 16, 4, 4, 4, 4, 6, 6, 8, 8, 9, 9, 9, 9, 11, 11, 11,…
    $ g_index        <int> 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ f_index        <int> 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1…
    $ gf_dist        <dbl> 3.4960832, 0.4588800, 0.7842695, 2.9300260, 2.3549229, …

Green-Red (3’-5’) proximity calculations. **Important: this calculations
find the closest Green/Red pairs (by calculating the the minimum
distance) on a per Red Spot basis**.

``` r
setkey(
  gr_dist_2D,
  row,
  column,
  experiment_key,
  field_index,
  time_point,
  cell_index,
  r_index
)

gr_dist_min_2D <- gr_dist_2D[, .SD[which.min(gr_dist), ], by = key(gr_dist_2D)]

glimpse(gr_dist_min_2D)
```

    Rows: 39,182
    Columns: 9
    $ row            <chr> "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", …
    $ column         <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11", "…
    $ experiment_key <chr> "hbec_myc", "hbec_myc", "hbec_myc", "hbec_myc", "hbec_m…
    $ field_index    <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ time_point     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ cell_index     <dbl> 2, 2, 4, 4, 5, 5, 6, 6, 9, 9, 10, 10, 12, 12, 14, 14, 1…
    $ r_index        <int> 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ g_index        <int> 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ gr_dist        <dbl> 0.51075083, 0.82461020, 0.53221283, 1.12539840, 0.98480…

Red-FarRed (5’-nRNA) proximity calculations. **Important: this
calculations find the closest Red/FarRed pairs (by calculating the the
minimum distance) on a per Red Spot basis**.

``` r
setkey(
  fr_dist_2D,
  row,
  column,
  experiment_key,
  field_index,
  cell_index,
  time_point,
  r_index
)

fr_dist_min_2D <- fr_dist_2D[, .SD[which.min(fr_dist), ], by = key(fr_dist_2D)]

glimpse(fr_dist_min_2D)
```

    Rows: 18,188
    Columns: 9
    $ row            <chr> "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", …
    $ column         <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11", "…
    $ experiment_key <chr> "hbec_myc", "hbec_myc", "hbec_myc", "hbec_myc", "hbec_m…
    $ field_index    <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2…
    $ cell_index     <dbl> 10, 10, 12, 12, 14, 14, 17, 17, 24, 24, 25, 25, 27, 27,…
    $ time_point     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ r_index        <int> 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ f_index        <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ fr_dist        <dbl> 4.60294986, 0.38693876, 1.96761886, 0.38552158, 0.13233…

Green-FarRed (3’-nRNA) proximity calculations. **Important: this
calculations find the closest Green/FarRed pairs (by calculating the the
minimum distance) on a per Green Spot basis**.

``` r
setkey(
  gf_dist_2D,
  row,
  column,
  experiment_key,
  field_index,
  cell_index,
  time_point,
  g_index
)

gf_dist_min_2D <- gf_dist_2D[, .SD[which.min(gf_dist), ], by = key(gf_dist_2D)]

glimpse(gf_dist_min_2D)
```

    Rows: 18,188
    Columns: 9
    $ row            <chr> "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", …
    $ column         <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11", "…
    $ experiment_key <chr> "hbec_myc", "hbec_myc", "hbec_myc", "hbec_myc", "hbec_m…
    $ field_index    <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2…
    $ cell_index     <dbl> 10, 10, 12, 12, 14, 14, 17, 17, 24, 24, 25, 25, 27, 27,…
    $ time_point     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ g_index        <int> 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ f_index        <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ gf_dist        <dbl> 4.1877128, 0.9431963, 2.6363399, 0.5037341, 0.5284627, …

Join GR (3’-5’) and FR (nRNA-5’) distances tables based on the red spot
index (5’). We use `full_join` because we want to keep pairs of
Red/Green probes that have no matching FarRed spots, or, in other words,
inactive alleles. These are labelled with `NA` for both the `f_index`
and for `fr_dist`.

``` r
gfr_dist_min_2D <- gr_dist_min_2D %>%
  full_join(
    fr_dist_min_2D,
    by = c(
      "row",
      "column",
      "experiment_key",
      "field_index",
      "cell_index",
      "time_point",
      "r_index"
    )
  ) %>%
  inner_join(md_tbl, by = c("row", "column", "experiment_key"))

glimpse(gfr_dist_min_2D)
```

    Rows: 39,182
    Columns: 15
    $ row            <chr> "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", …
    $ column         <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11", "…
    $ experiment_key <chr> "hbec_myc", "hbec_myc", "hbec_myc", "hbec_myc", "hbec_m…
    $ field_index    <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ time_point     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ cell_index     <dbl> 2, 2, 4, 4, 5, 5, 6, 6, 9, 9, 10, 10, 12, 12, 14, 14, 1…
    $ r_index        <int> 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ g_index        <int> 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ gr_dist        <dbl> 0.51075083, 0.82461020, 0.53221283, 1.12539840, 0.98480…
    $ f_index        <int> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1, 1, 1, 1, 1, …
    $ fr_dist        <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 4.6029499, 0.38…
    $ treat          <fct> Inter-TAD, Inter-TAD, Inter-TAD, Inter-TAD, Inter-TAD, …
    $ treat_conc     <fct> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    $ rep            <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", …
    $ experiment     <chr> "HBEC_MYC", "HBEC_MYC", "HBEC_MYC", "HBEC_MYC", "HBEC_M…

Generate a filtered dataset with only FarRed spots.

``` r
fr_spot_filt_tbl <- spot_filt_tbl %>%
  filter(channel == 4) %>%
  select(
    row,
    column,
    experiment_key,
    field_index,
    time_point,
    cell_index,
    f_index = spot_index,
    mean_intensity
  )

glimpse(fr_spot_filt_tbl)
```

    Rows: 13,900
    Columns: 8
    $ row            <chr> "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", …
    $ column         <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11", "…
    $ experiment_key <chr> "hbec_egfr", "hbec_egfr", "hbec_egfr", "hbec_egfr", "hb…
    $ field_index    <dbl> 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2…
    $ time_point     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ cell_index     <dbl> 16, 4, 4, 6, 8, 9, 9, 11, 11, 12, 14, 14, 20, 20, 23, 2…
    $ f_index        <int> 1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ mean_intensity <dbl> 437.5122, 163.0000, 332.4400, 264.8261, 327.7200, 273.2…

``` r
gfr_dist_min_2D <- gfr_dist_min_2D %>%
  #left_join(md_tbl, on = c("row", "column")) %>%
  left_join(
    fr_spot_filt_tbl,
    by = c(
      "row",
      "column",
      "experiment_key",
      "field_index",
      "time_point",
      "cell_index",
      "f_index"
    )
  )

glimpse(gfr_dist_min_2D)
```

    Rows: 39,182
    Columns: 16
    $ row            <chr> "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", …
    $ column         <chr> "11", "11", "11", "11", "11", "11", "11", "11", "11", "…
    $ experiment_key <chr> "hbec_myc", "hbec_myc", "hbec_myc", "hbec_myc", "hbec_m…
    $ field_index    <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ time_point     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
    $ cell_index     <dbl> 2, 2, 4, 4, 5, 5, 6, 6, 9, 9, 10, 10, 12, 12, 14, 14, 1…
    $ r_index        <int> 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ g_index        <int> 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1…
    $ gr_dist        <dbl> 0.51075083, 0.82461020, 0.53221283, 1.12539840, 0.98480…
    $ f_index        <int> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1, 1, 1, 1, 1, …
    $ fr_dist        <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 4.6029499, 0.38…
    $ treat          <fct> Inter-TAD, Inter-TAD, Inter-TAD, Inter-TAD, Inter-TAD, …
    $ treat_conc     <fct> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0…
    $ rep            <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", …
    $ experiment     <chr> "HBEC_MYC", "HBEC_MYC", "HBEC_MYC", "HBEC_MYC", "HBEC_M…
    $ mean_intensity <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 198.0769, 198.0…

### Red-FarRed distances Single Allele Density Plots

    [[1]]

![](output/gfr-fr-dist-density-pooled-1.pdf)


    [[2]]

![](output/gfr-fr-dist-density-pooled-2.pdf)


    [[3]]

![](output/gfr-fr-dist-density-pooled-3.pdf)


    [[4]]

![](output/gfr-fr-dist-density-pooled-4.pdf)

### Figure 2 (D-E)

    [[1]]

![](output/figure%202%20(D-E)-1.pdf)


    [[2]]

![](output/figure%202%20(D-E)-2.pdf)


    [[3]]

![](output/figure%202%20(D-E)-3.pdf)


    [[4]]

![](output/figure%202%20(D-E)-4.pdf)

### Numerical summaries (Mean, SD, Median, IQR)

Summarize **without** taking FarRed Status in account on a per well
basis.

``` r
well_2D <- gfr_dist_min_2D %>%
  group_by(row, column, experiment_key, treat, treat_conc) %>%
  summarize(across(
    c(gr_dist, fr_dist),
    list(
      mean = ~ mean(.x, na.rm = T),
      median = ~ median(.x, na.rm = T),
      sd = ~ sd(.x, na.rm = T),
      iqr = ~ IQR(.x, na.rm = TRUE)
    )
  ))

write_csv(well_2D, here("figure_2_d-e/output/well_2D.csv"))

knitr::kable(well_2D, digits = 2)
```

| row | column | experiment_key | treat | treat_conc | gr_dist_mean | gr_dist_median | gr_dist_sd | gr_dist_iqr | fr_dist_mean | fr_dist_median | fr_dist_sd | fr_dist_iqr |
|:---|:---|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 11 | hbec_myc | Inter-TAD | 0 | 0.74 | 0.71 | 0.42 | 0.50 | 2.82 | 0.54 | 3.13 | 5.17 |
| 3 | 13 | hbec_myc | Intra-TAD | 0 | 0.49 | 0.46 | 0.33 | 0.33 | 2.82 | 0.53 | 3.20 | 5.02 |
| 5 | 11 | hbec_egfr | Inter-TAD | 0 | 0.49 | 0.42 | 0.58 | 0.29 | 1.18 | 0.26 | 2.18 | 0.27 |
| 5 | 13 | hbec_egfr | Intra-TAD | 0 | 0.45 | 0.39 | 0.51 | 0.28 | 1.10 | 0.26 | 2.08 | 0.26 |
| 7 | 11 | hff_myc | Inter-TAD | 0 | 0.94 | 0.81 | 0.87 | 0.63 | 5.58 | 4.66 | 5.37 | 8.73 |
| 7 | 13 | hff_myc | Intra-TAD | 0 | 0.47 | 0.41 | 0.52 | 0.30 | 5.37 | 4.41 | 5.23 | 8.66 |
| 9 | 11 | hff_egfr | Inter-TAD | 0 | 0.80 | 0.56 | 1.20 | 0.48 | 2.11 | 0.27 | 3.87 | 0.58 |
| 9 | 13 | hff_egfr | Intra-TAD | 0 | 0.58 | 0.41 | 1.10 | 0.30 | 2.16 | 0.27 | 4.06 | 0.60 |

Summarize **without** taking FarRed Status in account by taking the
Mean, SD, Median, and IQR of wells per condition.

``` r
pool_2D <- well_2D %>%
  group_by(treat, experiment_key, treat_conc) %>%
  summarize(across(
    c(gr_dist_mean, gr_dist_median, fr_dist_mean, fr_dist_median),
    list(
      mean = ~ mean(.x, na.rm = T),
      sd = ~ sd(.x, na.rm = T),
      iqr = ~ IQR(.x, na.rm = TRUE)
    )
  ))

write_csv(pool_2D, here("figure_2_d-e/output/fr_status_pool_2D.csv"))

knitr::kable(pool_2D, digits = 2)
```

| treat | experiment_key | treat_conc | gr_dist_mean_mean | gr_dist_mean_sd | gr_dist_mean_iqr | gr_dist_median_mean | gr_dist_median_sd | gr_dist_median_iqr | fr_dist_mean_mean | fr_dist_mean_sd | fr_dist_mean_iqr | fr_dist_median_mean | fr_dist_median_sd | fr_dist_median_iqr |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Inter-TAD | hbec_egfr | 0 | 0.49 | NA | 0 | 0.42 | NA | 0 | 1.18 | NA | 0 | 0.26 | NA | 0 |
| Inter-TAD | hbec_myc | 0 | 0.74 | NA | 0 | 0.71 | NA | 0 | 2.82 | NA | 0 | 0.54 | NA | 0 |
| Inter-TAD | hff_egfr | 0 | 0.80 | NA | 0 | 0.56 | NA | 0 | 2.11 | NA | 0 | 0.27 | NA | 0 |
| Inter-TAD | hff_myc | 0 | 0.94 | NA | 0 | 0.81 | NA | 0 | 5.58 | NA | 0 | 4.66 | NA | 0 |
| Intra-TAD | hbec_egfr | 0 | 0.45 | NA | 0 | 0.39 | NA | 0 | 1.10 | NA | 0 | 0.26 | NA | 0 |
| Intra-TAD | hbec_myc | 0 | 0.49 | NA | 0 | 0.46 | NA | 0 | 2.82 | NA | 0 | 0.53 | NA | 0 |
| Intra-TAD | hff_egfr | 0 | 0.58 | NA | 0 | 0.41 | NA | 0 | 2.16 | NA | 0 | 0.27 | NA | 0 |
| Intra-TAD | hff_myc | 0 | 0.47 | NA | 0 | 0.41 | NA | 0 | 5.37 | NA | 0 | 4.41 | NA | 0 |

### Frequency of loops

Classify looping status Establish that: - green-red (3’-5’) spots who
are equal to or smaller than 2.5 microns in 2D are looped - green-red
(3’-5’) spots who are larger than 2.5 microns in 2D are unlooped

``` r
gfr_dist_min_2D <- gfr_dist_min_2D %>%
  mutate(
    cut_off = case_when(
      gr_dist <= 0.25 ~ "looped",
      gr_dist > 0.25 ~ "unlooped"
    )
  )
```

``` r
loop <- gfr_dist_min_2D %>%
  group_by(treat, experiment_key, cut_off) %>%
  summarize(n = n()) %>%
  mutate(freq = n / sum(n))

write_csv(loop, here("figure_2_d-e/output/loop_freq.csv"))

knitr::kable(loop, digits = 2)
```

| treat     | experiment_key | cut_off  |    n | freq |
|:----------|:---------------|:---------|-----:|-----:|
| Inter-TAD | hbec_egfr      | looped   | 1294 | 0.21 |
| Inter-TAD | hbec_egfr      | unlooped | 4910 | 0.79 |
| Inter-TAD | hbec_myc       | looped   |  812 | 0.08 |
| Inter-TAD | hbec_myc       | unlooped | 9192 | 0.92 |
| Inter-TAD | hff_egfr       | looped   |  155 | 0.14 |
| Inter-TAD | hff_egfr       | unlooped |  967 | 0.86 |
| Inter-TAD | hff_myc        | looped   |  128 | 0.06 |
| Inter-TAD | hff_myc        | unlooped | 1944 | 0.94 |
| Intra-TAD | hbec_egfr      | looped   | 1813 | 0.24 |
| Intra-TAD | hbec_egfr      | unlooped | 5873 | 0.76 |
| Intra-TAD | hbec_myc       | looped   | 1513 | 0.17 |
| Intra-TAD | hbec_myc       | unlooped | 7167 | 0.83 |
| Intra-TAD | hff_egfr       | looped   |  348 | 0.23 |
| Intra-TAD | hff_egfr       | unlooped | 1170 | 0.77 |
| Intra-TAD | hff_myc        | looped   |  397 | 0.21 |
| Intra-TAD | hff_myc        | unlooped | 1499 | 0.79 |

### Statistical tests

Run Mann Whitney Wilcoxon Parametric Test; Compare the DNA/RNA FISH and
MS2 groups

    # A tibble: 4 × 3
      experiment_key     n      pval
      <chr>          <int>     <dbl>
    1 hbec_egfr      13890 1.49e- 12
    2 hbec_myc       18684 0        
    3 hff_egfr        2640 1.08e- 39
    4 hff_myc         3968 6.32e-221

         Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     0.006863  0.369149  0.582058  0.685397  0.864389 19.807720 

         Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
     0.001442  0.278783  0.423525  0.479428  0.583650 17.133025 

    [1] 0.1585328

    [1] 0

Document the information about the analysis session

``` r
sessionInfo()
```

    R version 4.5.2 (2025-10-31)
    Platform: aarch64-apple-darwin20
    Running under: macOS Tahoe 26.2

    Matrix products: default
    BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
    LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: America/New_York
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices datasets  utils     methods   base     

    other attached packages:
     [1] FSA_0.10.1           ggpointdensity_0.2.1 ggthemes_5.2.0      
     [4] data.table_1.18.0    reshape2_1.4.5       fs_1.6.6            
     [7] SpatialTools_1.0.5   here_1.0.2           lubridate_1.9.4     
    [10] forcats_1.0.1        stringr_1.6.0        dplyr_1.1.4         
    [13] purrr_1.2.1          readr_2.1.6          tidyr_1.3.2         
    [16] tibble_3.3.1         ggplot2_4.0.1        tidyverse_2.0.0     

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6       xfun_0.56          lattice_0.22-7     tzdb_0.5.0        
     [5] vctrs_0.7.1        tools_4.5.2        generics_0.1.4     parallel_4.5.2    
     [9] pkgconfig_2.0.3    Matrix_1.7-4       RColorBrewer_1.1-3 S7_0.2.1          
    [13] spBayes_0.4-8      lifecycle_1.0.5    compiler_4.5.2     farver_2.1.2      
    [17] htmltools_0.5.9    yaml_2.3.12        Formula_1.2-5      pillar_1.11.1     
    [21] crayon_1.5.3       MASS_7.3-65        abind_1.4-8        tidyselect_1.2.1  
    [25] digest_0.6.39      stringi_1.8.7      labeling_0.4.3     magic_1.6-1       
    [29] rprojroot_2.1.1    fastmap_1.2.0      grid_4.5.2         cli_3.6.5         
    [33] magrittr_2.0.4     utf8_1.2.6         withr_3.0.2        scales_1.4.0      
    [37] sp_2.2-0           bit64_4.6.0-1      timechange_0.3.0   rmarkdown_2.30    
    [41] bit_4.6.0          hms_1.1.4          coda_0.19-4.1      evaluate_1.0.5    
    [45] knitr_1.51         rlang_1.1.7        Rcpp_1.1.1         glue_1.8.0        
    [49] renv_1.1.6         vroom_1.6.7        jsonlite_2.0.0     R6_2.6.1          
    [53] plyr_1.8.9        
