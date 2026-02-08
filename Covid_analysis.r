## COVID-19 cumulative rates and neighbourhood determinants
options(stringsAsFactors = FALSE)

# Libraries actually used
library(dplyr)
library(readr)
library(readxl)
library(writexl)
library(sf)
library(spdep)
library(spatialreg)
library(tmap)
library(car)
library(GWmodel)
library(spgwr)
library(randomForest)
library(ggplot2)
library(leaflet)
library(leaflet.extras)
library(viridisLite)
library(scales)

tmap_mode("plot")

safe_log10 <- function(x) ifelse(x > 0, log10(x), NA_real_)
model_metrics <- function(model, observed) {
  preds <- fitted(model)
  res   <- resid(model)
  tss   <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  rss   <- sum((observed - preds)^2, na.rm = TRUE)
  r2    <- 1 - (rss / tss)
  rse   <- sqrt(sum(res^2, na.rm = TRUE) / (length(res) - length(coef(model))))
  list(r_squared = r2, rse = rse)
}

# Config ---------------------------------------------------------------------
config_file <- if (file.exists("config.R")) "config.R" else "config.example.R"
source(config_file)

if (!exists("data_dir") || is.na(data_dir) || identical(data_dir, "")) {
  stop("Set data_dir in config.R or via environment variable COVID_DATA_DIR")
}

data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
message("Data directory: ", data_dir)

find_shapefile <- function(dir) {
  shp <- list.files(dir, pattern = "\\\\.shp$", full.names = TRUE, recursive = TRUE)
  message("Shapefile search found ", length(shp), " candidate(s).")
  if (length(shp) == 0) {
    message("Sample of files seen: ", paste(head(list.files(dir, recursive = TRUE), 10), collapse = "; "))
  }
  if (length(shp) == 0) return(NA_character_)
  shp[1]
}

# Paths
covid_path       <- file.path(data_dir, "ONSMapCOVID_EN.csv")
determinants_in  <- file.path(data_dir, "COVID DETERMINANTS IN ARC.xlsx")  # source
determinants_out <- file.path(data_dir, "covid_determinants_arc.xlsx")     # clean export
default_shp <- file.path(data_dir, "Ottawa_Neighbourhood_Study_(ONS)_-_Neighbourhood_Boundaries_Gen_2.shp")
shapefile_path   <- if (file.exists(default_shp)) default_shp else find_shapefile(data_dir)
message("Shapefile candidate: ", shapefile_path)

stopifnot(file.exists(covid_path), file.exists(determinants_in))
setwd(data_dir)

# Retry shapefile search from working directory if not found yet
if (is.na(shapefile_path) || !file.exists(shapefile_path)) {
  shapefile_path <- find_shapefile(getwd())
  message("Shapefile retry from wd: ", shapefile_path)
}
if (is.na(shapefile_path) || !file.exists(shapefile_path)) {
  fallback_shp <- file.path(data_dir, "Shape file", "shapefile", "ONS_Boundaries_Gen2.shp")
  if (file.exists(fallback_shp)) {
    shapefile_path <- fallback_shp
    message("Shapefile fallback used: ", shapefile_path)
  }
}

outputs_dir <- file.path(data_dir, "outputs")
dir.create(outputs_dir, showWarnings = FALSE, recursive = TRUE)

# Load and clean data
covid_raw <- read_csv(
  covid_path,
  col_types = cols(`Cumulative Rate Excluding Cases Linked to Outbreaks in LTCH & RH` = col_double())
) |>
  rename(
    ONS_ID     = `ONS ID`,
    covid_name = `ONS Neighbourhood Name`,
    covid_rate = `Cumulative Rate Excluding Cases Linked to Outbreaks in LTCH & RH`
  ) |>
  select(ONS_ID, covid_name, covid_rate)

determinants_raw <- read_excel(determinants_in)
if ("ONS_Neighbourhood_Name" %in% names(determinants_raw)) {
  determinants_raw <- determinants_raw |>
    rename(det_name = ONS_Neighbourhood_Name)
} else if (!"det_name" %in% names(determinants_raw)) {
  determinants_raw$det_name <- NA_character_
}

standardize_cols <- c(
  population_density = "Population_density",
  median_income_after_tax = "Median_income_AFT_TX",
  unemployment_rate = "Unemployment_rate",
  perc_public_transport_work = "Perc_pulibc_transport_work",
  perc_newcomers = "Perc_Newcomers",
  avg_household = "Avg_household",
  job_density = "Jobs_per_km2",
  perc_no_hs_diploma = "Perc_no_high_school_diploma",
  perc_bachelors_and_above = "Perc_Baechlors_and_above",
  perc_seniors = "Perc_Seniors"
)

for (new_nm in names(standardize_cols)) {
  old_nm <- standardize_cols[[new_nm]]
  if (old_nm %in% names(determinants_raw) && !(new_nm %in% names(determinants_raw))) {
    determinants_raw <- rename(determinants_raw, !!new_nm := all_of(old_nm))
  }
}

analysis_data <- determinants_raw |>
  left_join(covid_raw, by = "ONS_ID", suffix = c("", ".covid")) |>
  mutate(
    covid_rate = coalesce(covid_rate, covid_rate.covid),
    ONS_Name   = coalesce(det_name, covid_name)
  ) |>
  select(-det_name, -covid_name, -covid_rate.covid) |>
  filter(!is.na(covid_rate)) |>
  mutate(
    population_density_sqrt       = sqrt(population_density),
    unemployment_rate_log         = safe_log10(unemployment_rate),
    perc_newcomers_sqrt           = sqrt(perc_newcomers),
    avg_household_sq              = avg_household^2,
    job_density_log               = safe_log10(job_density),
    perc_no_hs_diploma_sqrt       = sqrt(perc_no_hs_diploma),
    perc_bachelors_and_above_sqrt = sqrt(perc_bachelors_and_above),
    perc_seniors_log              = safe_log10(perc_seniors),
    covid_rate_log                = safe_log10(covid_rate)
  )

# Export cleaned table for ArcGIS
write_xlsx(analysis_data, determinants_out)
message("Wrote cleaned determinants to: ", determinants_out)

# Quick descriptive plots (non-spatial)
covid_rate_plot <- ggplot(analysis_data, aes(x = covid_rate)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 25) +
  labs(title = "COVID cumulative rate distribution", x = "Cumulative rate", y = "Count") +
  theme_minimal()
ggsave(
  filename = file.path(outputs_dir, "covid_rate_hist.png"),
  plot = covid_rate_plot,
  width = 7,
  height = 5,
  dpi = 150
)

cor_mat <- analysis_data |>
  select(where(is.numeric)) |>
  cor(use = "pairwise.complete.obs")
cor_df <- as.data.frame(as.table(cor_mat))
names(cor_df) <- c("Var1", "Var2", "Correlation")
cor_plot <- ggplot(cor_df, aes(Var1, Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white", midpoint = 0, limits = c(-1, 1)) +
  coord_equal() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation heatmap (determinants)", x = NULL, y = NULL)
ggsave(
  filename = file.path(outputs_dir, "correlation_heatmap.png"),
  plot = cor_plot,
  width = 8,
  height = 7,
  dpi = 150
)

# Modeling data
covid_formula <- covid_rate_log ~ population_density_sqrt +
  median_income_after_tax +
  unemployment_rate_log +
  perc_public_transport_work +
  perc_newcomers_sqrt +
  avg_household_sq +
  job_density_log +
  perc_no_hs_diploma_sqrt +
  perc_bachelors_and_above_sqrt +
  perc_seniors_log

ols_model  <- lm(covid_formula, data = analysis_data)
step_model <- step(ols_model, trace = FALSE)
print(summary(step_model))
print(vif(step_model))

# OLS diagnostic plots (saved to outputs/)
png(file.path(outputs_dir, "ols_diagnostics.png"), width = 1200, height = 1200, res = 150)
par(mfrow = c(2, 2))
plot(step_model)
dev.off()

residuals_ols <- residuals(step_model)
print(shapiro.test(residuals_ols))
png(file.path(outputs_dir, "ols_residual_hist.png"), width = 800, height = 600, res = 150)
hist(residuals_ols, main = "Stepwise OLS residuals", xlab = "Residuals")
dev.off()

# Spatial analysis (runs only if shapefile is present)
if (!is.na(shapefile_path) && file.exists(shapefile_path)) {
  ons_shapes <- st_read(shapefile_path, quiet = TRUE) |> st_transform(32189)
  joined_spatial <- ons_shapes |>
    left_join(analysis_data, by = "ONS_ID") |>
    filter(!is.na(covid_rate_log))
  joined_spatial$residuals_ols <- residuals_ols[match(joined_spatial$ONS_ID, analysis_data$ONS_ID)]

  ons_nb <- poly2nb(joined_spatial, queen = TRUE)
  ons_wt <- nb2listw(ons_nb, style = "W", zero.policy = FALSE)

  moran_test <- moran.mc(joined_spatial$covid_rate_log, ons_wt, nsim = 9999, alternative = "greater")
  print(moran_test)

  ons_plot <- tm_shape(joined_spatial) + tm_polygons(col = "covid_rate_log")
  tmap_save(ons_plot, filename = file.path(outputs_dir, "covid_rate_log_map.png"))

  resid_plot <- tm_shape(joined_spatial) + tm_polygons(col = "residuals_ols")
  tmap_save(resid_plot, filename = file.path(outputs_dir, "ols_residuals_map.png"))

  sdem <- errorsarlm(step_model, listw = ons_wt, zero.policy = TRUE, etype = "emixed")
  sem  <- errorsarlm(step_model, listw = ons_wt, zero.policy = TRUE)
  sdm  <- lagsarlm(step_model, data = joined_spatial, listw = ons_wt, type = "mixed")
  sly  <- lagsarlm(step_model, data = joined_spatial, listw = ons_wt)
  print(summary(sdem)); print(summary(sem)); print(summary(sdm)); print(summary(sly, Nagelkerke = TRUE))
  print(model_metrics(sdem, joined_spatial$covid_rate_log))
  print(model_metrics(sem, joined_spatial$covid_rate_log))
  print(model_metrics(sdm, joined_spatial$covid_rate_log))
  print(model_metrics(sly, joined_spatial$covid_rate_log))
  lm.LMtests(step_model, ons_wt, test = "all")

  local_moran <- localmoran(joined_spatial$covid_rate_log, ons_wt)
  joined_spatial$quad_sig <- factor(
    ifelse(local_moran[, 5] < 0.05, "Significant", "Not significant"),
    levels = c("Significant", "Not significant")
  )
  quad_map <- tm_shape(joined_spatial) +
    tm_polygons(col = "quad_sig", palette = c("red", "grey90"), title = "Local Moran's I")
  tmap_save(quad_map, filename = file.path(outputs_dir, "local_moran_clusters.png"))

  # Moran scatterplot (log rate)
  moran_df <- data.frame(
    rate = joined_spatial$covid_rate_log,
    lag_rate = lag.listw(ons_wt, joined_spatial$covid_rate_log)
  )
  moran_plot <- ggplot(moran_df, aes(x = rate, y = lag_rate)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
    theme_minimal() +
    labs(title = "Moran scatterplot (COVID log rate)", x = "covid_rate_log", y = "Spatial lag")
  ggsave(file.path(outputs_dir, "moran_scatterplot.png"), moran_plot, width = 7, height = 5, dpi = 150)

  gwr_formula <- formula(step_model)
  joined_spatial_sp <- as_Spatial(joined_spatial)
  gwr_bw <- bw.gwr(gwr_formula, data = joined_spatial_sp, approach = "CV", kernel = "exponential")
  gwr_model <- gwr(gwr_formula, data = joined_spatial_sp, bandwidth = gwr_bw, gweight = gwr.Gauss, hatmatrix = TRUE)
  gwr_resid <- gwr_model$SDF@data$gwr.e
  print(shapiro.test(gwr_resid))
  print(gwr.morantest(gwr_model, ons_wt))

  mgwr_model <- gwr.multiscale(
    gwr_formula, data = joined_spatial_sp,
    approach = "CV", kernel = "bisquare", adaptive = TRUE, hatmatrix = TRUE, verbose = FALSE
  )
  mgwr_resid <- mgwr_model$SDF@data$residual
  print(shapiro.test(mgwr_resid))
  print(moran.mc(mgwr_resid, ons_wt, nsim = 999))

  # Interactive Leaflet map ---------------------------------------------------
  joined_spatial_wgs <- st_transform(joined_spatial, 4326)

  local_moran_df <- localmoran(joined_spatial$covid_rate_log, ons_wt)
  z_rate <- scale(joined_spatial$covid_rate_log)[, 1]
  z_lag  <- scale(lag.listw(ons_wt, joined_spatial$covid_rate_log))[, 1]
  cluster_cat <- case_when(
    z_rate >= 0 & z_lag >= 0 & local_moran_df[, 5] <= 0.05 ~ "High-High (hot spot)",
    z_rate <= 0 & z_lag <= 0 & local_moran_df[, 5] <= 0.05 ~ "Low-Low (cold spot)",
    z_rate >= 0 & z_lag <= 0 & local_moran_df[, 5] <= 0.05 ~ "High-Low (outlier)",
    z_rate <= 0 & z_lag >= 0 & local_moran_df[, 5] <= 0.05 ~ "Low-High (outlier)",
    TRUE ~ "Not significant"
  )
  joined_spatial_wgs$local_cluster <- factor(
    cluster_cat,
    levels = c("High-High (hot spot)", "Low-Low (cold spot)", "High-Low (outlier)", "Low-High (outlier)", "Not significant")
  )
  joined_spatial_wgs$local_moran_p <- local_moran_df[, 5]

  joined_spatial_wgs$pred_log <- as.numeric(predict(sdem, newdata = joined_spatial, listw = ons_wt))
  joined_spatial_wgs$pred_rate <- 10^joined_spatial_wgs$pred_log
  joined_spatial_wgs$resid_rate <- joined_spatial_wgs$covid_rate - joined_spatial_wgs$pred_rate
  joined_spatial_wgs$resid_log  <- joined_spatial_wgs$covid_rate_log - joined_spatial_wgs$pred_log

  risk_vars <- c("perc_no_hs_diploma", "perc_public_transport_work", "perc_seniors", "unemployment_rate")
  risk_scaled <- scale(st_drop_geometry(joined_spatial_wgs)[, risk_vars])
  income_scaled <- scale(joined_spatial_wgs$median_income_after_tax) * -1
  vuln_raw <- rowMeans(cbind(risk_scaled, income_scaled), na.rm = TRUE)
  joined_spatial_wgs$vulnerability_index <- rescale(
    vuln_raw,
    to = c(0, 100),
    from = range(vuln_raw, na.rm = TRUE)
  )

  joined_spatial_wgs$label_name <- if ("ONS_Name" %in% names(joined_spatial_wgs)) joined_spatial_wgs$ONS_Name else joined_spatial_wgs$ONS_ID

  make_bins <- function(x, n = 5) {
    qs <- quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE, type = 7)
    qs <- unique(qs)
    if (length(qs) < 3) qs <- pretty(x, n = n)
    qs
  }

  pal_rate <- colorNumeric(viridis(7, option = "C"), domain = joined_spatial_wgs$covid_rate, na.color = "transparent")
  pal_cluster <- colorFactor(
    palette = c(
      "High-High (hot spot)" = "#d73027",
      "Low-Low (cold spot)" = "#4575b4",
      "High-Low (outlier)" = "#fdae61",
      "Low-High (outlier)" = "#74add1",
      "Not significant" = "#d9d9d9"
    ),
    domain = joined_spatial_wgs$local_cluster,
    na.color = "transparent"
  )
  bins_nohs     <- make_bins(joined_spatial_wgs$perc_no_hs_diploma, 5)
  bins_transit  <- make_bins(joined_spatial_wgs$perc_public_transport_work, 5)
  bins_seniors  <- make_bins(joined_spatial_wgs$perc_seniors, 5)
  bins_newcomer <- make_bins(joined_spatial_wgs$perc_newcomers, 5)

  pal_nohs     <- colorBin("YlOrRd", domain = joined_spatial_wgs$perc_no_hs_diploma, bins = bins_nohs, na.color = "transparent")
  pal_transit  <- colorBin("YlOrRd", domain = joined_spatial_wgs$perc_public_transport_work, bins = bins_transit, na.color = "transparent")
  pal_seniors  <- colorBin("YlOrRd", domain = joined_spatial_wgs$perc_seniors, bins = bins_seniors, na.color = "transparent")
  pal_newcomer <- colorBin("YlOrRd", domain = joined_spatial_wgs$perc_newcomers, bins = bins_newcomer, na.color = "transparent")
  pal_resid    <- colorNumeric("RdBu", domain = joined_spatial_wgs$resid_rate, na.color = "transparent", reverse = TRUE)
  pal_vuln     <- colorNumeric(inferno(9), domain = joined_spatial_wgs$vulnerability_index, na.color = "transparent")

  popup_template <- sprintf(
    "<b>%s</b><br/>COVID rate: %s per 100k<br/>Predicted: %s<br/>Residual: %+0.1f<br/>%% no HS diploma: %0.1f%%<br/>%% transit to work: %0.1f%%<br/>%% age 65+: %0.1f%%<br/>Vulnerability index: %0.0f",
    joined_spatial_wgs$label_name,
    comma(joined_spatial_wgs$covid_rate, accuracy = 0.1),
    comma(joined_spatial_wgs$pred_rate, accuracy = 0.1),
    joined_spatial_wgs$resid_rate,
    joined_spatial_wgs$perc_no_hs_diploma,
    joined_spatial_wgs$perc_public_transport_work,
    joined_spatial_wgs$perc_seniors,
    joined_spatial_wgs$vulnerability_index
  )

  base_layers <- c("CartoDB.Positron", "Esri.WorldTopoMap", "Esri.WorldImagery")

  map <- leaflet(joined_spatial_wgs, options = leafletOptions(preferCanvas = TRUE)) %>%
    addProviderTiles(providers$CartoDB.Positron, group = "CartoDB.Positron") %>%
    addProviderTiles(providers$Esri.WorldTopoMap, group = "Esri.WorldTopoMap") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Esri.WorldImagery") %>%
    addPolygons(
      fillColor = ~pal_rate(covid_rate),
      color = "#585858", weight = 0.7, smoothFactor = 0.3,
      fillOpacity = 0.8,
      label = ~sprintf("%s: %s per 100k", label_name, comma(covid_rate, accuracy = 0.1)),
      popup = popup_template,
      group = "COVID-19 cumulative rate"
    ) %>%
    addLegend("topright", pal = pal_rate, values = joined_spatial_wgs$covid_rate, title = "COVID-19 cumulative rate (per 100k)", opacity = 0.9, group = "COVID-19 cumulative rate") %>%
    addPolygons(
      fillColor = ~pal_cluster(local_cluster),
      color = "#444444", weight = 0.6, smoothFactor = 0.3,
      fillOpacity = 0.8,
      label = ~sprintf("%s: %s (p = %.3f)", label_name, local_cluster, local_moran_p),
      popup = ~sprintf("<b>%s</b><br/>Local Moran's I cluster: %s<br/>p-value: %.3f", label_name, local_cluster, local_moran_p),
      group = "Local Moran's I clusters"
    ) %>%
    addLegend("topright", pal = pal_cluster, values = joined_spatial_wgs$local_cluster, title = "Local Moran's I clusters", opacity = 0.9, group = "Local Moran's I clusters") %>%
    addPolygons(
      fillColor = ~pal_nohs(perc_no_hs_diploma),
      color = "#555", weight = 0.6, smoothFactor = 0.3,
      fillOpacity = 0.8,
      label = ~sprintf("%s: %0.1f%% without high school diploma", label_name, perc_no_hs_diploma),
      popup = popup_template,
      group = "Percent without high school diploma"
    ) %>%
    addLegend("topright", pal = pal_nohs, values = joined_spatial_wgs$perc_no_hs_diploma, title = "Education disadvantage (percentage without high school diploma)", opacity = 0.9, group = "Percent without high school diploma") %>%
    addPolygons(
      fillColor = ~pal_newcomer(perc_newcomers),
      color = "#555", weight = 0.6, smoothFactor = 0.3,
      fillOpacity = 0.8,
      label = ~sprintf("%s: %0.1f%% recent immigrants", label_name, perc_newcomers),
      popup = popup_template,
      group = "Recent immigrants"
    ) %>%
    addLegend("topright", pal = pal_newcomer, values = joined_spatial_wgs$perc_newcomers, title = "Recent immigrants (percentage)", opacity = 0.9, group = "Recent immigrants") %>%
    addPolygons(
      fillColor = ~pal_transit(perc_public_transport_work),
      color = "#555", weight = 0.6, smoothFactor = 0.3,
      fillOpacity = 0.8,
      label = ~sprintf("%s: %0.1f%% transit commuters", label_name, perc_public_transport_work),
      popup = popup_template,
      group = "Percent transit commuters"
    ) %>%
    addLegend("topright", pal = pal_transit, values = joined_spatial_wgs$perc_public_transport_work, title = "Percentage commuting by transit", opacity = 0.9, group = "Percent transit commuters") %>%
    addPolygons(
      fillColor = ~pal_seniors(perc_seniors),
      color = "#555", weight = 0.6, smoothFactor = 0.3,
      fillOpacity = 0.8,
      label = ~sprintf("%s: %0.1f%% age 65 and over", label_name, perc_seniors),
      popup = popup_template,
      group = "Percent age 65 and over"
    ) %>%
    addLegend("topright", pal = pal_seniors, values = joined_spatial_wgs$perc_seniors, title = "Percentage age 65 and over", opacity = 0.9, group = "Percent age 65 and over") %>%
    addPolygons(
      fillColor = ~pal_resid(resid_rate),
      color = "#555", weight = 0.6, smoothFactor = 0.3,
      fillOpacity = 0.8,
      label = ~sprintf("%s residual: %+0.2f", label_name, resid_rate),
      popup = popup_template,
      group = "Model residuals (observed - predicted)"
    ) %>%
    addLegend("topright", pal = pal_resid, values = joined_spatial_wgs$resid_rate, title = "Residual (observed - predicted)", opacity = 0.9, group = "Model residuals (observed - predicted)") %>%
    addPolygons(
      fillColor = ~pal_vuln(vulnerability_index),
      color = "#444", weight = 0.6, smoothFactor = 0.3,
      fillOpacity = 0.8,
      label = ~sprintf("%s: vulnerability index %0.0f", label_name, vulnerability_index),
      popup = popup_template,
      group = "Socioeconomic vulnerability (0–100)"
    ) %>%
    addLegend("topright", pal = pal_vuln, values = joined_spatial_wgs$vulnerability_index, title = "Socioeconomic vulnerability (0–100)", opacity = 0.9, group = "Socioeconomic vulnerability (0–100)") %>%
    addLayersControl(
      baseGroups = base_layers,
      overlayGroups = c(
        "COVID-19 cumulative rate",
        "Local Moran's I clusters",
        "Percent without high school diploma",
        "Recent immigrants",
        "Percent transit commuters",
        "Percent age 65 and over",
        "Model residuals (observed - predicted)",
        "Socioeconomic vulnerability (0–100)"
      ),
      options = layersControlOptions(collapsed = FALSE)
    ) %>%
    hideGroup(c(
      "Model residuals (observed - predicted)",
      "Percent without high school diploma",
      "Recent immigrants",
      "Percent transit commuters",
      "Percent age 65 and over"
    )) %>%
    showGroup(c("COVID-19 cumulative rate", "Local Moran's I clusters", "Socioeconomic vulnerability (0–100)")) %>%
    addScaleBar(position = "bottomleft", options = scaleBarOptions(imperial = FALSE, updateWhenIdle = TRUE)) %>%
    addMeasure(primaryLengthUnit = "meters", primaryAreaUnit = "sqmeters") %>%
    addMiniMap(tiles = providers$CartoDB.Positron, toggleDisplay = TRUE, minimized = TRUE) %>%
    addEasyButton(
      easyButton(
        icon = "fa-list",
        title = "Show/hide all legends",
        onClick = JS("
          function(btn, map){
            var legs = document.getElementsByClassName('legend');
            for (var i = 0; i < legs.length; i++) {
              legs[i].style.display = (legs[i].style.display === 'none') ? 'block' : 'none';
            }
          }
        ")
      )
    ) %>%
    addResetMapButton()

  htmlwidgets::saveWidget(
    map,
    file = file.path(outputs_dir, "interactive_map.html"),
    selfcontained = FALSE,
    libdir = file.path(outputs_dir, "interactive_map_libs")
  )
  message("Interactive map saved to: ", file.path(outputs_dir, "interactive_map.html"))
} else {
  message("Shapefile not found; spatial analysis skipped.")
}

# Random forest
set.seed(123)
rf_df <- analysis_data |>
  select(
    covid_rate_log,
    population_density,
    median_income_after_tax,
    unemployment_rate,
    perc_public_transport_work,
    perc_newcomers,
    avg_household,
    job_density,
    perc_no_hs_diploma,
    perc_bachelors_and_above,
    perc_seniors
  )

train_idx <- sample(seq_len(nrow(rf_df)), size = floor(0.8 * nrow(rf_df)))
rf_train <- rf_df[train_idx, ]
rf_test  <- rf_df[-train_idx, ]

rf_model <- randomForest(
  covid_rate_log ~ ., data = rf_train,
  ntree = 200, importance = TRUE, proximity = TRUE
)
rf_predictions <- predict(rf_model, rf_test)
rf_diff <- rf_test$covid_rate_log - rf_predictions
cat("Random forest MSE:", mean((rf_predictions - rf_test$covid_rate_log)^2), "\n")
print(shapiro.test(rf_diff))

var_imp <- importance(rf_model, type = 1)[, 1]
top_vars <- sort(var_imp, decreasing = TRUE)[1:10]
var_imp_plot <- ggplot(
  data.frame(Variable = factor(names(top_vars), levels = names(top_vars)), Importance = top_vars),
  aes(x = Variable, y = Importance)
) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 variable importance (Random Forest)", x = NULL, y = "Mean decrease accuracy")

ggsave(
  filename = file.path(outputs_dir, "rf_variable_importance.png"),
  plot = var_imp_plot,
  width = 7,
  height = 5,
  dpi = 150
)
