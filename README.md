# Spatial Patterns and Socioeconomic Determinants of COVID-19 in Ottawa

Interactive GIS workflow (R + Leaflet) to map cumulative COVID-19 infection rates across Ottawa neighbourhoods, relate them to socioeconomic determinants, and expose modelled patterns through an explorable web map.

---

## What this project does
- Cleans and joins COVID-19 rates with neighbourhood determinants.
- Runs stepwise OLS, SAR error/lag/mixed models, GWR/MGWR, and a random forest baseline.
- Computes global and local Moran’s I (hot/cold spots and outliers).
- Publishes an interactive Leaflet map with togglable layers for observed rates, clusters, key determinants, predictions, residuals, and a composite vulnerability index.

---

## Inputs (as used in the script)
- `ONSMapCOVID_EN.csv` — cumulative COVID-19 rates by ONS neighbourhood (columns `ONS ID`, `ONS Neighbourhood Name`, `Cumulative Rate Excluding Cases Linked to Outbreaks in LTCH & RH`).
- `COVID DETERMINANTS IN ARC.xlsx` — socioeconomic determinants table (auto-copied from `covid_determinants_arc.xlsx` if not present).
- Shapefile: `Shape file/shapefile/ONS_Boundaries_Gen2.shp` (found automatically; any ONS Gen2 boundary shapefile works).
- Config: `config.R` points `data_dir` to the project folder by default.

### Determinant fields
`population_density`, `median_income_after_tax`, `unemployment_rate`, `perc_public_transport_work`, `perc_newcomers`, `avg_household`, `job_density`, `perc_no_hs_diploma`, `perc_bachelors_and_above`, `perc_seniors` (plus log/sqrt transforms created in-script).

---

## How to run
```powershell
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" Covid_analysis.r
```
Outputs are written to `outputs/` inside `data_dir`. Shapefile is required for spatial statistics and the web map; without it, spatial steps are skipped.

---

## Outputs generated
- Tables: `covid_determinants_arc.xlsx` (cleaned join for GIS).
- Static maps: `outputs/covid_rate_log_map.png`, `outputs/local_moran_clusters.png`, `outputs/ols_residuals_map.png`.
- Diagnostics: `outputs/correlation_heatmap.png`, `outputs/covid_rate_hist.png`, `outputs/moran_scatterplot.png`, `outputs/ols_diagnostics.png`, `outputs/ols_residual_hist.png`, `outputs/rf_variable_importance.png`.
- Interactive web map: `outputs/interactive_map.html` (+ assets in `outputs/interactive_map_libs/`).

---

## Interactive Leaflet map (layers & controls)
- Base layers: CartoDB.Positron, Esri.WorldTopoMap, Esri.WorldImagery.
- Overlays:
  - **COVID rate** (choropleth).
  - **Local Moran clusters** (High-High, Low-Low, High-Low, Low-High, Not significant).
  - **% no HS diploma**, **% transit commuters**, **% age 65+**.
  - **Predicted COVID rate** (SAR mixed predictions).
  - **Model residuals** (observed – predicted).
  - **Composite vulnerability index** (z-scored risk factors + inverse income, rescaled 0–100).
- Popups: neighbourhood name, observed rate, predicted rate, residual, determinants, vulnerability index.
- Controls: layer toggle, legend per layer, reset button, scale bar; polygons smoothed for cleaner rendering.

Open the map by double-clicking `outputs/interactive_map.html` (libs folder must stay beside it).

---

## Models and spatial stats
- Stepwise OLS on log(covid_rate) with key determinants.
- Spatial error, spatial lag, and spatial Durbin mixed models (spatialreg).
- Local Moran’s I for cluster/outlier detection (spdep).
- GWR and MGWR (GWmodel) for spatially varying effects.
- Random forest for variable importance benchmark.



