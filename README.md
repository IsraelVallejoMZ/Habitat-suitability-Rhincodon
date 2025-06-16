# Projected Distribution of *Rhincodon typus* Under Climate Change (SSP5-8.5, 2050)

This project models the ecological niche of the whale shark (*Rhincodon typus*) under current and future oceanographic conditions using MaxNet in R. It evaluates habitat stability and shifts under the SSP5-8.5 climate scenario.

## Objectives
- Download and clean GBIF occurrence records for *R. typus*.
- Integrate environmental variables from Bio-ORACLE (current and projected 2050).
- Train and evaluate a MaxNet model using presence-only data.
- Project future habitat suitability and quantify change.

## Environmental Layers   <!-- ⬅️ Nueva sección -->
Environmental predictors (sea-surface temperature, salinity, pH, chlorophyll-a, dissolved oxygen) were obtained from [Bio-ORACLE](https://www.bio-oracle.org) as `.nc` files for the year 2050 (scenario SSP5-8.5).

> **Note:** These rasters are large, so they are **not** included in this repository.  
> Download them manually here:
> - [Temperature 2050](https://www.bio-oracle.org/downloads/BO_envtemp_ssp585_2050.nc)
> - [Salinity 2050](https://www.bio-oracle.org/downloads/BO_salinity_ssp585_2050.nc)
> - [pH 2050](https://www.bio-oracle.org/downloads/BO_pH_ssp585_2050.nc)
> - [Chlorophyll-a 2050](https://www.bio-oracle.org/downloads/BO_chl_ssp585_2050.nc)
> - [Dissolved Oxygen 2050](https://www.bio-oracle.org/downloads/BO_o2_ssp585_2050.nc)

Place downloaded files in `data/env_layers_future_2050/` before running the script.

## Project Structure
Rhincodon_typus_Distribution_2050/
├── sdm_whaleshark.R
├── data/
│ ├── env_layers_current/
│ └── env_layers_future_2050/
├── output/
│ └── report.pdf
├── maps/
│ └── distribution_maps.png


## Requirements
```r
install.packages(c("rgbif", "dplyr", "sf", "terra", "maxnet", "ENMeval", "geodata", "ggplot2"))
How to Run
Run sdm_whaleshark.R to preprocess data, train the model, and generate outputs.

Maps and figures will be saved in /maps, and the final report in /report.
