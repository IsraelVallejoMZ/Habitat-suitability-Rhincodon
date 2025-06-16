# Projected Distribution of *Rhincodon typus* Under Climate Change (SSP5-8.5, 2050)

This project models the ecological niche of the whale shark (*Rhincodon typus*) under current and future oceanographic conditions using MaxNet in R. It evaluates habitat stability and shifts under climate scenario SSP5-8.5.

## Objectives

- Download and clean GBIF occurrence records for *R. typus*.
- Integrate environmental variables from Bio-ORACLE (current and projected 2050).
- Train and evaluate a MaxNet model using presence-only data.
- Project future habitat suitability and quantify change.

## Project Structure

```
Rhincodon_typus_Distribution_2050/
├── sdm_whaleshark.R
├── data/
│   ├── env_layers_current/
│   └── env_layers_future_2050/
├── output/
│   └── report.pdf
├── maps/
│   └── distribution_maps.png
```

## Requirements

```r
install.packages(c("rgbif", "dplyr", "sf", "terra", "maxnet", "ENMeval", "geodata", "ggplot2"))
```

## How to Run

1. Run `sdm_whaleshark.R` to preprocess data, train model, and visualize results.
2. Outputs are saved in `/maps` and `/output`.
