# ABAgSoilHydroTemp

`ABAgSoilHydroTemp` is a daily time-step soil water and soil temperature model for Alberta agricultural site-year-crop simulations. The model couples weather forcing, crop-stage evapotranspiration demand, root water uptake, vertical soil-water redistribution, snow accumulation and melt, freeze-thaw partitioning of soil water, and layered soil-temperature evolution.

The codebase is implemented as a Python package in [`src/ab_ag_soil_hydro_temp`](src/ab_ag_soil_hydro_temp) and is designed to run from plain text input tables into a single publication-ready output table.

The rendered scientific-methodology document that matches the current code is available at [`docs/scientific_methodology.md`](docs/scientific_methodology.md).

## Model Workflow

The model workflow is organized across the package modules as follows:

- [`main.py`](src/ab_ag_soil_hydro_temp/main.py): top-level orchestration, daily simulation loop, annual lower-boundary thermal forcing, carry-over of previous-year state, and output assembly
- [`io.py`](src/ab_ag_soil_hydro_temp/io.py): validation and ingestion of input files, soil-profile expansion to model layers, and output column ordering
- [`crop.py`](src/ab_ag_soil_hydro_temp/crop.py): Alberta township to latitude conversion, daylength, Blaney-Criddle daylight term, crop-stage timing, crop coefficient interpolation, and root weighting
- [`water.py`](src/ab_ag_soil_hydro_temp/water.py): hydraulic relationships, van Genuchten and Mualem functions, freeze-thaw water partitioning, snowmelt partitioning, Richards-type redistribution, and root water uptake
- [`heat.py`](src/ab_ag_soil_hydro_temp/heat.py): bulk thermal properties, enthalpy-temperature conversion, and daily soil heat diffusion
- [`common.py`](src/ab_ag_soil_hydro_temp/common.py): shared constants, snowpack utilities, surface temperature coupling, and depth parsing helpers

At a high level, each `site-year-crop` simulation proceeds as:

1. Read and validate weather, soil, crop, boundary, and layer-definition inputs.
2. Expand source soil horizons onto the model target layers.
3. Initialize layer water storage, snowpack state, and soil temperature.
4. For each day:
   - compute air-temperature metrics and crop demand
   - partition precipitation into rain and snow, and update snowpack
   - solve coupled soil-water redistribution, root uptake, and soil heat diffusion
   - write the daily state and diagnostic variables
5. Carry the end-of-year soil moisture, soil temperature, and snow depth to the next year for the same site.

## Running the Model

If the package is installed, the configured console script is:

```bash
ab-ag-soil-hydro-temp
```

Within this repository, the current workflow also works with:

```bash
python -m src.ab_ag_soil_hydro_temp
```

The default run reads the input files from [`data`](data) and writes the main output file to `outputs/ab_ag_soil_hydro_temp_model_outputs.txt`.

## Input Files

The model expects five input tables in the repository [`data`](data) directory.

### `daily_weather.txt`

Daily meteorological forcing for each `site-year-crop`.

| Column | Description |
| --- | --- |
| `site` | Site name |
| `year` | Simulation year |
| `crop` | Crop name |
| `date` | Daily date; supported formats are `dd-Mon-yy`, `dd-mm-YYYY`, and `YYYY-mm-dd` |
| `precp` | Daily precipitation, mm |
| `daily_max_air_temp` | Daily maximum air temperature, °C |
| `daily_min_air_temp` | Daily minimum air temperature, °C |

Notes:

- The model rejects rows where `daily_max_air_temp < daily_min_air_temp`.
- The weather table defines which `site-year-crop` combinations are simulated.

### `soil_data.txt`

Layered soil hydraulic and initial-state information for each `site-year-crop`.

| Column | Description |
| --- | --- |
| `site` | Site name |
| `year` | Simulation year |
| `crop` | Crop name |
| `depth` | Layer bottom depth or depth interval, cm |
| `porosity` | Saturated volumetric water content, m3 m-3 |
| `field_capacity` | Volumetric water content at field capacity, m3 m-3 |
| `wilting_point` | Volumetric water content at wilting point, m3 m-3 |
| `ksat_mm_h` | Saturated hydraulic conductivity, mm h-1 |
| `initial_soil_moisture` | Initial volumetric water content, m3 m-3 |
| `inflection_point_negMPa` | Retention-curve inflection suction magnitude, MPa |

Notes:

- In this file, `-1.0` is a sentinel, not a physical soil value.
- For `porosity`, `field_capacity`, `wilting_point`, `inflection_point_negMPa`, and `ksat_mm_h`, a negative value means "reuse the previously resolved value for the same site and depth".
- For `initial_soil_moisture`, a negative value means "seed this year from the previous available site-depth state". The code first uses the previous year's simulated end-of-year soil moisture if available; otherwise it falls back to the most recent non-negative input seed for that site and depth.
- `ksat_mm_h` has one additional first-year rule. If every `ksat_mm_h` value for a site-depth sequence is negative across the full record, the model estimates the first resolved `ksat_mm_h` from `porosity`, `field_capacity`, `wilting_point`, and `inflection_point_negMPa`, then carries that resolved value forward in later years.
- The model derives `theta_residual`, `van_genuchten_alpha_per_mpa`, `van_genuchten_n`, and `van_genuchten_m` from these inputs before simulation.

### `crop.txt`

Crop stage timing and crop coefficient information for each `site-year-crop`.

| Column | Description |
| --- | --- |
| `site` | Site name |
| `year` | Simulation year |
| `crop` | Crop name |
| `alberta_township` | Alberta township number used to estimate latitude |
| `season_start_date` | Start of crop season |
| `emergence_date` | Emergence date |
| `anthesis_date` | Anthesis date |
| `early_grain_fill_date` | Early grain-fill date |
| `maturity_date` | Maturity date |
| `harvest_date` | Harvest date |
| `kc_off_season` | Off-season crop coefficient |
| `kc_initial` | Initial crop coefficient |
| `kc_mid` | Mid-season crop coefficient |
| `kc_end` | End-season crop coefficient |

Notes:

- Stage dates are interpreted within the simulation year.
- Stage dates must be chronological for each `site-year-crop`.

### `boundary.txt`

Site-level boundary-condition parameters.

| Column | Description |
| --- | --- |
| `site` | Site name |
| `drainage` | Deep-drainage impedance factor, 0 to 1 |
| `initial_snowpack_depth_in_mm` | Initial snow depth for the first simulated year at the site, mm |

### `layer.txt`

Target output layers used by the model.

| Column | Description |
| --- | --- |
| `depth` | Layer bottom depth, cm |

Notes:

- Source soil horizons from `soil_data.txt` are expanded or mapped onto these target layers before the simulation starts.
- The current repository configuration uses target depths `1, 9, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105 cm`.

## Output File

The main output is written to:

- [`outputs/ab_ag_soil_hydro_temp_model_outputs.txt`](outputs/ab_ag_soil_hydro_temp_model_outputs.txt)

This file contains one row per simulated `site-year-crop-date`.

## Output Parameters

The output columns fall into two groups:

- scalar daily fields that apply to the whole profile or site-day
- repeated depth-specific fields for each target layer

### Scalar Daily Output Fields

| Column | Description |
| --- | --- |
| `site` | Site name |
| `year` | Simulation year |
| `crop` | Crop name |
| `date` | Simulation date |
| `alberta_township` | Township number from crop input |
| `latitude_deg` | Estimated site latitude, degrees |
| `daily_max_air_temp_c` | Daily maximum air temperature, °C |
| `daily_min_air_temp_c` | Daily minimum air temperature, °C |
| `air_temp_c` | Daily mean air temperature, °C |
| `daily_air_temp_range_c` | Daily air temperature range, °C |
| `surface_diurnal_thaw_bias_c` | Cold-condition thaw-bias term used in surface forcing, °C |
| `soil_surface_forcing_temp_c` | Effective surface temperature applied to the thermal solver, °C |
| `deep_boundary_temp_c` | Lower-boundary temperature forcing, °C |
| `lower_boundary_reference_depth_cm` | Effective lower-boundary reference depth, cm |
| `precip_mm` | Daily precipitation, mm |
| `rain_mm` | Daily rainfall portion of precipitation, mm |
| `snowfall_water_mm` | Daily snowfall water equivalent, mm |
| `snowmelt_mm` | Daily snowmelt, mm |
| `effective_precip_mm` | Rainfall plus snowmelt reaching the soil surface, mm |
| `snow_depth_mm` | Bulk snow depth, mm |
| `snow_water_equivalent_mm` | Bulk snow water equivalent, mm |
| `snow_density_kg_m3` | Bulk snow density, kg m-3 |
| `snowmelt_degree_c` | Degree-day thaw term used for melt, °C |
| `daylength_hours` | Daylength, h |
| `bc_p_factor_pct` | Daily Blaney-Criddle daylight percentage, % |
| `bc_temperature_term` | Blaney-Criddle temperature term |
| `rainfall_cloudiness_damping` | Rainfall damping factor used in ET forcing |
| `growth_stage` | Phenological stage label |
| `crop_coefficient` | Daily crop coefficient |
| `off_season_kc` | Off-season crop coefficient from input |
| `reference_et_mm` | Reference ET, mm |
| `seasonal_et_mm` | Crop-stage ET demand component, mm |
| `off_season_et_mm` | Off-season ET demand component, mm |
| `pet_mm` | Potential evapotranspiration, mm |
| `aet_mm` | Actual evapotranspiration, mm |
| `deep_drainage_impedence_factor` | Site drainage parameter from `boundary.txt` |
| `deep_drainage_mm` | Deep drainage from the bottom layer, mm |
| `coupled_solver_iterations` | Total nonlinear iterations used that day |
| `coupled_substeps_used` | Number of accepted substeps used that day |
| `newton_iterations_used` | Newton iterations accepted by the coupled solver |
| `picard_fallback_iterations_used` | Picard-style fallback iterations used |
| `damped_newton_iterations_used` | Damped Newton iterations used |
| `worst_substep_residual_norm` | Maximum nonlinear residual norm among accepted substeps |

### Depth-Specific Output Families

For each target depth token such as `1`, `9`, `15`, `25`, ..., `105`, the output contains the following repeated fields:

| Pattern | Description |
| --- | --- |
| `soil_temp_{depth}_c` | Soil temperature at the layer, °C |
| `theta_{depth}` | Total volumetric water content, m3 m-3 |
| `theta_liquid_{depth}` | Liquid volumetric water content, m3 m-3 |
| `theta_ice_{depth}` | Ice volumetric water content, m3 m-3 |
| `theta_residual_{depth}` | Residual volumetric water content used by the retention curve, m3 m-3 |
| `storage_{depth}_mm` | Total water storage in the layer, mm |
| `liquid_storage_{depth}_mm` | Liquid water storage in the layer, mm |
| `ice_storage_{depth}_mm` | Ice storage in the layer, mm |
| `storage_{depth}_at_porosity_mm` | Storage at porosity, mm |
| `storage_{depth}_at_field_capacity_mm` | Storage at field capacity, mm |
| `storage_{depth}_at_wilting_point_mm` | Storage at wilting point, mm |
| `storage_{depth}_at_residual_mm` | Storage at residual water content, mm |
| `plant_available_storage_{depth}_mm` | Water storage above wilting point, mm |
| `et_uptake_{depth}_mm` | Root uptake from the layer, mm |
| `freeze_thaw_damping_{depth}` | Freeze-thaw damping index used diagnostically in the thermal solution |
| `thermal_diffusivity_{depth}_m2_s` | Soil thermal diffusivity, m2 s-1 |
| `apparent_latent_heat_capacity_{depth}_j_m3_k` | Apparent latent heat capacity, J m-3 K-1 |
| `depressed_freezing_temp_{depth}_c` | Freezing-point depression temperature, °C |
| `ksat_{depth}_m_per_s` | Saturated hydraulic conductivity, m s-1 |
| `van_genuchten_alpha_{depth}_per_mpa` | van Genuchten alpha parameter, MPa-1 |
| `van_genuchten_n_{depth}` | van Genuchten n parameter |
| `van_genuchten_m_{depth}` | van Genuchten m parameter |

### Interlayer Flux Fields

For each layer except the deepest one, the output also includes:

| Pattern | Description |
| --- | --- |
| `vertical_flux_{depth}_to_next_mm` | Net vertical water flux from this layer to the next deeper layer, mm |
