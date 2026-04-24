# ABAgSoilHydroTemp Scientific Methodology

## Overview

`ABAgSoilHydroTemp` is a daily time-step, layered soil water and soil temperature model for agricultural sites in Alberta. The current implementation couples:

- daily weather forcing from `daily_weather.txt`
- crop phenology and stage-based crop coefficients from `crop.txt`
- layered soil hydraulic and initial-condition inputs from `soil_data.txt`
- site-specific drainage and initial snowpack settings from `boundary.txt`
- user-defined target output layers from `layer.txt`

The model is implemented in `src/ab_ag_soil_hydro_temp` and writes one daily output record per `site-year-crop-date`.

## Model Structure

### Spatial representation

The model runs on a one-dimensional soil column. Source soil layers provided in `soil_data.txt` are expanded onto the target layer bottoms specified in `layer.txt`. In the current repository, the target depths are:

`1, 9, 15, 25, 35, 45, 55, 65, 75, 85, 95, and 105 cm`

Each target layer stores:

- total volumetric water content
- liquid water and ice fractions
- soil temperature
- layer hydraulic parameters
- layer thermal properties

### Temporal representation

The model runs at a daily time step, but the coupled water-snow-temperature solution is internally sub-stepped over the day. This allows precipitation partitioning, snow accumulation and melt, liquid-water movement, evapotranspiration uptake, freeze-thaw partitioning, and heat diffusion to interact within the same daily update.

## Inputs and State Initialization

### Weather forcing

`daily_weather.txt` provides daily precipitation, daily maximum air temperature, and daily minimum air temperature for each `site-year-crop` combination. Dates are parsed from `dd-Mon-yy`, `dd-mm-YYYY`, or `YYYY-mm-dd`.

The weather file determines the simulation calendar. If a `site-year-crop` appears in weather but not in soil or crop inputs, the run is skipped for that combination and reported as a warning.

### Crop development and evapotranspiration

`crop.txt` provides:

- Alberta township, which is converted to approximate latitude
- crop stage dates from season start through harvest
- crop coefficients for off-season, initial, mid-season, and end-season conditions

The code converts township number to latitude using the Dominion Land Survey township spacing from the Alberta south border. Daylength is then computed astronomically from latitude and day of year.

Potential evapotranspiration is based on a daily Blaney-Criddle formulation:

`PET = Kc * p * (0.46 * Tmean + 8)`

where:

- `p` is the daily percentage of annual daytime hours
- `Tmean` is daily mean air temperature
- `Kc` is interpolated from the crop growth stage schedule

The crop coefficient is held low from seeding to emergence, increases to a peak by anthesis, stays high through early grain fill, declines toward maturity, and then reverts to the off-season value outside the crop season.

### Soil hydraulic inputs and the meaning of `-1.0`

`soil_data.txt` provides `porosity`, `field_capacity`, `wilting_point`, `ksat_mm_h`, `initial_soil_moisture`, and `inflection_point_negMPa` by `site-year-crop-depth`.

For this model, `-1.0` is not a physical soil value. It is a sentinel that means "reuse the previously resolved value for the same site and depth."

The current code uses that convention in two distinct ways:

1. For hydraulic descriptors:
   `porosity`, `field_capacity`, `wilting_point`, `inflection_point_negMPa`, and `ksat_mm_h`

   A negative value means the model carries forward the previous year's resolved hydraulic value for the same site and depth. This is handled in `io.py` before the simulation starts.

2. For initial soil water content:
   `initial_soil_moisture`

   A negative value means the model seeds the current year from the previous available state at that site and depth. The code first tries to use the previous year's simulated end-of-year soil water content. If that does not yet exist, it falls back to the most recent non-negative input seed encountered for that site-depth sequence.

This means `-1.0` in `initial_soil_moisture` is a carry-over instruction, not missing data in a generic sense.

### Special handling of `ksat_mm_h`

`ksat_mm_h` also uses `-1.0` as the carry-forward sentinel. In addition, the code contains a special first-year rule:

- if every `ksat_mm_h` entry for a given site-depth sequence is negative across all years, the first year's `ksat_mm_h` is estimated from `porosity`, `field_capacity`, `wilting_point`, and `inflection_point_negMPa`
- this estimate is obtained by fitting constrained van Genuchten parameters and then inverting the Mualem conductivity relation so that conductivity at field capacity is `0.1 mm day-1`
- that resolved `ksat_mm_h` is then reused for later years at the same site and depth

As a result, `-1.0` in `ksat_mm_h` means either:

- carry the previous resolved value forward, or
- for the first resolved year only, estimate `Ksat` from the other hydraulic anchors when the full site-depth series is negative

### Derived hydraulic parameters

The model does not require `theta_residual`, `alpha`, `n`, or `m` as direct inputs. Instead, the code derives them from:

- porosity
- field capacity
- wilting point
- inflection-point suction magnitude

Specifically:

- `theta_residual`
- `van_genuchten_alpha_per_mpa`
- `van_genuchten_n`
- `van_genuchten_m`

are solved once for the first resolved year at each site-depth and then reused in subsequent years.

## Water Balance Formulation

### Precipitation partitioning and snow

Daily precipitation is partitioned into rain and snowfall using air temperature and the snow-insulated surface temperature. Snow is represented as a conceptual five-layer pack. Fresh snow is added to the top of the pack, snow ages through time, and the pack conductivity and density evolve with age and overburden.

Snowmelt is computed with a degree-day approach. The melt factor varies with daylength, so the potential daily melt factor increases during longer days.

The water available to infiltrate the soil surface on a given day is:

`effective_precipitation = rainfall + snowmelt`

### Soil water redistribution

Vertical water movement is represented as a simplified one-dimensional Richards-type redistribution scheme. For each layer interface, the model:

- computes liquid-water suction using the van Genuchten retention curve
- computes unsaturated conductivity using the Mualem relation
- reduces conductivity when pore water is frozen
- limits fluxes by donor supply and receiver pore space

Deep drainage is removed from the bottom layer only when water exceeds field capacity, and the drainage rate is capped by both the unsaturated conductivity and the site drainage impedance factor from `boundary.txt`.

### Root water uptake

Potential transpiration demand is distributed through the profile using stage-dependent rooting depth and depth weighting. Actual uptake is constrained by:

- liquid water availability above wilting point
- moisture stress between wilting point and field capacity
- frozen-water limitation

This makes actual evapotranspiration smaller than potential evapotranspiration under dry or frozen soil conditions.

## Thermal Formulation

### Soil heat storage and diffusion

The temperature model is written in terms of volumetric enthalpy. Each layer has:

- sensible heat storage from the mineral matrix, liquid water, ice, and air
- latent heat associated with phase change

Thermal diffusivity is estimated from a simplified soil thermal conductivity formulation that depends on porosity, total water content, and ice content. The heat equation is then solved implicitly across the layered soil column each day.

### Freeze-thaw coupling

Freeze-thaw partitioning is based on the relation between capillary suction and freezing-point depression. For a given total water content and temperature, the code computes:

- liquid water content
- ice content
- depressed freezing temperature

This allows the model to represent latent heat buffering and the reduction in hydraulic conductivity associated with frozen pore space.

### Lower thermal boundary

The lower boundary temperature is not fixed as a constant. Instead, the model fits an annual harmonic to the daily mean air temperatures for the simulated year, damps and phase-shifts that harmonic with depth, and applies it at an effective boundary depth below the model column. This gives a seasonally varying lower boundary that is smoother and lagged relative to air temperature.

## Year-to-Year Carryover

The implementation carries the following states from one year to the next for the same site and crop:

- end-of-year soil water content by layer
- end-of-year soil temperature by layer
- end-of-year snow depth

This is important because the provided input tables are organized as sequential site-year simulations rather than isolated single years.

## Model Performance and Validation

The current implementation was tested against site data collected across 25 sites in Alberta from 2008 to 2012. Weather, soil temperature and moisture, and site properties such as soil texture and vegetation were measured as part of the provincial monitoring network documented through Alberta Climate Information Service:

<https://acis.alberta.ca/>

The evaluation plots in `docs/plots` compare modelled and measured soil moisture and soil temperature through time at multiple depths for each site. The model is intentionally simple: it uses a daily time step, a one-dimensional soil column, a Blaney-Criddle evapotranspiration formulation, simplified snow physics, and a practical hydraulic parameterization driven by a small set of soil descriptors. Even with that level of simplification, the model performs reasonably well as a screening and comparative tool across sites, seasons, and depths.

Performance should therefore be interpreted as:

- strong enough for broad site-year pattern reproduction, seasonal timing, and cross-site comparison
- weaker for event-scale processes, lateral flow, residue effects, preferential flow, and other processes not explicitly resolved

## Validation Figures

### Figure 1: Andrew AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 1: Andrew AGDM moisture](plots/andrew_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 2: Andrew AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 2: Andrew AGDM temperature](plots/andrew_agdm_soil_temperature_timeseries.png)

### Figure 3: Atmore AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 3: Atmore AGDM moisture](plots/atmore_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 4: Atmore AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 4: Atmore AGDM temperature](plots/atmore_agdm_soil_temperature_timeseries.png)

### Figure 5: Barnwell AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 5: Barnwell AGDM moisture](plots/barnwell_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 6: Barnwell AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 6: Barnwell AGDM temperature](plots/barnwell_agdm_soil_temperature_timeseries.png)

### Figure 7: Breton Plots comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 7: Breton Plots moisture](plots/breton_plots_model_vs_observed_moisture_timeseries.png)

### Figure 8: Breton Plots comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 8: Breton Plots temperature](plots/breton_plots_soil_temperature_timeseries.png)

### Figure 9: Brocket AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 9: Brocket AGDM moisture](plots/brocket_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 10: Brocket AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 10: Brocket AGDM temperature](plots/brocket_agdm_soil_temperature_timeseries.png)

### Figure 11: Champion AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 11: Champion AGDM moisture](plots/champion_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 12: Champion AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 12: Champion AGDM temperature](plots/champion_agdm_soil_temperature_timeseries.png)

### Figure 13: Cleardale AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 13: Cleardale AGDM moisture](plots/cleardale_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 14: Cleardale AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 14: Cleardale AGDM temperature](plots/cleardale_agdm_soil_temperature_timeseries.png)

### Figure 15: Del Bonita AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 15: Del Bonita AGDM moisture](plots/del_bonita_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 16: Del Bonita AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 16: Del Bonita AGDM temperature](plots/del_bonita_agdm_soil_temperature_timeseries.png)

### Figure 17: Fairview AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 17: Fairview AGDM moisture](plots/fairview_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 18: Fairview AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 18: Fairview AGDM temperature](plots/fairview_agdm_soil_temperature_timeseries.png)

### Figure 19: High Prairie AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 19: High Prairie AGDM moisture](plots/high_prairie_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 20: High Prairie AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 20: High Prairie AGDM temperature](plots/high_prairie_agdm_soil_temperature_timeseries.png)

### Figure 21: Hussar AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 21: Hussar AGDM moisture](plots/hussar_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 22: Hussar AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 22: Hussar AGDM temperature](plots/hussar_agdm_soil_temperature_timeseries.png)

### Figure 23: Killam AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 23: Killam AGDM moisture](plots/killam_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 24: Killam AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 24: Killam AGDM temperature](plots/killam_agdm_soil_temperature_timeseries.png)

### Figure 25: Lethbridge CDA comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 25: Lethbridge CDA moisture](plots/lethbridge_cda_model_vs_observed_moisture_timeseries.png)

### Figure 26: Lethbridge CDA comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 26: Lethbridge CDA temperature](plots/lethbridge_cda_soil_temperature_timeseries.png)

### Figure 27: Morrin AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 27: Morrin AGDM moisture](plots/morrin_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 28: Morrin AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 28: Morrin AGDM temperature](plots/morrin_agdm_soil_temperature_timeseries.png)

### Figure 29: Mundare AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 29: Mundare AGDM moisture](plots/mundare_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 30: Mundare AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 30: Mundare AGDM temperature](plots/mundare_agdm_soil_temperature_timeseries.png)

### Figure 31: Olds College AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 31: Olds College AGDM moisture](plots/olds_college_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 32: Olds College AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 32: Olds College AGDM temperature](plots/olds_college_agdm_soil_temperature_timeseries.png)

### Figure 33: Oliver AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 33: Oliver AGDM moisture](plots/oliver_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 34: Oliver AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 34: Oliver AGDM temperature](plots/oliver_agdm_soil_temperature_timeseries.png)

### Figure 35: Oyen AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 35: Oyen AGDM moisture](plots/oyen_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 36: Oyen AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 36: Oyen AGDM temperature](plots/oyen_agdm_soil_temperature_timeseries.png)

### Figure 37: Peoria AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 37: Peoria AGDM moisture](plots/peoria_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 38: Peoria AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 38: Peoria AGDM temperature](plots/peoria_agdm_soil_temperature_timeseries.png)

### Figure 39: Schuler AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 39: Schuler AGDM moisture](plots/schuler_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 40: Schuler AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 40: Schuler AGDM temperature](plots/schuler_agdm_soil_temperature_timeseries.png)

### Figure 41: Smoky Lake AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 41: Smoky Lake AGDM moisture](plots/smoky_lake_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 42: Smoky Lake AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 42: Smoky Lake AGDM temperature](plots/smoky_lake_agdm_soil_temperature_timeseries.png)

### Figure 43: Stettler AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 43: Stettler AGDM moisture](plots/stettler_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 44: Stettler AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 44: Stettler AGDM temperature](plots/stettler_agdm_soil_temperature_timeseries.png)

### Figure 45: Tomahawk AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 45: Tomahawk AGDM moisture](plots/tomahawk_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 46: Tomahawk AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 46: Tomahawk AGDM temperature](plots/tomahawk_agdm_soil_temperature_timeseries.png)

### Figure 47: Vermilion AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 47: Vermilion AGDM moisture](plots/vermilion_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 48: Vermilion AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 48: Vermilion AGDM temperature](plots/vermilion_agdm_soil_temperature_timeseries.png)

### Figure 49: Wrentham AGDM comparison of modelled vs. measured soil moisture at different depth over 2008-2012

![Figure 49: Wrentham AGDM moisture](plots/wrentham_agdm_model_vs_observed_moisture_timeseries.png)

### Figure 50: Wrentham AGDM comparison of modelled vs. measured soil temperature at different depths over 2008-2012

![Figure 50: Wrentham AGDM temperature](plots/wrentham_agdm_soil_temperature_timeseries.png)

## Notes on Scope

This methodology reflects the current Python implementation in `src/ab_ag_soil_hydro_temp`, not the earlier Word document. In particular, it matches the present code paths for:

- daily Blaney-Criddle evapotranspiration
- snowpack layering and degree-day melt
- van Genuchten and Mualem hydraulic functions
- freeze-thaw partitioning through freezing-point depression
- implicit daily heat diffusion in enthalpy form
- year-to-year carryover of soil water, soil temperature, snow depth, and negative-value input sentinels
