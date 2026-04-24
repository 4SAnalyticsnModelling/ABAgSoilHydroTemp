# ABAgSoilHydroTemp Methodology and Performance Notes

## 1. Purpose and Scope

`ABAgSoilHydroTemp` is a one-dimensional daily time-step model for agricultural soils in Alberta. It simulates soil temperature, liquid water, ice content, evapotranspiration uptake, snowpack effects, vertical redistribution of water, and deep drainage for each site-year-crop combination. The model is process-based but intentionally practical, combining physically motivated equations with pragmatic simplifications suitable for data-limited agricultural applications.

## 2. Model Structure

The model workflow is organized across the package modules `main.py`, `io.py`, `crop.py`, `water.py`, `heat.py`, and `common.py`. Together these modules implement daily weather forcing, crop demand, snow accumulation and melt, freeze-thaw partitioning, soil-water redistribution, and soil heat diffusion within a vertically layered soil column.

## 3. Spatial Representation and Layers

The model is one-dimensional in the vertical direction. Each site is represented by a layered soil column. Source soil properties from `soil_data.txt` are mapped to target layers defined in `layer.txt`. Each target layer stores bottom depth, thickness, midpoint depth, hydraulic parameters, and dynamic thermal and hydrologic state variables.

## 4. Weather Forcing

Daily weather forcing is provided as precipitation, maximum air temperature, and minimum air temperature. The model computes mean daily air temperature and temperature range using the standard daily summaries shown below.

```text
T_mean = (T_max + T_min) / 2
```

where `T_mean` is mean daily air temperature (⁰C), `T_max` is daily maximum air temperature (⁰C), and `T_min` is daily minimum air temperature (⁰C).

```text
T_range = T_max - T_min
```

where `T_range` is daily air-temperature range (⁰C). The code also derives a cold-condition thaw-bias term from daily temperature range to reduce over-severe daily freezing in conditions with daytime thaw and nighttime cooling.

## 5. Latitude, Daylength, and Reference Evapotranspiration (ET)

Latitude is estimated from Alberta township number using Dominion Land Survey geometry. Daily daylength is then computed from standard FAO-56 astronomical relationships (Allen et al., 1998). Reference evapotranspiration is represented using a daily Blaney-Criddle form.

```text
ET_o = p * (0.46 * T_mean + 8)
```

where `ET_o` is reference evapotranspiration (mm d⁻¹), `p` is the daily percentage contribution of that day to annual daylight hours (percent or fraction converted consistently in the implementation), and `T_mean` is mean daily air temperature (⁰C). This follows the Blaney-Criddle approach described by FAO (Allen et al., 1998).

## 6. Crop Coefficient and Crop Stages

Potential evapotranspiration is computed by scaling reference ET with a daily crop coefficient `K_c`. The value of `K_c` changes with crop stage using dates and coefficients supplied in `crop.txt`.

```text
PET = K_c * ET_o
```

where `PET` is potential evapotranspiration (mm d⁻¹), `K_c` is crop coefficient (dimensionless), and `ET_o` is reference evapotranspiration (mm d⁻¹). The implemented stages are off season, seeding to emergence, emergence to anthesis, anthesis peak, anthesis to early grain fill, grain fill to maturity, and maturity to harvest.

## 7. Rooting Depth and Root Weights

Rooting depth is stage dependent, and root uptake weights are distributed with depth using an exponentially decaying active-root term plus a weaker background-access term. This means the model places more uptake capacity in shallow and intermediate layers but still allows some deep access when roots are established.

## 8. Soil Hydraulic Parameterization

The hydraulic framework is based on the van Genuchten retention model and Mualem conductivity model (van Genuchten, 1980; Mualem, 1976). The soil input file provides porosity, field capacity, wilting point, saturated hydraulic conductivity, initial soil moisture, and an inflection-point suction. From these values the model derives residual water content and the van Genuchten shape parameters.

The liquid-water retention equation is implemented in the following form.

```text
theta = theta_r + (theta_s - theta_r) * [1 + (alpha * |psi|)^n]^(-m)
```

where `theta` is liquid volumetric water content (m³ m⁻³), `theta_r` is residual volumetric water content (m³ m⁻³), `theta_s` is saturated volumetric water content or porosity (m³ m⁻³), `alpha` is the van Genuchten inverse air-entry parameter (MPa⁻¹), `psi` is matric suction magnitude (MPa), `n` is a shape parameter (dimensionless), and `m = 1 - 1/n` (dimensionless).

## 9. Freeze-Thaw Water Partitioning

A key feature of the model is that total water in each layer is partitioned into liquid water and ice depending on temperature and suction. The freezing-point depression logic follows the generalized Clapeyron concept used in freezing-soil models such as Dall'Amico et al. (2011).

For freezing-point depression, the model uses relative total water content rather than liquid-only saturation.

```text
S_total = (theta_total - theta_r) / (theta_s - theta_r)
```

where `S_total` is relative total water saturation (dimensionless), `theta_total` is total volumetric water content including liquid plus ice (m³ m⁻³), `theta_r` is residual volumetric water content (m³ m⁻³), and `theta_s` is porosity (m³ m⁻³).

The depressed freezing temperature is then computed from suction head `h` using a Clapeyron-type relation (Dall'Amico et al., 2011).

```text
T_f,dep = T_0 - (g * T_0 / L_f) * h
```

where `T_f,dep` is depressed freezing temperature (K), `T_0` is the freezing point of pure water (K), `g` is gravitational acceleration (m s⁻²), `L_f` is latent heat of fusion (J kg⁻¹), and `h` is suction head (m).

This freezing threshold is smoothed and converted into a state-dependent freezing suction before liquid water is obtained from the retention equation. Ice is then assigned as the remainder of total water.

```text
theta_ice = theta_total - theta_liquid
```

where `theta_ice` is volumetric ice content (m³ m⁻³) and `theta_liquid` is volumetric liquid water content (m³ m⁻³).

## 10. Unsaturated Hydraulic Conductivity

Hydraulic conductivity is based on the Mualem conductivity formulation (Mualem, 1976) combined with the van Genuchten effective saturation function (van Genuchten, 1980):

```text
K = K_sat * K_r(S_e) * f_ice
```

where `K` is unsaturated hydraulic conductivity (m s⁻¹), `K_sat` is saturated hydraulic conductivity (m s⁻¹), `K_r` is relative conductivity (dimensionless) as a function of effective saturation `S_e`, and `f_ice` is the ice-impedance factor. The code computes `K_r` with the Mualem-van Genuchten form and scales it by an exponential ice reduction term.

## 11. Precipitation Partitioning, Snow Accumulation, and Snowmelt

Daily precipitation is partitioned into rainfall and snowfall using a mixed phase temperature derived from air temperature and the snow-surface temperature. Snow is tracked through bulk depth, snow water equivalent, density, and a five-layer snowpack scheme. Snowmelt is represented with a degree-day-type relationship adjusted by daylength.

```text
M = F_m * D_thaw
```

where `M` is snowmelt (mm d⁻¹ or substep-equivalent mm), `F_m` is melt factor (mm ⁰C⁻¹ d⁻¹), and `D_thaw` is thaw-degree term (⁰C). The effective liquid input reaching the soil is:

```text
P_eff = P_rain + M
```

where `P_eff` is effective precipitation (mm), `P_rain` is rainfall (mm), and `M` is snowmelt (mm).

## 12. Snow Insulation and Surface Temperature Coupling

Snow reduces thermal coupling between air and the topsoil. The model estimates snow thermal conductivity from density and age, converts the snowpack to a thermal resistance, and combines that with the upper-soil resistance. Deeper and less conductive snow insulates the soil and buffers it from atmospheric cold.

## 13. Vertical Soil-Water Redistribution

The model resolves water redistribution using multiple substeps per day. For each adjacent pair of layers it computes a hydraulic gradient from gravity and suction difference and then applies a conductivity-limited flux. A simplified Richards-type relation is used.

```text
nablaH = 1 + (psi_lower - psi_upper) / delta_z
```

where `nablaH` is total hydraulic gradient (dimensionless), `psi_lower` and `psi_upper` are matric suction heads of the lower and upper layers (m), and `delta_z` is distance between layer centers (m).

```text
q = K_interface * nablaH * delta_t
```

where `q` is vertical water flux over the timestep (m or converted to mm in the code), `K_interface` is interface hydraulic conductivity (m s⁻¹), `nablaH` is total hydraulic gradient (dimensionless), and `delta_t` is timestep length (s). The actual flux is then limited by donor liquid water, receiver pore space, and additional practical caps used to improve numerical stability.

## 14. Root Water Uptake

Root water uptake is computed from potential plant demand, root distribution, liquid water availability, and moisture-stress scaling. A layer contributes to root uptake only if roots are present, liquid water above wilting point is available, and freezing has not strongly reduced access to that water. The implementation allocates demand by root weights, scales it by relative extractable water and frozen fraction, and then redistributes unmet demand to layers that still have removable liquid water.

## 15. Thermal Properties

The thermal module represents volumetric sensible heat capacity as the sum of mineral, liquid-water, ice, and air contributions. This follows standard porous-media heat-capacity accounting.

```text
C = (1 - phi) * C_s + theta_l * C_w + theta_i * C_i + (phi - theta_total) * C_a
```

where `C` is volumetric sensible heat capacity (J m⁻³ K⁻¹), `phi` is porosity (m³ m⁻³), `C_s` is volumetric heat capacity of the mineral fraction (J m⁻³ K⁻¹), `C_w` is volumetric heat capacity of liquid water (J m⁻³ K⁻¹), `C_i` is volumetric heat capacity of ice (J m⁻³ K⁻¹), `C_a` is volumetric heat capacity of air (J m⁻³ K⁻¹), `theta_l` is liquid volumetric water content (m³ m⁻³), `theta_i` is ice volumetric water content (m³ m⁻³), and `theta_total` is total volumetric water content (m³ m⁻³).

Thermal diffusivity is then computed as:

```text
alpha_th = lambda / C
```

where `alpha_th` is thermal diffusivity (m² s⁻¹), `lambda` is bulk thermal conductivity (W m⁻¹ K⁻¹), and `C` is volumetric sensible heat capacity (J m⁻³ K⁻¹).

## 16. Enthalpy Formulation and Temperature Inversion

The thermal state is advanced through volumetric enthalpy. For each layer, the implementation uses the following enthalpy relation:

```text
H = C * T - L_v * theta_i
```

where `H` is volumetric enthalpy (J m⁻³), `C` is volumetric sensible heat capacity (J m⁻³ K⁻¹), `T` is soil temperature (K), `L_v` is volumetric latent heat of fusion (J m⁻³), and `theta_i` is volumetric ice content (m³ m⁻³).

Because `theta_i` depends on temperature through the freezing curve, converting between `H` and `T` is nonlinear. The model therefore uses an iterative inversion with a safeguarded Newton update and numerical estimate of `dH/dT`.

## 17. Daily Soil Heat Diffusion

The daily thermal solve is implicit in temperature space for stability, after which the updated temperature is converted back to enthalpy. The energy balance of each layer is:

```text
C * delta_z * (T_new - T_old) / delta_t = Q_in - Q_out
```

where `C` is volumetric sensible heat capacity (J m⁻³ K⁻¹), `delta_z` is layer thickness (m), `T_new` and `T_old` are new and old temperatures (⁰C), `delta_t` is timestep length (s), and `Q_in` and `Q_out` are conductive heat terms (W m⁻² converted consistently in the discretization).

## 18. Lower Thermal Boundary Condition

The lower thermal boundary is not fixed. Instead, the model fits a first annual harmonic to mean daily air temperature and damps it with depth. This follows the standard annual thermal-wave concept for soils (e.g., Andujar Marquez et al., 2016).

```text
T(z, t) = T_mean + b + D(z) * [A_cos * cos(theta - L(z)) + A_sin * sin(theta - L(z))]
```

where `T(z, t)` is boundary temperature at depth `z` and time `t` (⁰C), `T_mean` is annual mean air temperature (⁰C), `D(z)` is depth-dependent damping factor (dimensionless), `A_cos` and `A_sin` are harmonic coefficients (⁰C), `theta` is annual phase angle (radians), `L(z)` is depth-dependent phase lag (radians), and `b` is an empirical deep-boundary temperature offset (⁰C). In the code, `D(z) = exp(-z / z_d)` and `L(z) = z / z_d`, with a fixed calibrated damping depth and offset for the lower boundary evaluation.

## 19. Daily Solver Sequence

For each site-year-crop combination, the implemented daily sequence is: compute air-temperature metrics, compute daylength and crop demand, partition precipitation and update snowpack, run substeps of redistribution and root uptake, solve the thermal step, recompute freeze-thaw diagnostics, and write daily outputs.

## 20. Carry-Over Between Years

The model carries forward end-of-year soil moisture, soil temperature, and snow depth to initialize the next year for the same site. This means the simulations are linked across years and are not independent one-year runs unless the user resets those states manually in the input files.

## 21. Inputs and Outputs

Inputs include daily meteorological forcing, layer-wise soil hydraulic properties, crop-stage timing and crop coefficients, boundary drainage and initial snow depth, and the target output-layer structure. Outputs include daily soil temperature, total water content, liquid water, ice content, water storage, plant-available storage, evapotranspiration uptake by layer, interlayer vertical fluxes, deep drainage, thermal diffusivity, latent heat diagnostics, depressed freezing temperature, and fitted hydraulic parameters for each layer.

## 22. Model Performance

Model performance was assessed by comparing modeled and observed soil temperature and soil moisture time series at the same depths and dates. The most informative diagnostics are seasonal timing, freeze-up and thaw timing, amplitude damping with depth, and whether the model reproduces persistence of wet and dry periods. For soil temperature, performance metrics included root mean square error (⁰C), and coefficient of determination (R²). For soil moisture, index of model agreement (Willmott's d) and root mean square error (m³ m⁻³) were used as perfromance validation metrics. Graphical diagnostics were also used as summary metrics. 

The current implementation was tested against site data collected across 25 sites in Alberta from 2008 to 2012. Weather, soil temperature and moisture, and site properties such as soil texture and vegetation were measured as part of the provincial monitoring network led by Alberta Climate Information Service:

<https://acis.alberta.ca/>

The evaluation plots in `plots` compare modelled and measured soil moisture and soil temperature through time at multiple depths for each site. The model is intentionally simple: it uses a daily time step, a one-dimensional soil column, a Blaney-Criddle evapotranspiration formulation, simplified snow physics, and a practical hydraulic parameterization driven by a small set of soil descriptors. Even with that level of simplification, the model performs reasonably well as a screening and comparative tool across sites, seasons, and depths.

Performance should therefore be interpreted as:

- strong enough for broad site-year pattern reproduction, seasonal timing, and cross-site comparison
- weaker for event-scale processes, lateral flow, residue effects, preferential flow, and other processes not explicitly resolved

## 23. Practical Interpretation of the Model

The model should be understood as a process-informed simulation framework rather than a fully resolved mechanistic cryo-hydrologic model. It captures the dominant controls relevant to agricultural soils in a way that is computationally manageable and transparent. Its strongest value lies in linking crop demand, snow insulation, freeze-thaw water partitioning, and heat diffusion within one daily framework, while remaining simple enough for multi-site, multi-year application.

## 24. Validation Figures

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

## 25. References

Allen, R.G., Pereira, L.S., Raes, D., Smith, M. 1998. Crop evapotranspiration: Guidelines for computing crop water requirements. FAO Irrigation and Drainage Paper 56.

van Genuchten, M.T. 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils. Soil Science Society of America Journal 44: 892-898.

Mualem, Y. 1976. A new model for predicting the hydraulic conductivity of unsaturated porous media. Water Resources Research 12: 513-522.

Dall'Amico, M., Endrizzi, S., Gruber, S., Rigon, R. 2011. A robust and energy-conserving model of freezing variably-saturated soil. The Cryosphere 5: 469-484.

Campbell, G.S. 1974. A simple method for determining unsaturated conductivity from moisture retention data. Soil Science 117: 311-314.

Andujar Marquez, J.M., et al. 2016. Soil temperature prediction models and thermal damping-depth concepts for subsurface processes. Sensors 16(3):306.
