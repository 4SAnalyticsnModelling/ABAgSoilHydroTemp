from __future__ import annotations

"""
Shared utility functions for snow, depth parsing, and soil-surface coupling.

References used for equations in this project:

1. Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998).
   Crop evapotranspiration - Guidelines for computing crop water requirements.
   FAO Irrigation and Drainage Paper 56.
   https://www.fao.org/4/X0490E/X0490E00.htm

2. FAO training material on the Blaney-Criddle method:
   "The Blaney-Criddle formula: ETo = p (0.46 Tmean + 8)".
   https://www.fao.org/3/S2022E/s2022e07.htm

3. Alberta Energy Regulator DLS help:
   Alberta townships are numbered from the 49th parallel and are 6 miles
   north-south.
   https://dds.aer.ca/eas.net/webhelp/entering_dls_locations_.htm

4. Campbell, G.S. (1974). A Simple Method for Determining Unsaturated
   Conductivity from Moisture Retention Data. Soil Science 117(6):311-314.
   DOI: 10.1097/00010694-197406000-00001

5. Soil-temperature damping depth follows the 1-D heat diffusion solution
   where damping depth depends on thermal diffusivity and forcing period,
   e.g. Andujar Marquez et al. (2016), Sensors 16(3):306.
   DOI: 10.3390/s16030306

6. For soil thermal conductivity, this project uses a simplified geometric-mean
   mixing approach motivated by Johansen-style and de Vries-style soil thermal
   property models, with quartz-rich solids assumed more conductive than other
   minerals. See discussion and Johansen equations summarized in:
   https://pmc.ncbi.nlm.nih.gov/articles/PMC11283485/

"""

import numpy as np


ALBERTA_SOUTH_BORDER_LATITUDE_DEG = 49.0
TOWNSHIP_HEIGHT_MILES = 6.0
MILES_PER_LATITUDE_DEGREE = 69.172
HYDROLOGY_SUBSTEPS_PER_DAY = 24
DEFAULT_TARGET_DEPTHS = [
    7.5,
    15.0,
    30.0,
    45.0,
    60.0,
    75.0,
    90.0,
]
ICE_DENSITY_KG_M3 = 917.0
SNOW_DENSITY_MIN_KG_M3 = 80.0
SNOW_DENSITY_MAX_KG_M3 = 450.0
SNOW_THERMAL_CONDUCTIVITY_MIN_W_M_K = 0.08
SNOW_THERMAL_CONDUCTIVITY_MAX_W_M_K = 0.85
SNOW_EFFECTIVE_RESISTANCE_FACTOR = 0.9
SNOW_LAYER_COUNT = 5
SNOW_MIN_LAYER_SWE_MM = 1.0e-4
SNOW_OLD_LAYER_AGE_SCALE_DAYS = 10.0


def clamp(x: float, lo: float, hi: float) -> float:
    """Return ``x`` limited to the closed interval ``[lo, hi]``."""
    return max(lo, min(hi, x))


def compute_bulk_snow_density_kg_m3(
    swe_mm: float,
    snow_depth_mm: float,
) -> float:
    """Compute bulk snow density from snow-water equivalent and depth."""
    swe_mm = max(float(swe_mm), 0.0)
    snow_depth_mm = max(float(snow_depth_mm), 0.0)
    if snow_depth_mm <= 1.0e-9:
        return SNOW_DENSITY_MIN_KG_M3

    density = 1000.0 * (swe_mm / snow_depth_mm)
    return clamp(density, SNOW_DENSITY_MIN_KG_M3, SNOW_DENSITY_MAX_KG_M3)


def compute_snow_thermal_conductivity_w_m_k(
    snow_density_kg_m3: float,
    age_days: float = 0.0,
    fresh_fraction: float = 0.0,
) -> float:
    """Estimate bulk snow thermal conductivity from snow density and age."""
    snow_density_kg_m3 = clamp(
        float(snow_density_kg_m3),
        SNOW_DENSITY_MIN_KG_M3,
        SNOW_DENSITY_MAX_KG_M3,
    )
    age_days = max(float(age_days), 0.0)
    fresh_fraction = clamp(float(fresh_fraction), 0.0, 1.0)

    density_ratio = snow_density_kg_m3 / ICE_DENSITY_KG_M3
    conductivity = 0.021 + (2.5 * density_ratio**2)

    # Fresh low-density snow has weak bonds and is especially insulating.
    # Aging, settling, and bond growth gradually increase the conductivity.
    age_factor = 1.0 - np.exp(-age_days / SNOW_OLD_LAYER_AGE_SCALE_DAYS)
    conductivity *= (0.80 + 0.35 * age_factor)
    conductivity *= (1.0 - 0.18 * fresh_fraction)

    return clamp(
        conductivity,
        SNOW_THERMAL_CONDUCTIVITY_MIN_W_M_K,
        SNOW_THERMAL_CONDUCTIVITY_MAX_W_M_K,
    )


def initialize_snowpack_layers(
    swe_mm: float,
    snow_depth_mm: float,
    air_temp_c: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Initialize the conceptual multilayer snowpack state."""
    swe_mm = max(float(swe_mm), 0.0)
    snow_depth_mm = max(float(snow_depth_mm), 0.0)
    if swe_mm <= 1.0e-9 or snow_depth_mm <= 1.0e-9:
        return (
            np.zeros(SNOW_LAYER_COUNT, dtype=float),
            np.zeros(SNOW_LAYER_COUNT, dtype=float),
        )

    layer_swe_mm = np.full(SNOW_LAYER_COUNT, swe_mm / SNOW_LAYER_COUNT, dtype=float)
    # Unknown inherited pack is treated as pre-existing, mildly aged snow.
    layer_age_days = np.linspace(8.0, 20.0, SNOW_LAYER_COUNT, dtype=float)
    return layer_swe_mm, layer_age_days


def _regrid_snowpack_layers(
    layer_swe_mm: np.ndarray,
    layer_age_days: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Repartition snowpack SWE and age into the fixed conceptual layers."""
    layer_swe_mm = np.asarray(layer_swe_mm, dtype=float)
    layer_age_days = np.asarray(layer_age_days, dtype=float)
    valid = layer_swe_mm > SNOW_MIN_LAYER_SWE_MM
    if not np.any(valid):
        return (
            np.zeros(SNOW_LAYER_COUNT, dtype=float),
            np.zeros(SNOW_LAYER_COUNT, dtype=float),
        )

    swe = layer_swe_mm[valid]
    age = layer_age_days[valid]
    total_swe_mm = float(np.sum(swe))
    target_edges = np.linspace(0.0, total_swe_mm, SNOW_LAYER_COUNT + 1)

    new_swe = np.zeros(SNOW_LAYER_COUNT, dtype=float)
    new_age_mass = np.zeros(SNOW_LAYER_COUNT, dtype=float)

    cumulative = 0.0
    source_idx = 0
    remaining_in_source = float(swe[source_idx])

    for layer_idx in range(SNOW_LAYER_COUNT):
        target_amount = float(target_edges[layer_idx + 1] - target_edges[layer_idx])
        while target_amount > 1.0e-12 and source_idx < len(swe):
            moved = min(target_amount, remaining_in_source)
            new_swe[layer_idx] += moved
            new_age_mass[layer_idx] += moved * age[source_idx]
            target_amount -= moved
            remaining_in_source -= moved
            cumulative += moved
            if remaining_in_source <= 1.0e-12:
                source_idx += 1
                if source_idx < len(swe):
                    remaining_in_source = float(swe[source_idx])

    new_age = np.divide(
        new_age_mass,
        np.maximum(new_swe, SNOW_MIN_LAYER_SWE_MM),
        out=np.zeros_like(new_age_mass),
        where=new_swe > SNOW_MIN_LAYER_SWE_MM,
    )
    return new_swe, new_age


def compute_snowpack_bulk_properties(
    layer_swe_mm: np.ndarray | None,
    layer_age_days: np.ndarray | None,
    air_temp_c: float,
    surface_temp_c: float,
) -> tuple[float, float, float, float]:
    """Aggregate layered snow state into bulk SWE, depth, density, and conductivity."""
    if layer_swe_mm is None or layer_age_days is None:
        return 0.0, 0.0, SNOW_DENSITY_MIN_KG_M3, SNOW_THERMAL_CONDUCTIVITY_MIN_W_M_K

    layer_swe_mm = np.asarray(layer_swe_mm, dtype=float)
    layer_age_days = np.asarray(layer_age_days, dtype=float)
    if layer_swe_mm.size == 0 or np.sum(layer_swe_mm) <= 1.0e-9:
        return 0.0, 0.0, SNOW_DENSITY_MIN_KG_M3, SNOW_THERMAL_CONDUCTIVITY_MIN_W_M_K

    air_temp_c = float(air_temp_c)
    surface_temp_c = float(surface_temp_c)
    total_swe_mm = float(np.sum(layer_swe_mm))
    total_depth_m = 0.0
    total_resistance = 0.0
    total_age_mass = 0.0

    overburden_swe_mm = 0.0
    for swe_i, age_i in zip(layer_swe_mm, layer_age_days):
        if swe_i <= SNOW_MIN_LAYER_SWE_MM:
            continue

        fresh_density = compute_fresh_snow_density_kg_m3(air_temp_c=air_temp_c)
        age_factor = 1.0 - np.exp(-max(float(age_i), 0.0) / SNOW_OLD_LAYER_AGE_SCALE_DAYS)
        warm_factor = clamp((max(surface_temp_c, -15.0) + 15.0) / 20.0, 0.0, 1.0)
        overburden_factor = clamp(overburden_swe_mm / 120.0, 0.0, 1.0)
        target_old_density = 220.0 + 120.0 * warm_factor + 70.0 * overburden_factor
        density_i = fresh_density + age_factor * (target_old_density - fresh_density)
        density_i = clamp(density_i, SNOW_DENSITY_MIN_KG_M3, SNOW_DENSITY_MAX_KG_M3)

        fresh_fraction = np.exp(-max(float(age_i), 0.0) / 2.5)
        conductivity_i = compute_snow_thermal_conductivity_w_m_k(
            snow_density_kg_m3=density_i,
            age_days=float(age_i),
            fresh_fraction=float(fresh_fraction),
        )
        layer_depth_m = max((swe_i / 1000.0) * (1000.0 / max(density_i, 1.0)), 0.0)

        total_depth_m += layer_depth_m
        total_resistance += layer_depth_m / max(conductivity_i, 1.0e-6)
        total_age_mass += swe_i * age_i
        overburden_swe_mm += swe_i

    if total_depth_m <= 1.0e-12:
        return 0.0, 0.0, SNOW_DENSITY_MIN_KG_M3, SNOW_THERMAL_CONDUCTIVITY_MIN_W_M_K

    effective_conductivity = total_depth_m / max(total_resistance, 1.0e-12)
    bulk_density = 1000.0 * (total_swe_mm / max(total_depth_m * 1000.0, 1.0e-9))
    return (
        total_swe_mm,
        total_depth_m * 1000.0,
        clamp(bulk_density, SNOW_DENSITY_MIN_KG_M3, SNOW_DENSITY_MAX_KG_M3),
        clamp(
            effective_conductivity,
            SNOW_THERMAL_CONDUCTIVITY_MIN_W_M_K,
            SNOW_THERMAL_CONDUCTIVITY_MAX_W_M_K,
        ),
    )


def evolve_snowpack_layers(
    layer_swe_mm: np.ndarray | None,
    layer_age_days: np.ndarray | None,
    snowfall_water_mm: float,
    snowmelt_mm: float,
    dt_seconds: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Advance the conceptual snowpack by snowfall accumulation and melt removal."""
    snowfall_water_mm = max(float(snowfall_water_mm), 0.0)
    snowmelt_mm = max(float(snowmelt_mm), 0.0)
    dt_days = max(float(dt_seconds), 1.0) / 86400.0

    if layer_swe_mm is None or layer_age_days is None:
        layer_swe_mm = np.zeros(SNOW_LAYER_COUNT, dtype=float)
        layer_age_days = np.zeros(SNOW_LAYER_COUNT, dtype=float)
    else:
        layer_swe_mm = np.asarray(layer_swe_mm, dtype=float).copy()
        layer_age_days = np.asarray(layer_age_days, dtype=float).copy()

    if layer_swe_mm.size != SNOW_LAYER_COUNT:
        layer_swe_mm, layer_age_days = _regrid_snowpack_layers(layer_swe_mm, layer_age_days)

    layer_age_days = np.where(layer_swe_mm > SNOW_MIN_LAYER_SWE_MM, layer_age_days + dt_days, 0.0)

    if snowfall_water_mm > 1.0e-9:
        layer_swe_mm = np.concatenate(([snowfall_water_mm], layer_swe_mm))
        layer_age_days = np.concatenate(([0.0], layer_age_days))

    melt_remaining_mm = snowmelt_mm
    for idx in range(len(layer_swe_mm)):
        if melt_remaining_mm <= 1.0e-12:
            break
        removable = min(layer_swe_mm[idx], melt_remaining_mm)
        layer_swe_mm[idx] -= removable
        melt_remaining_mm -= removable

    return _regrid_snowpack_layers(layer_swe_mm, layer_age_days)


def compute_snow_surface_temperature_c(
    air_temp_c: float,
    top_soil_temp_c: float,
    swe_mm: float,
    snow_depth_mm: float,
    top_soil_thermal_conductivity_w_m_k: float,
    top_soil_half_thickness_m: float,
    snow_layer_swe_mm: np.ndarray | None = None,
    snow_layer_age_days: np.ndarray | None = None,
) -> float:
    """Estimate the effective snow or soil surface temperature."""
    air_temp_c = float(air_temp_c)
    top_soil_temp_c = float(top_soil_temp_c)
    top_soil_half_thickness_m = max(float(top_soil_half_thickness_m), 1.0e-6)
    top_soil_thermal_conductivity_w_m_k = max(float(top_soil_thermal_conductivity_w_m_k), 0.05)

    soil_resistance_m2_k_w = top_soil_half_thickness_m / top_soil_thermal_conductivity_w_m_k

    if snow_layer_swe_mm is not None and snow_layer_age_days is not None:
        layer_swe_mm = np.asarray(snow_layer_swe_mm, dtype=float)
        layer_age_days = np.asarray(snow_layer_age_days, dtype=float)
        pack_swe_mm, pack_depth_mm, _, pack_k = compute_snowpack_bulk_properties(
            layer_swe_mm=layer_swe_mm,
            layer_age_days=layer_age_days,
            air_temp_c=air_temp_c,
            surface_temp_c=top_soil_temp_c,
        )
        swe_mm = pack_swe_mm
        snow_depth_mm = pack_depth_mm
        snow_thermal_conductivity_w_m_k = pack_k
    else:
        snow_depth_mm = max(float(snow_depth_mm), 0.0)
        if snow_depth_mm <= 1.0e-9:
            return air_temp_c
        snow_density_kg_m3 = compute_bulk_snow_density_kg_m3(swe_mm=swe_mm, snow_depth_mm=snow_depth_mm)
        snow_thermal_conductivity_w_m_k = compute_snow_thermal_conductivity_w_m_k(
            snow_density_kg_m3=snow_density_kg_m3
        )

    snow_depth_m = max(float(snow_depth_mm), 0.0) / 1000.0
    if snow_depth_m <= 1.0e-9 or swe_mm <= 1.0e-9:
        return air_temp_c

    # Daily forcing and bulk snow depth overstate insulation when wind pumping,
    # patchiness, and unresolved ventilation are ignored.
    snow_resistance_m2_k_w = (
        SNOW_EFFECTIVE_RESISTANCE_FACTOR
        * snow_depth_m
        / max(snow_thermal_conductivity_w_m_k, 1.0e-6)
    )
    total_resistance_m2_k_w = snow_resistance_m2_k_w + soil_resistance_m2_k_w
    if total_resistance_m2_k_w <= 1.0e-9:
        return air_temp_c

    heat_flux_w_m2 = (air_temp_c - top_soil_temp_c) / total_resistance_m2_k_w
    surface_temp_c = air_temp_c - (heat_flux_w_m2 * snow_resistance_m2_k_w)
    return surface_temp_c


def compute_fresh_snow_density_kg_m3(air_temp_c: float) -> float:
    """Estimate fresh-snow density from air temperature."""
    air_temp_c = float(air_temp_c)
    density = 67.9 + (51.3 * np.exp(air_temp_c / 2.6))
    return clamp(density, SNOW_DENSITY_MIN_KG_M3, 220.0)


def evolve_snowpack_depth_mm(
    swe_mm: float,
    snow_depth_mm: float,
    snowfall_water_mm: float,
    snowmelt_mm: float,
    air_temp_c: float,
    surface_temp_c: float,
    dt_seconds: float,
) -> tuple[float, float, float]:
    """Update bulk snow depth and density after snowfall and melt."""
    swe_mm = max(float(swe_mm), 0.0)
    snow_depth_mm = max(float(snow_depth_mm), 0.0)
    snowfall_water_mm = max(float(snowfall_water_mm), 0.0)
    snowmelt_mm = max(float(snowmelt_mm), 0.0)
    dt_seconds = max(float(dt_seconds), 1.0)
    air_temp_c = float(air_temp_c)
    surface_temp_c = float(surface_temp_c)

    if swe_mm <= 1.0e-9 and snowfall_water_mm <= 1.0e-9:
        return 0.0, 0.0, SNOW_DENSITY_MIN_KG_M3

    current_density_kg_m3 = compute_bulk_snow_density_kg_m3(swe_mm=swe_mm, snow_depth_mm=snow_depth_mm)
    current_depth_mm = 1000.0 * swe_mm / current_density_kg_m3 if swe_mm > 1.0e-9 else 0.0
    fresh_density_kg_m3 = compute_fresh_snow_density_kg_m3(air_temp_c=air_temp_c)
    fresh_snow_depth_mm = 1000.0 * snowfall_water_mm / fresh_density_kg_m3

    swe_after_snow_mm = swe_mm + snowfall_water_mm
    depth_after_snow_mm = current_depth_mm + fresh_snow_depth_mm
    if swe_after_snow_mm <= 1.0e-9 or depth_after_snow_mm <= 1.0e-9:
        return 0.0, 0.0, SNOW_DENSITY_MIN_KG_M3

    pre_melt_density_kg_m3 = compute_bulk_snow_density_kg_m3(
        swe_mm=swe_after_snow_mm,
        snow_depth_mm=depth_after_snow_mm,
    )
    final_swe_mm = max(swe_after_snow_mm - snowmelt_mm, 0.0)
    if final_swe_mm <= 1.0e-9:
        return 0.0, 0.0, SNOW_DENSITY_MIN_KG_M3

    final_depth_mm = 1000.0 * final_swe_mm / pre_melt_density_kg_m3
    final_density_kg_m3 = compute_bulk_snow_density_kg_m3(
        swe_mm=final_swe_mm,
        snow_depth_mm=final_depth_mm,
    )
    return final_swe_mm, final_depth_mm, final_density_kg_m3


def normalize_depth(depth_value: object) -> str:
    """Normalize a depth token from the input files."""
    return str(depth_value).strip().replace(" ", "")


def parse_depth_to_cm(depth_value: object) -> tuple[float, float]:
    """Parse a depth interval like ``'15-25'`` into top and bottom depth in cm."""
    s = normalize_depth(depth_value)
    a, b = s.split("-")
    return float(a), float(b)


def parse_bottom_depth_to_cm(depth_value: object) -> float:
    """Parse a depth token into a layer bottom depth in cm."""
    s = normalize_depth(depth_value)
    if "-" in s:
        _, bottom_cm = parse_depth_to_cm(s)
        return bottom_cm
    return float(s)


def format_bottom_depth_cm(depth_cm: float) -> str:
    """Format a bottom depth for stable column naming."""
    depth_cm = float(depth_cm)
    if depth_cm.is_integer():
        return str(int(depth_cm))
    return f"{depth_cm:g}"


def depth_to_safe_name(depth_cm: float) -> str:
    """Convert a numeric depth into a filesystem- and column-safe token."""
    return format_bottom_depth_cm(depth_cm).replace(".", "p")


def parse_weather_dates(date_series: object) -> np.ndarray:
    """Parse supported daily weather date formats into pandas timestamps."""
    import pandas as pd

    supported_formats = ["%d-%b-%y", "%d-%m-%Y", "%Y-%m-%d"]
    last_error = None

    for date_format in supported_formats:
        try:
            return pd.to_datetime(date_series, format=date_format)
        except (TypeError, ValueError) as exc:
            last_error = exc

    raise ValueError(
        "Could not parse daily_weather.txt dates. Supported formats are "
        "'dd-Mon-yy', 'dd-mm-YYYY', and 'YYYY-mm-dd'."
    ) from last_error


def depth_contains(
    outer_top: float,
    outer_bottom: float,
    inner_top: float,
    inner_bottom: float,
    tol: float,
) -> bool:
    """Return whether an outer layer fully contains an inner layer."""
    return (outer_top <= inner_top + tol) and (outer_bottom >= inner_bottom - tol)
