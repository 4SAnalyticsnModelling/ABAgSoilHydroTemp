from __future__ import annotations

"""Hydraulic, freeze-thaw, and snowmelt utilities for the coupled model.

References:
- van Genuchten (1980), closed-form soil-water retention equation.
- Mualem (1976), hydraulic conductivity model based on pore-size distribution.
- Dall'Amico et al. (2011), freezing soil formulation linking suction and phase
  partition through the generalized Clapeyron relation.
"""

import numpy as np

from .common import (
    clamp,
    compute_snow_surface_temperature_c,
    compute_snowpack_bulk_properties,
    evolve_snowpack_layers,
)


FREEZING_TEMPERATURE_K = 273.15
LATENT_HEAT_FUSION_J_KG = 3.34e5
GRAVITY_M_S2 = 9.80665
WATER_DENSITY_KG_M3 = 1000.0
ICE_IMPEDANCE_SHAPE_FACTOR = 8.0
FREEZE_CURVE_SMOOTHING_C = 0.15



def estimate_van_genuchten_parameters(
    theta_sat: np.ndarray | float,
    theta_fc: np.ndarray | float,
    theta_wp: np.ndarray | float,
    inflection_point_neg_mpa: np.ndarray | float,
    psi_fc_mpa: float = 0.033,
    psi_wp_mpa: float = 1.5,
    n_upper: float = 10.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Estimate van Genuchten theta_r, alpha, and n from theta_s, theta_fc, theta_wp,
    and the magnitude of the retention-curve inflection point in MPa.

    alpha is taken as the inverse of the inflection-point suction. For that fixed alpha,
    theta_r and n are selected to best satisfy the field-capacity and wilting-point
    anchors while enforcing n < n_upper.
    """
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_fc = np.asarray(theta_fc, dtype=float)
    theta_wp = np.asarray(theta_wp, dtype=float)
    inflection_point_neg_mpa = np.asarray(inflection_point_neg_mpa, dtype=float)

    theta_r_out = np.zeros_like(theta_sat, dtype=float)
    alpha_vg_out = np.zeros_like(theta_sat, dtype=float)
    n_vg_out = np.zeros_like(theta_sat, dtype=float)

    for idx in np.ndindex(theta_sat.shape):
        theta_s_i = max(float(theta_sat[idx]), 1.0e-6)
        theta_fc_i = float(theta_fc[idx])
        theta_wp_i = float(theta_wp[idx])
        inflection_i = max(float(inflection_point_neg_mpa[idx]), 1.0e-6)
        alpha_i = 1.0 / inflection_i

        n_grid = np.linspace(1.01, max(min(float(n_upper), 10.0), 1.02), 1500)
        m_grid = 1.0 - (1.0 / n_grid)

        se_fc = np.clip(
            (1.0 + (max(alpha_i, 1.0e-9) * max(psi_fc_mpa, 1.0e-9)) ** n_grid) ** (-m_grid),
            1.0e-8,
            0.999999,
        )
        se_wp = np.clip(
            (1.0 + (max(alpha_i, 1.0e-9) * max(psi_wp_mpa, 1.0e-9)) ** n_grid) ** (-m_grid),
            1.0e-8,
            0.999999,
        )
        theta_r_grid = (theta_fc_i - theta_s_i * se_fc) / np.maximum(1.0 - se_fc, 1.0e-9)
        theta_wp_pred = theta_r_grid + (theta_s_i - theta_r_grid) * se_wp

        valid = (
            (theta_r_grid >= 0.0)
            & (theta_r_grid < theta_wp_i - 1.0e-4)
            & (theta_r_grid < theta_fc_i - 1.0e-4)
            & (theta_fc_i < theta_s_i - 1.0e-5)
        )

        penalty = np.abs(theta_wp_pred - theta_wp_i)
        penalty += 1.0e3 * (~valid)
        best_idx = int(np.argmin(penalty))

        if not bool(valid[best_idx]):
            raise ValueError(
                "Could not estimate van Genuchten parameters from porosity, field capacity, "
                f"wilting point, and inflection point. Values were theta_sat={theta_s_i}, "
                f"theta_fc={theta_fc_i}, theta_wp={theta_wp_i}, inflection_point_neg_mpa={inflection_i}."
            )

        theta_r_out[idx] = float(theta_r_grid[best_idx])
        alpha_vg_out[idx] = float(alpha_i)
        n_vg_out[idx] = float(n_grid[best_idx])

    return theta_r_out, alpha_vg_out, n_vg_out


def tune_inflection_point_neg_mpa(
    theta_sat: float,
    theta_fc: float,
    theta_wp: float,
    initial_inflection_point_neg_mpa: float | None = None,
    psi_fc_mpa: float = 0.033,
    psi_wp_mpa: float = 1.5,
    n_upper: float = 10.0,
) -> tuple[float, float, float, float]:
    """
    Tune the inflection-point suction so the constrained van Genuchten fit remains
    physically consistent and satisfies n < n_upper.
    """
    if initial_inflection_point_neg_mpa is None or initial_inflection_point_neg_mpa <= 0.0:
        initial_inflection_point_neg_mpa = np.sqrt(max(psi_fc_mpa, 1.0e-6) * max(psi_wp_mpa, 1.0e-6))

    search_grid = np.unique(
        np.concatenate(
            (
                np.logspace(-3, np.log10(3.0), 240),
                np.asarray([initial_inflection_point_neg_mpa], dtype=float),
            )
        )
    )

    best = None
    for inflection in search_grid:
        try:
            theta_r, alpha_vg, n_vg = estimate_van_genuchten_parameters(
                theta_sat=theta_sat,
                theta_fc=theta_fc,
                theta_wp=theta_wp,
                inflection_point_neg_mpa=inflection,
                psi_fc_mpa=psi_fc_mpa,
                psi_wp_mpa=psi_wp_mpa,
                n_upper=n_upper,
            )
        except ValueError:
            continue

        theta_r_val = float(np.asarray(theta_r).item())
        alpha_val = float(np.asarray(alpha_vg).item())
        n_val = float(np.asarray(n_vg).item())
        m_val = float(1.0 - (1.0 / max(n_val, 1.000001)))

        se_wp = float((1.0 + (alpha_val * psi_wp_mpa) ** n_val) ** (-m_val))
        theta_wp_pred = theta_r_val + (theta_sat - theta_r_val) * se_wp
        wp_error = abs(theta_wp_pred - theta_wp)
        regularization = 1.0e-4 * abs(np.log(max(inflection, 1.0e-9) / max(initial_inflection_point_neg_mpa, 1.0e-9)))
        objective = wp_error + regularization

        if best is None or objective < best[0]:
            best = (objective, float(inflection), theta_r_val, alpha_val, n_val)

    if best is None:
        raise ValueError(
            "Could not tune a feasible inflection-point suction for the provided "
            f"theta_sat={theta_sat}, theta_fc={theta_fc}, theta_wp={theta_wp}."
        )

    _, inflection_best, theta_r_best, alpha_best, n_best = best
    return inflection_best, theta_r_best, alpha_best, n_best




def compute_van_genuchten_raw_effective_saturation(
    theta_liquid: np.ndarray | float,
    theta_sat: np.ndarray | float,
    theta_r: np.ndarray | float,
) -> np.ndarray:
    """Compute effective saturation for the van Genuchten retention curve."""
    theta_liquid = np.asarray(theta_liquid, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_r = np.asarray(theta_r, dtype=float)

    return np.clip(
        (theta_liquid - theta_r) / np.maximum(theta_sat - theta_r, 1e-9),
        0.0,
        1.0,
    )


def compute_van_genuchten_liquid_water_content(
    suction_mpa: np.ndarray | float,
    theta_sat: np.ndarray | float,
    theta_r: np.ndarray | float,
    alpha_vg_mpa_inv: np.ndarray | float,
    n_vg: np.ndarray | float,
    m_vg: np.ndarray | float,
) -> np.ndarray:
    """Evaluate liquid water content from matric suction via van Genuchten."""
    suction_mpa = np.maximum(np.asarray(suction_mpa, dtype=float), 0.0)
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_r = np.asarray(theta_r, dtype=float)
    alpha_vg_mpa_inv = np.asarray(alpha_vg_mpa_inv, dtype=float)
    n_vg = np.asarray(n_vg, dtype=float)
    m_vg = np.asarray(m_vg, dtype=float)

    effective_saturation = (
        1.0 + (np.maximum(alpha_vg_mpa_inv, 1e-9) * suction_mpa) ** n_vg
    ) ** (-m_vg)
    return theta_r + (theta_sat - theta_r) * effective_saturation


def compute_freezing_state(
    theta_total: np.ndarray | float,
    theta_sat: np.ndarray | float,
    theta_r: np.ndarray | float,
    alpha_vg_mpa_inv: np.ndarray | float,
    n_vg: np.ndarray | float,
    m_vg: np.ndarray | float,
    soil_temp_c: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Mechanistic freeze-thaw partition following Dall'Amico et al. (2011).

    Total water is mapped to an unfrozen reference suction head through the
    van Genuchten retention relation. Under freezing conditions, Clapeyron
    depression increases the liquid-water suction, and the residual between
    total and liquid water is assigned to ice.
    """
    theta_total = np.asarray(theta_total, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_r = np.asarray(theta_r, dtype=float)
    alpha_vg_mpa_inv = np.asarray(alpha_vg_mpa_inv, dtype=float)
    n_vg = np.asarray(n_vg, dtype=float)
    m_vg = np.asarray(m_vg, dtype=float)
    soil_temp_c = np.asarray(soil_temp_c, dtype=float)

    theta_total = np.clip(theta_total, 0.0, theta_sat)
    base_suction_head_m = compute_total_water_freezing_pressure_head_m(
        theta_total=theta_total,
        theta_sat=theta_sat,
        theta_r=theta_r,
        alpha_vg_mpa_inv=alpha_vg_mpa_inv,
        n_vg=n_vg,
        m_vg=m_vg,
    )
    depressed_freezing_temp_k = FREEZING_TEMPERATURE_K - (
        GRAVITY_M_S2 * FREEZING_TEMPERATURE_K / LATENT_HEAT_FUSION_J_KG
    ) * base_suction_head_m
    soil_temp_k = np.asarray(soil_temp_c, dtype=float) + FREEZING_TEMPERATURE_K

    freeze_temp_offset_c = depressed_freezing_temp_k - soil_temp_k
    freeze_smoothing_k = max(FREEZE_CURVE_SMOOTHING_C, 1e-6)
    freeze_argument = np.clip(freeze_temp_offset_c / freeze_smoothing_k, -60.0, 60.0)
    frozen_fraction = 1.0 / (1.0 + np.exp(-freeze_argument))
    freezing_increment_m = (
        LATENT_HEAT_FUSION_J_KG
        / np.maximum(GRAVITY_M_S2 * depressed_freezing_temp_k, 1e-9)
    ) * np.maximum(freeze_temp_offset_c, 0.0)
    freezing_suction_head_m = base_suction_head_m + frozen_fraction * freezing_increment_m

    suction_mpa = (
        np.maximum(freezing_suction_head_m, 0.0)
        * WATER_DENSITY_KG_M3
        * GRAVITY_M_S2
        / 1.0e6
    )
    theta_liquid = compute_van_genuchten_liquid_water_content(
        suction_mpa=suction_mpa,
        theta_sat=theta_sat,
        theta_r=theta_r,
        alpha_vg_mpa_inv=alpha_vg_mpa_inv,
        n_vg=n_vg,
        m_vg=m_vg,
    )
    theta_liquid = np.minimum(theta_liquid, theta_total)
    theta_ice = np.maximum(theta_total - theta_liquid, 0.0)

    return (
        theta_liquid,
        theta_ice,
        base_suction_head_m,
        freezing_suction_head_m,
        depressed_freezing_temp_k - FREEZING_TEMPERATURE_K,
    )


def compute_van_genuchten_pressure_head_m(
    theta_liquid: np.ndarray | float,
    theta_sat: np.ndarray | float,
    theta_r: np.ndarray | float,
    alpha_vg_mpa_inv: np.ndarray | float,
    n_vg: np.ndarray | float,
    m_vg: np.ndarray | float,
) -> np.ndarray:
    """Invert the van Genuchten relation from liquid water content to suction head."""
    theta_liquid = np.asarray(theta_liquid, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_r = np.asarray(theta_r, dtype=float)
    alpha_vg_mpa_inv = np.asarray(alpha_vg_mpa_inv, dtype=float)
    n_vg = np.asarray(n_vg, dtype=float)
    m_vg = np.asarray(m_vg, dtype=float)
    effective_saturation = compute_van_genuchten_raw_effective_saturation(
        theta_liquid=theta_liquid,
        theta_sat=theta_sat,
        theta_r=theta_r,
    )
    effective_saturation = np.maximum(effective_saturation, 1.0e-9)
    suction_mpa = np.zeros_like(effective_saturation, dtype=float)
    unsat_mask = effective_saturation < 0.999999
    suction_mpa[unsat_mask] = (
        (
            effective_saturation[unsat_mask] ** (-1.0 / np.maximum(m_vg[unsat_mask], 1e-9))
            - 1.0
        ) ** (1.0 / np.maximum(n_vg[unsat_mask], 1e-9))
    ) / np.maximum(alpha_vg_mpa_inv[unsat_mask], 1e-9)

    pressure_head_m = suction_mpa * 1.0e6 / (WATER_DENSITY_KG_M3 * GRAVITY_M_S2)
    return np.clip(pressure_head_m, 0.0, 1.0e4)


def compute_total_water_freezing_pressure_head_m(
    theta_total: np.ndarray | float,
    theta_sat: np.ndarray | float,
    theta_r: np.ndarray | float,
    alpha_vg_mpa_inv: np.ndarray | float,
    n_vg: np.ndarray | float,
    m_vg: np.ndarray | float,
) -> np.ndarray:
    """
    Total-water suction used only for freezing-point depression.

    Flow and uptake keep using the liquid-water suction relation. Freezing-point
    depression is based on the amount of pore water present as liquid plus ice.
    """
    theta_total = np.asarray(theta_total, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_r = np.asarray(theta_r, dtype=float)
    alpha_vg_mpa_inv = np.asarray(alpha_vg_mpa_inv, dtype=float)
    n_vg = np.asarray(n_vg, dtype=float)
    m_vg = np.asarray(m_vg, dtype=float)

    total_saturation = np.clip(
        (theta_total - theta_r) / np.maximum(theta_sat - theta_r, 1.0e-9),
        1.0e-9,
        1.0,
    )
    suction_mpa = np.zeros_like(total_saturation, dtype=float)
    unsat_mask = total_saturation < 0.999999
    suction_mpa[unsat_mask] = (
        (
            total_saturation[unsat_mask]
            ** (-1.0 / np.maximum(m_vg[unsat_mask], 1.0e-9))
            - 1.0
        ) ** (1.0 / np.maximum(n_vg[unsat_mask], 1.0e-9))
    ) / np.maximum(alpha_vg_mpa_inv[unsat_mask], 1.0e-9)

    pressure_head_m = suction_mpa * 1.0e6 / (WATER_DENSITY_KG_M3 * GRAVITY_M_S2)
    return np.clip(pressure_head_m, 0.0, 1.0e4)


def compute_mualem_unsat_conductivity_m_per_s(
    theta_liquid: np.ndarray | float,
    theta_sat: np.ndarray | float,
    theta_r: np.ndarray | float,
    theta_total: np.ndarray | float,
    theta_ice: np.ndarray | float,
    m_vg: np.ndarray | float,
    ksat_m_per_s: np.ndarray | float,
    tortuosity_power: float = 0.5,
) -> np.ndarray:
    """Compute unsaturated hydraulic conductivity with ice impedance."""
    theta_liquid = np.asarray(theta_liquid, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_r = np.asarray(theta_r, dtype=float)
    theta_total = np.asarray(theta_total, dtype=float)
    theta_ice = np.asarray(theta_ice, dtype=float)
    m_vg = np.asarray(m_vg, dtype=float)
    ksat_m_per_s = np.asarray(ksat_m_per_s, dtype=float)
    effective_saturation = compute_van_genuchten_raw_effective_saturation(
        theta_liquid=theta_liquid,
        theta_sat=theta_sat,
        theta_r=theta_r,
    )
    effective_saturation = np.clip(effective_saturation, 0.0, 1.0)
    m_safe = np.maximum(m_vg, 1e-9)
    relative_conductivity = (
        effective_saturation ** tortuosity_power
    ) * (1.0 - (1.0 - effective_saturation ** (1.0 / m_safe)) ** m_safe) ** 2
    conductivity = (
        ksat_m_per_s
        * relative_conductivity
    )
    liquid_plus_ice = np.maximum(theta_total, 1e-9)
    ice_fraction = np.clip(theta_ice / liquid_plus_ice, 0.0, 1.0)
    ice_impedance = np.exp(-ICE_IMPEDANCE_SHAPE_FACTOR * ice_fraction)
    conductivity = conductivity * ice_impedance
    conductivity = np.where(effective_saturation >= 0.999999, ksat_m_per_s, conductivity)
    return np.clip(conductivity, 0.0, ksat_m_per_s)


def compute_root_water_uptake(
    pet_mm: float,
    growth_stage_progress: float,
    root_weights: np.ndarray,
    theta_liquid: np.ndarray,
    theta_total: np.ndarray | None,
    theta_ice: np.ndarray | None,
    theta_sat: np.ndarray,
    theta_fc: np.ndarray,
    theta_stress_floor: np.ndarray,
    theta_wp: np.ndarray,
    thickness_mm: np.ndarray,
    theta_r: np.ndarray,
    alpha_vg_mpa_inv: np.ndarray,
    n_vg: np.ndarray,
    m_vg: np.ndarray,
    ksat_m_per_s: np.ndarray,
    dt_seconds: float,
) -> tuple[np.ndarray, float]:
    """
    Layer-wise root uptake using root distribution and soil-water stress.

    Potential transpiration is allocated across the rooting zone according to
    root weights and a Feddes-style water-stress factor based on plant-
    extractable liquid water. Uptake is limited by removable liquid water
    above wilting point and redistributed iteratively to layers that still
    have available supply.
    """
    theta_liquid = np.asarray(theta_liquid, dtype=float)
    theta_total = (
        np.asarray(theta_total, dtype=float)
        if theta_total is not None
        else np.maximum(theta_liquid, 0.0)
    )
    theta_ice = (
        np.asarray(theta_ice, dtype=float)
        if theta_ice is not None
        else np.zeros_like(theta_liquid)
    )
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_fc = np.asarray(theta_fc, dtype=float)
    theta_stress_floor = np.asarray(theta_stress_floor, dtype=float)
    theta_wp = np.asarray(theta_wp, dtype=float)
    root_weights = np.asarray(root_weights, dtype=float)
    thickness_mm = np.asarray(thickness_mm, dtype=float)
    theta_r = np.asarray(theta_r, dtype=float)
    alpha_vg_mpa_inv = np.asarray(alpha_vg_mpa_inv, dtype=float)
    n_vg = np.asarray(n_vg, dtype=float)
    m_vg = np.asarray(m_vg, dtype=float)
    ksat_m_per_s = np.asarray(ksat_m_per_s, dtype=float)
    growth_stage_progress = clamp(float(growth_stage_progress), 0.0, 1.0)

    relative_extractable_water = np.clip(
        (theta_liquid - theta_wp)
        / np.maximum(theta_fc - theta_wp, 1.0e-9),
        0.0,
        1.0,
    )
    moisture_stress = relative_extractable_water
    moisture_stress = np.where(theta_liquid <= theta_wp, 0.0, moisture_stress)
    frozen_stress = np.clip(1.0 - (theta_ice / np.maximum(theta_total, 1.0e-9)), 0.0, 1.0)
    moisture_stress = moisture_stress * frozen_stress
    removable_liquid_mm = np.maximum((theta_liquid - theta_wp) * thickness_mm, 0.0)
    activity_scale = 0.25 + 0.75 * np.sin(0.5 * np.pi * growth_stage_progress)
    demand_weights = np.maximum(root_weights, 0.0) * moisture_stress * activity_scale
    layer_capacity_mm = removable_liquid_mm.copy()

    total_capacity_mm = float(np.sum(layer_capacity_mm))
    if pet_mm <= 0.0 or total_capacity_mm <= 1e-12:
        return np.zeros_like(layer_capacity_mm), 0.0

    actual_uptake_mm = min(float(pet_mm), total_capacity_mm)
    uptake_by_layer_mm = np.zeros_like(layer_capacity_mm)
    remaining_demand_mm = actual_uptake_mm
    remaining_capacity_mm = layer_capacity_mm.copy()
    remaining_weights = demand_weights.copy()

    for _ in range(len(layer_capacity_mm)):
        weight_sum = float(np.sum(remaining_weights))
        if remaining_demand_mm <= 1e-12 or weight_sum <= 1e-12:
            break

        proposed_uptake_mm = remaining_demand_mm * (remaining_weights / weight_sum)
        accepted_uptake_mm = np.minimum(proposed_uptake_mm, remaining_capacity_mm)
        uptake_by_layer_mm += accepted_uptake_mm
        remaining_demand_mm -= float(np.sum(accepted_uptake_mm))
        remaining_capacity_mm = np.maximum(remaining_capacity_mm - accepted_uptake_mm, 0.0)
        remaining_weights = np.where(remaining_capacity_mm > 1.0e-12, remaining_weights, 0.0)

    return uptake_by_layer_mm, actual_uptake_mm


def compute_substep_coupled_snow_balance(
    precip_mm: float,
    swe_mm: float,
    snow_depth_mm: float,
    snow_layer_swe_mm: np.ndarray | None,
    snow_layer_age_days: np.ndarray | None,
    air_temp_c: float,
    top_soil_temp_c: float,
    top_soil_thermal_conductivity_w_m_k: float,
    top_soil_half_thickness_m: float,
    dt_seconds: float,
    melt_factor_mm_per_degC_day: float,
) -> tuple[float, float, float, float, float, float, float, np.ndarray, np.ndarray]:
    """
    State-dependent substep snow partition and melt with an explicit five-layer pack.

    Fresh snowfall is added as a new top layer and aging snow becomes denser and
    more conductive than newly fallen snow. The effective conductivity passed to
    the surface-energy partition is the harmonic aggregate of the layer stack.
    """
    precip_mm = max(float(precip_mm), 0.0)
    swe_mm = max(float(swe_mm), 0.0)
    snow_depth_mm = max(float(snow_depth_mm), 0.0)
    air_temp_c = float(air_temp_c)
    top_soil_temp_c = float(top_soil_temp_c)
    top_soil_thermal_conductivity_w_m_k = max(float(top_soil_thermal_conductivity_w_m_k), 0.05)
    top_soil_half_thickness_m = max(float(top_soil_half_thickness_m), 1.0e-6)
    dt_seconds = max(float(dt_seconds), 1.0)
    melt_factor_mm_per_degC_day = max(float(melt_factor_mm_per_degC_day), 0.0)

    surface_temp_c = compute_snow_surface_temperature_c(
        air_temp_c=air_temp_c,
        top_soil_temp_c=top_soil_temp_c,
        swe_mm=swe_mm,
        snow_depth_mm=snow_depth_mm,
        top_soil_thermal_conductivity_w_m_k=top_soil_thermal_conductivity_w_m_k,
        top_soil_half_thickness_m=top_soil_half_thickness_m,
        snow_layer_swe_mm=snow_layer_swe_mm,
        snow_layer_age_days=snow_layer_age_days,
    )
    phase_temp_c = 0.65 * air_temp_c + 0.35 * surface_temp_c
    rain_fraction = clamp(0.5 + 0.5 * phase_temp_c, 0.0, 1.0)

    rain_mm = precip_mm * rain_fraction
    snowfall_water_mm = precip_mm - rain_mm

    thaw_degree_c = max(0.5 * (max(air_temp_c, 0.0) + max(surface_temp_c, 0.0)), 0.0)
    melt_potential_mm = melt_factor_mm_per_degC_day * thaw_degree_c * (dt_seconds / 86400.0)
    available_swe_mm = swe_mm + snowfall_water_mm
    snowmelt_mm = min(available_swe_mm, melt_potential_mm)

    updated_layer_swe_mm, updated_layer_age_days = evolve_snowpack_layers(
        layer_swe_mm=snow_layer_swe_mm,
        layer_age_days=snow_layer_age_days,
        snowfall_water_mm=snowfall_water_mm,
        snowmelt_mm=snowmelt_mm,
        dt_seconds=dt_seconds,
    )
    updated_swe_mm, updated_snow_depth_mm, updated_snow_density_kg_m3, _ = compute_snowpack_bulk_properties(
        layer_swe_mm=updated_layer_swe_mm,
        layer_age_days=updated_layer_age_days,
        air_temp_c=air_temp_c,
        surface_temp_c=surface_temp_c,
    )

    return (
        rain_mm,
        snowfall_water_mm,
        snowmelt_mm,
        thaw_degree_c,
        updated_swe_mm,
        updated_snow_depth_mm,
        updated_snow_density_kg_m3,
        updated_layer_swe_mm,
        updated_layer_age_days,
    )


def apply_simple_richards_flow(
    storage_mm: np.ndarray,
    thickness_mm: np.ndarray,
    mid_cm: np.ndarray,
    soil_temp_c: np.ndarray,
    theta_sat: np.ndarray,
    theta_fc: np.ndarray,
    theta_r: np.ndarray,
    alpha_vg_mpa_inv: np.ndarray,
    n_vg: np.ndarray,
    m_vg: np.ndarray,
    ksat_m_per_s: np.ndarray,
    deep_drainage_impedence_factor: float,
    dt_seconds: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Simplified 1-D Richards-type redistribution using layer-based conductivities.

    Interface conductivity uses the harmonic mean of adjacent saturated
    conductivities combined with the arithmetic mean of adjacent relative
    conductivities, which softens dry-layer choking compared with applying
    a harmonic mean directly to the unsaturated conductivities.
    """
    storage_mm = np.asarray(storage_mm, dtype=float).copy()
    thickness_mm = np.asarray(thickness_mm, dtype=float)
    mid_cm = np.asarray(mid_cm, dtype=float)
    soil_temp_c = np.asarray(soil_temp_c, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_fc = np.asarray(theta_fc, dtype=float)
    theta_r = np.asarray(theta_r, dtype=float)
    alpha_vg_mpa_inv = np.asarray(alpha_vg_mpa_inv, dtype=float)
    n_vg = np.asarray(n_vg, dtype=float)
    m_vg = np.asarray(m_vg, dtype=float)
    ksat_m_per_s = np.asarray(ksat_m_per_s, dtype=float)
    deep_drainage_impedence_factor = float(np.clip(deep_drainage_impedence_factor, 0.0, 1.0))
    dt_seconds = max(float(dt_seconds), 1.0)

    n_layers = len(storage_mm)
    interlayer_flux_mm = np.zeros(max(n_layers - 1, 0), dtype=float)
    for j in range(n_layers - 1):
        theta_total = np.clip(storage_mm / thickness_mm, theta_r, theta_sat)
        theta_liquid, theta_ice, _, _, _ = compute_freezing_state(
            theta_total=theta_total,
            theta_sat=theta_sat,
            theta_r=theta_r,
            alpha_vg_mpa_inv=alpha_vg_mpa_inv,
            n_vg=n_vg,
            m_vg=m_vg,
            soil_temp_c=soil_temp_c,
        )

        suction_upper_m = float(
            compute_van_genuchten_pressure_head_m(
                theta_liquid=theta_liquid[j],
                theta_sat=theta_sat[j],
                theta_r=theta_r[j],
                alpha_vg_mpa_inv=alpha_vg_mpa_inv[j],
                n_vg=n_vg[j],
                m_vg=m_vg[j],
            )
        )
        suction_lower_m = float(
            compute_van_genuchten_pressure_head_m(
                theta_liquid=theta_liquid[j + 1],
                theta_sat=theta_sat[j + 1],
                theta_r=theta_r[j + 1],
                alpha_vg_mpa_inv=alpha_vg_mpa_inv[j + 1],
                n_vg=n_vg[j + 1],
                m_vg=m_vg[j + 1],
            )
        )

        conductivity_upper = float(
            compute_mualem_unsat_conductivity_m_per_s(
                theta_liquid=theta_liquid[j],
                theta_sat=theta_sat[j],
                theta_r=theta_r[j],
                theta_total=theta_total[j],
                theta_ice=theta_ice[j],
                m_vg=m_vg[j],
                ksat_m_per_s=ksat_m_per_s[j],
            )
        )
        conductivity_lower = float(
            compute_mualem_unsat_conductivity_m_per_s(
                theta_liquid=theta_liquid[j + 1],
                theta_sat=theta_sat[j + 1],
                theta_r=theta_r[j + 1],
                theta_total=theta_total[j + 1],
                theta_ice=theta_ice[j + 1],
                m_vg=m_vg[j + 1],
                ksat_m_per_s=ksat_m_per_s[j + 1],
            )
        )
        relative_conductivity_upper = conductivity_upper / max(ksat_m_per_s[j], 1.0e-18)
        relative_conductivity_lower = conductivity_lower / max(ksat_m_per_s[j + 1], 1.0e-18)
        relative_conductivity_interface = 0.5 * (
            relative_conductivity_upper + relative_conductivity_lower
        )

        upper_half_m = max(0.5 * thickness_mm[j] / 1000.0, 1.0e-6)
        lower_half_m = max(0.5 * thickness_mm[j + 1] / 1000.0, 1.0e-6)
        hydraulic_resistance = (
            upper_half_m / max(ksat_m_per_s[j], 1.0e-18)
        ) + (
            lower_half_m / max(ksat_m_per_s[j + 1], 1.0e-18)
        )
        interface_ksat = (upper_half_m + lower_half_m) / max(hydraulic_resistance, 1.0e-12)
        interface_conductivity = interface_ksat * np.clip(relative_conductivity_interface, 0.0, 1.0)

        center_spacing_m = max((mid_cm[j + 1] - mid_cm[j]) / 100.0, 1.0e-6)
        total_gradient = 1.0 + ((suction_lower_m - suction_upper_m) / center_spacing_m)
        total_gradient = float(np.clip(total_gradient, -2.5, 2.5))

        potential_flux_mm = interface_conductivity * total_gradient * dt_seconds * 1000.0

        if potential_flux_mm >= 0.0:
            donor_liquid_mm = max(
                theta_liquid[j] * thickness_mm[j] - theta_r[j] * thickness_mm[j],
                0.0,
            )
            receiver_space_mm = max(
                theta_sat[j + 1] * thickness_mm[j + 1] - storage_mm[j + 1],
                0.0,
            )
            actual_flux_mm = min(
                potential_flux_mm,
                donor_liquid_mm,
                receiver_space_mm,
                0.35 * max(donor_liquid_mm, 0.0),
                0.35 * max(receiver_space_mm, 0.0),
            )
            storage_mm[j] -= actual_flux_mm
            storage_mm[j + 1] += actual_flux_mm
            interlayer_flux_mm[j] = actual_flux_mm
        else:
            donor_liquid_mm = max(
                theta_liquid[j + 1] * thickness_mm[j + 1] - theta_r[j + 1] * thickness_mm[j + 1],
                0.0,
            )
            receiver_space_mm = max(
                theta_sat[j] * thickness_mm[j] - storage_mm[j],
                0.0,
            )
            actual_flux_mm = min(
                -potential_flux_mm,
                donor_liquid_mm,
                receiver_space_mm,
                0.35 * max(donor_liquid_mm, 0.0),
                0.35 * max(receiver_space_mm, 0.0),
            )
            storage_mm[j] += actual_flux_mm
            storage_mm[j + 1] -= actual_flux_mm
            interlayer_flux_mm[j] = -actual_flux_mm
    theta_total = storage_mm / thickness_mm
    theta_liquid, theta_ice, _, _, _ = compute_freezing_state(
        theta_total=theta_total,
        theta_sat=theta_sat,
        theta_r=theta_r,
        alpha_vg_mpa_inv=alpha_vg_mpa_inv,
        n_vg=n_vg,
        m_vg=m_vg,
        soil_temp_c=soil_temp_c,
    )
    bottom_conductivity = float(
        compute_mualem_unsat_conductivity_m_per_s(
            theta_liquid=theta_liquid[-1],
            theta_sat=theta_sat[-1],
            theta_r=theta_r[-1],
            theta_total=theta_total[-1],
            theta_ice=theta_ice[-1],
            m_vg=m_vg[-1],
            ksat_m_per_s=ksat_m_per_s[-1],
        )
    )
    donor_liquid_bottom_mm = max(
        theta_liquid[-1] * thickness_mm[-1] - theta_fc[-1] * thickness_mm[-1],
        0.0,
    )
    deep_drainage_mm = min(
        bottom_conductivity * dt_seconds * 1000.0,
        donor_liquid_bottom_mm,
        deep_drainage_impedence_factor * max(donor_liquid_bottom_mm, 0.0),
    )
    storage_mm[-1] -= deep_drainage_mm

    for j in range(n_layers):
        saturation_excess_mm = max(storage_mm[j] - theta_sat[j] * thickness_mm[j], 0.0)
        if saturation_excess_mm <= 0.0:
            continue

        storage_mm[j] -= saturation_excess_mm
        if j < n_layers - 1:
            storage_mm[j + 1] += saturation_excess_mm
            interlayer_flux_mm[j] += saturation_excess_mm
        else:
            deep_drainage_mm += saturation_excess_mm

    return storage_mm, interlayer_flux_mm, float(deep_drainage_mm)


def partition_precipitation_and_snowmelt(
    precip_mm: float,
    swe_mm: float,
    daily_max_air_temp_c: float,
    daily_min_air_temp_c: float,
    daily_mean_air_temp_c: float,
    melt_factor_mm_per_degC_day: float,
) -> tuple[float, float, float, float]:
    """Partition precipitation into rain and snow and estimate degree-day melt."""
    precip_mm = max(float(precip_mm), 0.0)
    swe_mm = max(float(swe_mm), 0.0)
    daily_max_air_temp_c = float(daily_max_air_temp_c)
    daily_min_air_temp_c = float(daily_min_air_temp_c)
    daily_mean_air_temp_c = float(daily_mean_air_temp_c)

    if daily_max_air_temp_c <= 0.0:
        rain_fraction = 0.0
    elif daily_min_air_temp_c >= 0.0:
        rain_fraction = 1.0
    else:
        rain_fraction = clamp(
            daily_max_air_temp_c / max(daily_max_air_temp_c - daily_min_air_temp_c, 1e-6),
            0.0,
            1.0,
        )

    rain_mm = precip_mm * rain_fraction
    snowfall_water_mm = precip_mm - rain_mm

    thaw_degree_c = max(0.5 * (max(daily_max_air_temp_c, 0.0) + max(daily_mean_air_temp_c, 0.0)), 0.0)
    snowmelt_mm = min(swe_mm + snowfall_water_mm, melt_factor_mm_per_degC_day * thaw_degree_c)

    return rain_mm, snowfall_water_mm, snowmelt_mm, thaw_degree_c
