from __future__ import annotations

"""Thermal and enthalpy utilities for the soil temperature solver.

References:
- Johansen-style bulk thermal conductivity concepts for porous media.
- Dall'Amico et al. (2011) style freeze-thaw enthalpy partitioning via the
  freezing curve implemented in :mod:`ab_ag_soil_hydro_temp.water`.
"""

import numpy as np

from .common import compute_snow_surface_temperature_c
from .water import LATENT_HEAT_FUSION_J_KG, compute_freezing_state

WATER_VOLUMETRIC_HEAT_CAPACITY_J_M3_K = 4.18e6
ICE_VOLUMETRIC_HEAT_CAPACITY_J_M3_K = 1.93e6
SOLID_VOLUMETRIC_HEAT_CAPACITY_J_M3_K = 2.0e6
AIR_VOLUMETRIC_HEAT_CAPACITY_J_M3_K = 1.2e3
LATENT_HEAT_VOLUMETRIC_J_M3 = 1000.0 * LATENT_HEAT_FUSION_J_KG
ENTHALPY_INVERSION_TOL_C = 5.0e-4
ENTHALPY_INVERSION_MAX_ITERATIONS = 12
FREEZE_ZONE_ENTHALPY_RESIDUAL_TOL_J_M3 = 10.0
BACKGROUND_ENTHALPY_RESIDUAL_TOL_J_M3 = 50.0
MINERAL_THERMAL_CONDUCTIVITY_W_M_K = 2.5
WATER_THERMAL_CONDUCTIVITY_W_M_K = 0.57
ICE_THERMAL_CONDUCTIVITY_W_M_K = 2.24
AIR_THERMAL_CONDUCTIVITY_W_M_K = 0.025
PARTICLE_DENSITY_KG_M3 = 2650.0


def compute_volumetric_sensible_heat_capacity_j_m3_k(
    theta_total: np.ndarray | float,
    theta_liquid: np.ndarray | float,
    theta_ice: np.ndarray | float,
    porosity: np.ndarray | float,
) -> np.ndarray:
    """Compute bulk sensible heat capacity for the soil matrix and pore space."""
    theta_total = np.asarray(theta_total, dtype=float)
    theta_liquid = np.asarray(theta_liquid, dtype=float)
    theta_ice = np.asarray(theta_ice, dtype=float)
    porosity = np.clip(np.asarray(porosity, dtype=float), 0.05, 0.80)

    solids_fraction = np.clip(1.0 - porosity, 1.0e-6, 1.0)
    air_fraction = np.clip(porosity - theta_total, 0.0, porosity)

    return (
        solids_fraction * SOLID_VOLUMETRIC_HEAT_CAPACITY_J_M3_K
        + theta_liquid * WATER_VOLUMETRIC_HEAT_CAPACITY_J_M3_K
        + theta_ice * ICE_VOLUMETRIC_HEAT_CAPACITY_J_M3_K
        + air_fraction * AIR_VOLUMETRIC_HEAT_CAPACITY_J_M3_K
    )


def compute_soil_thermal_diffusivity_m2_s(
    theta_total: np.ndarray,
    theta_liquid: np.ndarray,
    theta_ice: np.ndarray,
    porosity: np.ndarray,
) -> np.ndarray:
    """Estimate bulk soil thermal diffusivity from composition and saturation."""
    theta_total = np.asarray(theta_total, dtype=float)
    theta_liquid = np.asarray(theta_liquid, dtype=float)
    theta_ice = np.asarray(theta_ice, dtype=float)
    porosity = np.clip(np.asarray(porosity, dtype=float), 0.05, 0.80)
    solids_fraction = np.clip(1.0 - porosity, 1.0e-6, 1.0)
    saturation = np.clip(theta_total / np.maximum(porosity, 1.0e-9), 0.0, 1.0)
    frozen_fraction = np.clip(theta_ice / np.maximum(theta_total, 1.0e-9), 0.0, 1.0)
    bulk_density_kg_m3 = solids_fraction * PARTICLE_DENSITY_KG_M3
    dry_thermal_conductivity = (
        0.135 * bulk_density_kg_m3 + 64.7
    ) / np.maximum(PARTICLE_DENSITY_KG_M3 - (0.947 * bulk_density_kg_m3), 1.0)
    dry_thermal_conductivity = np.clip(dry_thermal_conductivity, 0.05, 1.5)

    pore_fluid_thermal_conductivity = (
        WATER_THERMAL_CONDUCTIVITY_W_M_K ** np.clip(porosity - theta_ice, 0.0, porosity)
    ) * (
        ICE_THERMAL_CONDUCTIVITY_W_M_K ** np.clip(theta_ice, 0.0, porosity)
    ) * (
        AIR_THERMAL_CONDUCTIVITY_W_M_K ** np.clip(porosity - theta_total, 0.0, porosity)
    )
    saturated_thermal_conductivity = (
        MINERAL_THERMAL_CONDUCTIVITY_W_M_K ** solids_fraction
    ) * pore_fluid_thermal_conductivity

    ker_unfrozen = np.where(
        saturation > 0.1,
        np.log10(np.maximum(saturation, 1.0e-6)) + 1.0,
        0.0,
    )
    ker_frozen = saturation
    ker = np.clip(
        (1.0 - frozen_fraction) * ker_unfrozen + frozen_fraction * ker_frozen,
        0.0,
        1.0,
    )
    thermal_conductivity = dry_thermal_conductivity + (
        ker * (saturated_thermal_conductivity - dry_thermal_conductivity)
    )
    volumetric_heat_capacity = compute_volumetric_sensible_heat_capacity_j_m3_k(
        theta_total=theta_total,
        theta_liquid=theta_liquid,
        theta_ice=theta_ice,
        porosity=porosity,
    )
    thermal_diffusivity = thermal_conductivity / np.maximum(volumetric_heat_capacity, 1.0)
    return np.clip(thermal_diffusivity, 0.05e-6, 1.50e-6)


def compute_volumetric_enthalpy_j_m3(
    soil_temp_c: np.ndarray | float,
    theta_total: np.ndarray | float,
    theta_liquid: np.ndarray | float,
    theta_ice: np.ndarray | float,
    porosity: np.ndarray | float,
) -> np.ndarray:
    """Convert temperature and phase partition into volumetric enthalpy."""
    soil_temp_c = np.asarray(soil_temp_c, dtype=float)
    theta_ice = np.asarray(theta_ice, dtype=float)
    sensible_heat_capacity = compute_volumetric_sensible_heat_capacity_j_m3_k(
        theta_total=theta_total,
        theta_liquid=theta_liquid,
        theta_ice=theta_ice,
        porosity=porosity,
    )
    return sensible_heat_capacity * soil_temp_c - (LATENT_HEAT_VOLUMETRIC_J_M3 * theta_ice)


def invert_enthalpy_to_temperature(
    enthalpy_j_m3: np.ndarray | float,
    theta_total: np.ndarray | float,
    theta_sat: np.ndarray | float,
    theta_residual: np.ndarray | float,
    alpha_vg_mpa_inv: np.ndarray | float,
    n_vg: np.ndarray | float,
    m_vg: np.ndarray | float,
    porosity: np.ndarray | float,
    initial_temp_c: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Invert volumetric enthalpy to temperature and phase partition."""
    enthalpy_j_m3 = np.asarray(enthalpy_j_m3, dtype=float)
    theta_total = np.asarray(theta_total, dtype=float)
    theta_sat = np.asarray(theta_sat, dtype=float)
    theta_residual = np.asarray(theta_residual, dtype=float)
    alpha_vg_mpa_inv = np.asarray(alpha_vg_mpa_inv, dtype=float)
    n_vg = np.asarray(n_vg, dtype=float)
    m_vg = np.asarray(m_vg, dtype=float)
    porosity = np.asarray(porosity, dtype=float)
    trial_temp_c = np.asarray(initial_temp_c, dtype=float).copy()

    enthalpy_safe = np.asarray(enthalpy_j_m3, dtype=float)
    for _ in range(ENTHALPY_INVERSION_MAX_ITERATIONS):
        theta_liquid, theta_ice, _, _, depressed_freezing_temp_c = compute_freezing_state(
            theta_total=theta_total,
            theta_sat=theta_sat,
            theta_r=theta_residual,
            alpha_vg_mpa_inv=alpha_vg_mpa_inv,
            n_vg=n_vg,
            m_vg=m_vg,
            soil_temp_c=trial_temp_c,
        )
        sensible_heat_capacity = compute_volumetric_sensible_heat_capacity_j_m3_k(
            theta_total=theta_total,
            theta_liquid=theta_liquid,
            theta_ice=theta_ice,
            porosity=porosity,
        )
        trial_enthalpy = sensible_heat_capacity * trial_temp_c - (LATENT_HEAT_VOLUMETRIC_J_M3 * theta_ice)
        residual = trial_enthalpy - enthalpy_safe
        residual_tolerance = np.where(
            np.abs(trial_temp_c - depressed_freezing_temp_c) <= 0.75,
            FREEZE_ZONE_ENTHALPY_RESIDUAL_TOL_J_M3,
            BACKGROUND_ENTHALPY_RESIDUAL_TOL_J_M3,
        )
        if bool(np.all(np.abs(residual) <= residual_tolerance)):
            return trial_temp_c, theta_liquid, theta_ice, depressed_freezing_temp_c

        delta_temp_c = 0.02
        theta_liquid_plus, theta_ice_plus, _, _, _ = compute_freezing_state(
            theta_total=theta_total,
            theta_sat=theta_sat,
            theta_r=theta_residual,
            alpha_vg_mpa_inv=alpha_vg_mpa_inv,
            n_vg=n_vg,
            m_vg=m_vg,
            soil_temp_c=trial_temp_c + delta_temp_c,
        )
        sensible_heat_capacity_plus = compute_volumetric_sensible_heat_capacity_j_m3_k(
            theta_total=theta_total,
            theta_liquid=theta_liquid_plus,
            theta_ice=theta_ice_plus,
            porosity=porosity,
        )
        enthalpy_plus = (
            sensible_heat_capacity_plus * (trial_temp_c + delta_temp_c)
            - (LATENT_HEAT_VOLUMETRIC_J_M3 * theta_ice_plus)
        )
        d_enthalpy_d_temp = (enthalpy_plus - trial_enthalpy) / delta_temp_c
        d_enthalpy_d_temp = np.where(
            np.abs(d_enthalpy_d_temp) < 1.0e3,
            np.sign(d_enthalpy_d_temp) * 1.0e3 + (d_enthalpy_d_temp == 0.0) * 1.0e3,
            d_enthalpy_d_temp,
        )
        temp_step = -residual / d_enthalpy_d_temp
        freeze_zone_mask = np.abs(trial_temp_c - depressed_freezing_temp_c) <= 0.5
        temp_step = np.where(freeze_zone_mask, np.clip(temp_step, -0.25, 0.25), np.clip(temp_step, -2.0, 2.0))
        trial_temp_c = trial_temp_c + temp_step
        if float(np.max(np.abs(temp_step))) <= ENTHALPY_INVERSION_TOL_C:
            break

    theta_liquid, theta_ice, _, _, depressed_freezing_temp_c = compute_freezing_state(
        theta_total=theta_total,
        theta_sat=theta_sat,
        theta_r=theta_residual,
        alpha_vg_mpa_inv=alpha_vg_mpa_inv,
        n_vg=n_vg,
        m_vg=m_vg,
        soil_temp_c=trial_temp_c,
    )
    return trial_temp_c, theta_liquid, theta_ice, depressed_freezing_temp_c


def solve_enthalpy_diffusion_step(
    previous_enthalpy_j_m3: np.ndarray,
    previous_soil_temp_c: np.ndarray,
    theta_total: np.ndarray,
    theta_liquid: np.ndarray,
    theta_ice: np.ndarray,
    porosity: np.ndarray,
    theta_residual: np.ndarray,
    alpha_vg_mpa_inv: np.ndarray,
    n_vg: np.ndarray,
    m_vg: np.ndarray,
    thickness_mm: np.ndarray,
    swe_mm: float,
    snow_depth_mm: float,
    snow_layer_swe_mm: np.ndarray | None,
    snow_layer_age_days: np.ndarray | None,
    air_temp_c: float,
    deep_boundary_temp_c: float,
    lower_boundary_distance_m: float,
    dt_seconds: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Advance soil heat diffusion over one timestep."""
    previous_enthalpy_j_m3 = np.asarray(previous_enthalpy_j_m3, dtype=float)
    previous_soil_temp_c = np.asarray(previous_soil_temp_c, dtype=float)
    thickness_m = np.asarray(thickness_mm, dtype=float) / 1000.0
    theta_total = np.asarray(theta_total, dtype=float)
    theta_liquid = np.asarray(theta_liquid, dtype=float)
    theta_ice = np.asarray(theta_ice, dtype=float)
    porosity = np.asarray(porosity, dtype=float)
    theta_residual = np.asarray(theta_residual, dtype=float)
    alpha_vg_mpa_inv = np.asarray(alpha_vg_mpa_inv, dtype=float)
    n_vg = np.asarray(n_vg, dtype=float)
    m_vg = np.asarray(m_vg, dtype=float)
    lower_boundary_distance_m = max(float(lower_boundary_distance_m), 1.0e-6)

    volumetric_heat_capacity = compute_volumetric_sensible_heat_capacity_j_m3_k(
        theta_total=theta_total,
        theta_liquid=theta_liquid,
        theta_ice=theta_ice,
        porosity=porosity,
    )
    thermal_diffusivity = compute_soil_thermal_diffusivity_m2_s(
        theta_total=theta_total,
        theta_liquid=theta_liquid,
        theta_ice=theta_ice,
        porosity=porosity,
    )
    thermal_conductivity = thermal_diffusivity * volumetric_heat_capacity

    surface_temp_c = compute_snow_surface_temperature_c(
        air_temp_c=float(air_temp_c),
        top_soil_temp_c=float(previous_soil_temp_c[0]),
        swe_mm=float(swe_mm),
        snow_depth_mm=float(snow_depth_mm),
        top_soil_thermal_conductivity_w_m_k=float(max(thermal_conductivity[0], 0.05)),
        top_soil_half_thickness_m=float(max(0.5 * thickness_m[0], 1.0e-6)),
        snow_layer_swe_mm=snow_layer_swe_mm,
        snow_layer_age_days=snow_layer_age_days,
    )

    n_layers = len(previous_soil_temp_c)
    matrix = np.zeros((n_layers, n_layers), dtype=float)
    rhs = (
        volumetric_heat_capacity
        * previous_soil_temp_c
        * thickness_m
        / max(dt_seconds, 1.0e-9)
    )

    for i in range(n_layers):
        storage_coeff = (
            volumetric_heat_capacity[i]
            * thickness_m[i]
            / max(dt_seconds, 1.0e-9)
        )
        matrix[i, i] += storage_coeff

        if i == 0:
            upper_distance_m = max(0.5 * thickness_m[i], 1.0e-6)
            upper_conductance = thermal_conductivity[i] / upper_distance_m
            matrix[i, i] += upper_conductance
            rhs[i] += upper_conductance * surface_temp_c

        if i < n_layers - 1:
            center_distance_m = max(0.5 * (thickness_m[i] + thickness_m[i + 1]), 1.0e-6)
            interface_conductivity = 2.0 / max(
                (1.0 / max(thermal_conductivity[i], 1.0e-9))
                + (1.0 / max(thermal_conductivity[i + 1], 1.0e-9)),
                1.0e-9,
            )
            conductance = interface_conductivity / center_distance_m
            matrix[i, i] += conductance
            matrix[i + 1, i + 1] += conductance
            matrix[i, i + 1] -= conductance
            matrix[i + 1, i] -= conductance

        if i == n_layers - 1:
            lower_conductance = thermal_conductivity[i] / lower_boundary_distance_m
            matrix[i, i] += lower_conductance
            rhs[i] += lower_conductance * float(deep_boundary_temp_c)

    updated_soil_temp_c = np.linalg.solve(matrix, rhs)
    updated_theta_liquid, updated_theta_ice, _, _, _ = compute_freezing_state(
        theta_total=theta_total,
        theta_sat=porosity,
        theta_r=theta_residual,
        alpha_vg_mpa_inv=alpha_vg_mpa_inv,
        n_vg=n_vg,
        m_vg=m_vg,
        soil_temp_c=updated_soil_temp_c,
    )
    updated_enthalpy_j_m3 = compute_volumetric_enthalpy_j_m3(
        soil_temp_c=updated_soil_temp_c,
        theta_total=theta_total,
        theta_liquid=updated_theta_liquid,
        theta_ice=updated_theta_ice,
        porosity=porosity,
    )
    return updated_enthalpy_j_m3, thermal_diffusivity


def compute_daily_air_temperature_metrics(
    daily_max_air_temp_c: float,
    daily_min_air_temp_c: float,
) -> tuple[float, float, float]:
    """Derive daily mean temperature, range, and a cold-condition thaw bias."""
    daily_max_air_temp_c = float(daily_max_air_temp_c)
    daily_min_air_temp_c = float(daily_min_air_temp_c)
    if daily_max_air_temp_c < daily_min_air_temp_c:
        raise ValueError(
            "daily_max_air_temp must be >= daily_min_air_temp for each weather row."
        )

    daily_mean_air_temp_c = 0.5 * (daily_max_air_temp_c + daily_min_air_temp_c)
    daily_temp_range_c = daily_max_air_temp_c - daily_min_air_temp_c
    nocturnal_cooling_c = 0.18 * daily_temp_range_c
    cold_condition_factor = np.clip((2.0 - daily_mean_air_temp_c) / 10.0, 0.0, 1.0)
    thaw_bias_c = -float(np.clip(nocturnal_cooling_c * cold_condition_factor, 0.0, 3.0))

    return daily_mean_air_temp_c, daily_temp_range_c, thaw_bias_c
