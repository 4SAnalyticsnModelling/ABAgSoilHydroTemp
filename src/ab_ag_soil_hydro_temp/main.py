from __future__ import annotations

"""Top-level model orchestration for daily Alberta site simulations.

This module coordinates:
- annual lower-boundary temperature forcing from a damped harmonic fit
- daily snow, water-flow, evapotranspiration, and thermal coupling
- assembly of publication-ready output tables
"""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .common import (
    HYDROLOGY_SUBSTEPS_PER_DAY,
    depth_to_safe_name,
    initialize_snowpack_layers,
    parse_weather_dates,
)
from .crop import (
    compute_blaney_criddle_p_factor,
    compute_daylength_hours,
    compute_season_progress,
    build_crop_lookup,
    compute_root_weights,
    compute_rooting_depth_cm,
    compute_stage_crop_coefficient,
)
from .heat import (
    compute_volumetric_enthalpy_j_m3,
    compute_soil_thermal_diffusivity_m2_s,
    compute_volumetric_sensible_heat_capacity_j_m3_k,
    compute_daily_air_temperature_metrics,
    invert_enthalpy_to_temperature,
    solve_enthalpy_diffusion_step,
)
from .io import (
    build_boundary_lookup,
    build_soil_lookup,
    load_output_layer_depths,
    read_delimited_text_table,
    reorder_output_columns,
)
from .water import (
    apply_simple_richards_flow,
    compute_freezing_state,
    compute_root_water_uptake,
    compute_substep_coupled_snow_balance,
)

MAX_COUPLED_NONLINEAR_ITERATIONS = 3
PICARD_ITERATIONS_BEFORE_NEWTON = 1
STORAGE_CONVERGENCE_TOL_MM = 0.05
TEMPERATURE_CONVERGENCE_TOL_C = 0.05
SWE_CONVERGENCE_TOL_MM = 0.02
NEWTON_STORAGE_PERTURBATION_MM = 0.10
NEWTON_ENTHALPY_PERTURBATION_J_M3 = 1.5e5
NEWTON_SWE_PERTURBATION_MM = 0.05
FREEZE_ZONE_OSCILLATION_BAND_C = 0.35
MAX_FREEZE_ZONE_ENTHALPY_STEP_J_M3 = 2.5e5
MAX_STORAGE_NEWTON_STEP_MM = 5.0
ANNUAL_SECONDS = 365.2425 * 86400.0
MIN_DEEP_THERMAL_DIFFUSIVITY_M2_S = 0.15e-6
MAX_DEEP_THERMAL_DIFFUSIVITY_M2_S = 1.20e-6
BOUNDARY_SOIL_TEMP_AMPLITUDE_IN_CM = 187.75174655555554
BOUNDARY_SOIL_TEMP_OFFSET_IN_DEG_C = 3.508180888888889


def fit_annual_air_temperature_harmonic(
    dates: pd.Series,
    daily_mean_air_temp_c: np.ndarray,
    days_in_year: int,
) -> tuple[float, float, float]:
    """Fit the annual harmonic used for lower thermal boundary forcing."""
    theta = 2.0 * np.pi * ((pd.to_datetime(dates).dt.dayofyear.to_numpy(dtype=float) - 1.0) / float(days_in_year))
    design_matrix = np.column_stack(
        (
            np.ones_like(theta),
            np.cos(theta),
            np.sin(theta),
        )
    )
    coefficients, _, _, _ = np.linalg.lstsq(design_matrix, np.asarray(daily_mean_air_temp_c, dtype=float), rcond=None)
    return float(coefficients[0]), float(coefficients[1]), float(coefficients[2])


def compute_harmonic_soil_temperature_c(
    mean_temp_c: float,
    cosine_coeff_c: float,
    sine_coeff_c: float,
    depth_cm: np.ndarray | float,
    day_of_year: int,
    days_in_year: int,
    thermal_diffusivity_m2_s: float,
) -> np.ndarray:
    """Compute harmonic soil temperature at depth from annual air forcing."""
    depth_m = np.asarray(depth_cm, dtype=float) / 100.0
    thermal_diffusivity_m2_s = float(
        np.clip(
            thermal_diffusivity_m2_s,
            MIN_DEEP_THERMAL_DIFFUSIVITY_M2_S,
            MAX_DEEP_THERMAL_DIFFUSIVITY_M2_S,
        )
    )
    angular_frequency_rad_s = 2.0 * np.pi / ANNUAL_SECONDS
    damping_depth_m = np.sqrt((2.0 * thermal_diffusivity_m2_s) / max(angular_frequency_rad_s, 1.0e-12))
    damping_depth_m = max(float(damping_depth_m), 0.05)
    phase_lag_rad = depth_m / damping_depth_m
    damping = np.exp(-depth_m / damping_depth_m)
    theta = 2.0 * np.pi * ((float(day_of_year) - 1.0) / float(days_in_year))
    return mean_temp_c + damping * (
        cosine_coeff_c * np.cos(theta - phase_lag_rad)
        + sine_coeff_c * np.sin(theta - phase_lag_rad)
    )


def compute_damped_lower_boundary_temperature_c(
    mean_air_temp_c: float,
    cosine_coeff_c: float,
    sine_coeff_c: float,
    depth_cm: float,
    day_of_year: int,
    days_in_year: int,
) -> float:
    """Evaluate the damped annual lower-boundary temperature."""
    depth_cm = max(float(depth_cm), 0.0)
    damping_depth_cm = max(float(BOUNDARY_SOIL_TEMP_AMPLITUDE_IN_CM), 1.0e-6)
    damping = np.exp(-depth_cm / damping_depth_cm)
    phase_lag_rad = depth_cm / damping_depth_cm
    theta = 2.0 * np.pi * ((float(day_of_year) - 1.0) / float(days_in_year))
    return float(
        mean_air_temp_c
        + float(BOUNDARY_SOIL_TEMP_OFFSET_IN_DEG_C)
        + damping
        * (
            cosine_coeff_c * np.cos(theta - phase_lag_rad)
            + sine_coeff_c * np.sin(theta - phase_lag_rad)
        )
    )


def compute_annual_damping_depth_cm(
    thermal_diffusivity_m2_s: float,
) -> float:
    """Convert thermal diffusivity into annual damping depth."""
    thermal_diffusivity_m2_s = float(
        np.clip(
            thermal_diffusivity_m2_s,
            MIN_DEEP_THERMAL_DIFFUSIVITY_M2_S,
            MAX_DEEP_THERMAL_DIFFUSIVITY_M2_S,
        )
    )
    angular_frequency_rad_s = 2.0 * np.pi / ANNUAL_SECONDS
    damping_depth_m = np.sqrt((2.0 * thermal_diffusivity_m2_s) / max(angular_frequency_rad_s, 1.0e-12))
    return max(float(damping_depth_m * 100.0), 5.0)


def compute_lower_boundary_reference_depth_cm(
    column_bottom_depth_cm: float,
    thermal_diffusivity_m2_s: float,
) -> float:
    """Choose the effective lower-boundary reference depth for the soil column."""
    column_bottom_depth_cm = max(float(column_bottom_depth_cm), 1.0)
    annual_damping_depth_cm = compute_annual_damping_depth_cm(
        thermal_diffusivity_m2_s=thermal_diffusivity_m2_s,
    )
    return max(
        column_bottom_depth_cm + annual_damping_depth_cm,
        2.0 * annual_damping_depth_cm,
    )


def compute_degree_day_melt_factor_mm_per_degC_day(
    daylength_hours: float,
) -> float:
    """Compute the daylength-adjusted snowmelt factor."""
    daylength_hours = float(np.clip(daylength_hours, 0.0, 24.0))
    normalized_daylength = daylength_hours / 24.0
    return 2.0 + (4.0 * normalized_daylength)


def solve_coupled_day(
    storage_mm: np.ndarray,
    soil_temp_c: np.ndarray,
    thickness_mm: np.ndarray,
    mid_cm: np.ndarray,
    porosity: np.ndarray,
    fc: np.ndarray,
    wp: np.ndarray,
    theta_residual: np.ndarray,
    alpha_vg: np.ndarray,
    n_vg: np.ndarray,
    m_vg: np.ndarray,
    ksat_m_per_s: np.ndarray,
    precip_mm: float,
    pet_mm: float,
    season_progress: float,
    root_weights: np.ndarray,
    swe_mm: float,
    snow_depth_mm: float,
    snow_layer_swe_mm: np.ndarray,
    snow_layer_age_days: np.ndarray,
    air_temp_c: float,
    deep_boundary_temp_c: float,
    lower_boundary_distance_m: float,
    deep_drainage_impedence_factor: float,
    melt_factor_mm_per_degC_day: float,
    max_iterations: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    float,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    np.ndarray,
    np.ndarray,
    int,
    int,
    int,
    int,
    int,
    float,
]:
    """Solve one day of coupled water, snow, and temperature dynamics."""
    theta_total_initial = np.clip(np.asarray(storage_mm, dtype=float) / np.asarray(thickness_mm, dtype=float), theta_residual, porosity)
    theta_liquid_initial, theta_ice_initial, _, _, _ = compute_freezing_state(
        theta_total=theta_total_initial,
        theta_sat=porosity,
        theta_r=theta_residual,
        alpha_vg_mpa_inv=alpha_vg,
        n_vg=n_vg,
        m_vg=m_vg,
        soil_temp_c=soil_temp_c,
    )
    current_enthalpy_j_m3 = compute_volumetric_enthalpy_j_m3(
        soil_temp_c=soil_temp_c,
        theta_total=theta_total_initial,
        theta_liquid=theta_liquid_initial,
        theta_ice=theta_ice_initial,
        porosity=porosity,
    )

    def advance_one_substep(
        current_storage_mm: np.ndarray,
        current_enthalpy_j_m3: np.ndarray,
        current_soil_temp_c: np.ndarray,
        current_swe_mm: float,
        current_snow_depth_mm: float,
        current_snow_layer_swe_mm: np.ndarray,
        current_snow_layer_age_days: np.ndarray,
        dt_seconds: float,
        precip_mm: float,
        pet_substep_mm: float,
    ) -> tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        float,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        float,
        float,
        float,
        float,
        float,
        float,
        float,
        np.ndarray,
        np.ndarray,
        int,
        int,
        int,
        float,
    ]:
        storage_min_mm = theta_residual * thickness_mm
        storage_max_mm = porosity * thickness_mm
        n_layers_local = len(current_storage_mm)
        storage_slice = slice(0, n_layers_local)
        enthalpy_slice = slice(n_layers_local, 2 * n_layers_local)
        swe_index = 2 * n_layers_local
        best_state_norm = np.inf
        best_payload = None
        newton_accept_count = 0
        picard_only_count = 0
        damped_newton_count = 0

        def pack_state(
            packed_storage_mm: np.ndarray,
            packed_enthalpy_j_m3: np.ndarray,
            packed_swe_mm: float,
        ) -> np.ndarray:
            return np.concatenate(
                (
                    np.asarray(packed_storage_mm, dtype=float),
                    np.asarray(packed_enthalpy_j_m3, dtype=float),
                    np.array([max(float(packed_swe_mm), 0.0)], dtype=float),
                )
            )

        def unpack_state(
            state_vector: np.ndarray,
        ) -> tuple[np.ndarray, np.ndarray, float]:
            state_vector = np.asarray(state_vector, dtype=float)
            unpacked_storage_mm = np.clip(state_vector[storage_slice], storage_min_mm, storage_max_mm)
            unpacked_enthalpy_j_m3 = state_vector[enthalpy_slice].copy()
            unpacked_swe_mm = max(float(state_vector[swe_index]), 0.0)
            return unpacked_storage_mm, unpacked_enthalpy_j_m3, unpacked_swe_mm

        def compute_scaled_state_norm(residual_vector: np.ndarray) -> float:
            residual_vector = np.asarray(residual_vector, dtype=float)
            storage_norm = float(np.max(np.abs(residual_vector[storage_slice])))
            enthalpy_norm = float(np.max(np.abs(residual_vector[enthalpy_slice])))
            swe_norm = float(abs(residual_vector[swe_index]))
            return max(
                storage_norm / max(STORAGE_CONVERGENCE_TOL_MM, 1e-9),
                enthalpy_norm / max(NEWTON_ENTHALPY_PERTURBATION_J_M3, 1.0),
                swe_norm / max(SWE_CONVERGENCE_TOL_MM, 1e-9),
            )

        def evaluate_state(
            state_vector: np.ndarray,
        ) -> tuple[np.ndarray, np.ndarray, tuple[object, ...], np.ndarray]:
            guess_storage_mm, guess_enthalpy_j_m3, guess_swe_mm = unpack_state(state_vector)
            guess_theta_total = np.clip(guess_storage_mm / thickness_mm, theta_residual, porosity)
            guess_soil_temp_c, _, _, _ = invert_enthalpy_to_temperature(
                enthalpy_j_m3=guess_enthalpy_j_m3,
                theta_total=guess_theta_total,
                theta_sat=porosity,
                theta_residual=theta_residual,
                alpha_vg_mpa_inv=alpha_vg,
                n_vg=n_vg,
                m_vg=m_vg,
                porosity=porosity,
                initial_temp_c=current_soil_temp_c,
            )
            guess_theta_liquid_pre, guess_theta_ice_pre, _, _, _ = compute_freezing_state(
                theta_total=guess_theta_total,
                theta_sat=porosity,
                theta_r=theta_residual,
                alpha_vg_mpa_inv=alpha_vg,
                n_vg=n_vg,
                m_vg=m_vg,
                soil_temp_c=guess_soil_temp_c,
            )
            guess_heat_capacity = compute_volumetric_sensible_heat_capacity_j_m3_k(
                theta_total=guess_theta_total,
                theta_liquid=guess_theta_liquid_pre,
                theta_ice=guess_theta_ice_pre,
                porosity=porosity,
            )
            guess_thermal_diffusivity = compute_soil_thermal_diffusivity_m2_s(
                theta_total=guess_theta_total,
                theta_liquid=guess_theta_liquid_pre,
                theta_ice=guess_theta_ice_pre,
                porosity=porosity,
            )
            top_soil_thermal_conductivity_w_m_k = float(
                max(guess_thermal_diffusivity[0] * guess_heat_capacity[0], 0.05)
            )
            (
                rain_mm,
                snowfall_water_mm,
                snowmelt_mm,
                thaw_degree_c,
                trial_swe_mm,
                trial_snow_depth_mm,
                trial_snow_density_kg_m3,
                trial_snow_layer_swe_mm,
                trial_snow_layer_age_days,
            ) = compute_substep_coupled_snow_balance(
                precip_mm=precip_mm,
                swe_mm=guess_swe_mm,
                snow_depth_mm=current_snow_depth_mm,
                snow_layer_swe_mm=current_snow_layer_swe_mm,
                snow_layer_age_days=current_snow_layer_age_days,
                air_temp_c=air_temp_c,
                top_soil_temp_c=guess_soil_temp_c[0],
                top_soil_thermal_conductivity_w_m_k=top_soil_thermal_conductivity_w_m_k,
                top_soil_half_thickness_m=float(max(0.5 * thickness_mm[0] / 1000.0, 1.0e-6)),
                dt_seconds=dt_seconds,
                melt_factor_mm_per_degC_day=melt_factor_mm_per_degC_day,
            )
            effective_precip_mm = rain_mm + snowmelt_mm

            trial_storage_mm = np.asarray(current_storage_mm, dtype=float).copy()
            trial_storage_mm[0] += max(float(effective_precip_mm), 0.0)
            trial_storage_mm, interlayer_flux_mm, deep_drainage_mm = apply_simple_richards_flow(
                storage_mm=trial_storage_mm,
                thickness_mm=thickness_mm,
                mid_cm=mid_cm,
                soil_temp_c=guess_soil_temp_c,
                theta_sat=porosity,
                theta_fc=fc,
                theta_r=theta_residual,
                alpha_vg_mpa_inv=alpha_vg,
                n_vg=n_vg,
                m_vg=m_vg,
                ksat_m_per_s=ksat_m_per_s,
                deep_drainage_impedence_factor=deep_drainage_impedence_factor,
                dt_seconds=dt_seconds,
            )

            theta_total_trial = np.clip(trial_storage_mm / thickness_mm, theta_residual, porosity)
            theta_liquid_trial, theta_ice_trial, _, _, _ = compute_freezing_state(
                theta_total=theta_total_trial,
                theta_sat=porosity,
                theta_r=theta_residual,
                alpha_vg_mpa_inv=alpha_vg,
                n_vg=n_vg,
                m_vg=m_vg,
                soil_temp_c=guess_soil_temp_c,
            )
            aet_by_layer_mm, _ = compute_root_water_uptake(
                pet_mm=pet_substep_mm,
                growth_stage_progress=season_progress,
                root_weights=root_weights,
                theta_liquid=theta_liquid_trial,
                theta_total=theta_total_trial,
                theta_ice=theta_ice_trial,
                theta_sat=porosity,
                theta_fc=fc,
                theta_stress_floor=theta_residual,
                theta_wp=wp,
                thickness_mm=thickness_mm,
                theta_r=theta_residual,
                alpha_vg_mpa_inv=alpha_vg,
                n_vg=n_vg,
                m_vg=m_vg,
                ksat_m_per_s=ksat_m_per_s,
                dt_seconds=dt_seconds,
            )
            trial_storage_mm = np.clip(trial_storage_mm - aet_by_layer_mm, storage_min_mm, storage_max_mm)
            theta_total_trial = np.clip(trial_storage_mm / thickness_mm, theta_residual, porosity)
            theta_liquid_trial, theta_ice_trial, _, _, _ = compute_freezing_state(
                theta_total=theta_total_trial,
                theta_sat=porosity,
                theta_r=theta_residual,
                alpha_vg_mpa_inv=alpha_vg,
                n_vg=n_vg,
                m_vg=m_vg,
                soil_temp_c=guess_soil_temp_c,
            )
            trial_enthalpy_j_m3, thermal_diffusivity = solve_enthalpy_diffusion_step(
                previous_enthalpy_j_m3=current_enthalpy_j_m3,
                previous_soil_temp_c=current_soil_temp_c,
                theta_total=theta_total_trial,
                theta_liquid=theta_liquid_trial,
                theta_ice=theta_ice_trial,
                porosity=porosity,
                theta_residual=theta_residual,
                alpha_vg_mpa_inv=alpha_vg,
                n_vg=n_vg,
                m_vg=m_vg,
                thickness_mm=thickness_mm,
                swe_mm=trial_swe_mm,
                snow_depth_mm=trial_snow_depth_mm,
                snow_layer_swe_mm=trial_snow_layer_swe_mm,
                snow_layer_age_days=trial_snow_layer_age_days,
                air_temp_c=air_temp_c,
                deep_boundary_temp_c=deep_boundary_temp_c,
                lower_boundary_distance_m=lower_boundary_distance_m,
                dt_seconds=dt_seconds,
            )
            trial_soil_temp_c, theta_liquid_trial, theta_ice_trial, depressed_freezing_temp_c = invert_enthalpy_to_temperature(
                enthalpy_j_m3=trial_enthalpy_j_m3,
                theta_total=theta_total_trial,
                theta_sat=porosity,
                theta_residual=theta_residual,
                alpha_vg_mpa_inv=alpha_vg,
                n_vg=n_vg,
                m_vg=m_vg,
                porosity=porosity,
                initial_temp_c=guess_soil_temp_c,
            )
            latent_heat_capacity = (
                np.maximum(theta_ice_trial, 0.0)
                * 1000.0
                * 3.34e5
                / np.maximum(np.abs(depressed_freezing_temp_c - trial_soil_temp_c) + 0.1, 0.1)
            )
            candidate_state = pack_state(trial_storage_mm, trial_enthalpy_j_m3, trial_swe_mm)
            residual = candidate_state - state_vector
            payload = (
                trial_storage_mm.copy(),
                trial_soil_temp_c.copy(),
                theta_liquid_trial.copy(),
                theta_ice_trial.copy(),
                interlayer_flux_mm.copy(),
                float(deep_drainage_mm),
                aet_by_layer_mm.copy(),
                thermal_diffusivity.copy(),
                latent_heat_capacity.copy(),
                float(trial_swe_mm),
                float(rain_mm),
                float(snowfall_water_mm),
                float(snowmelt_mm),
                float(thaw_degree_c),
                float(trial_snow_depth_mm),
                float(trial_snow_density_kg_m3),
                np.asarray(trial_snow_layer_swe_mm, dtype=float).copy(),
                np.asarray(trial_snow_layer_age_days, dtype=float).copy(),
            )
            return residual, candidate_state, payload, depressed_freezing_temp_c.copy()

        current_state = pack_state(current_storage_mm, current_enthalpy_j_m3, current_swe_mm)

        for nonlinear_iteration in range(1, MAX_COUPLED_NONLINEAR_ITERATIONS + 1):
            residual, candidate_state, payload, depressed_freezing_temp_c = evaluate_state(current_state)
            state_norm = compute_scaled_state_norm(residual)
            if state_norm < best_state_norm:
                best_state_norm = state_norm
                best_payload = payload

            if (
                float(np.max(np.abs(residual[storage_slice]))) <= STORAGE_CONVERGENCE_TOL_MM
                and float(np.max(np.abs(residual[enthalpy_slice]))) <= NEWTON_ENTHALPY_PERTURBATION_J_M3
                and float(abs(residual[swe_index])) <= SWE_CONVERGENCE_TOL_MM
            ):
                assert best_payload is not None
                return (
                    best_payload[0],
                    best_payload[1],
                    best_payload[2],
                    best_payload[3],
                    best_payload[4],
                    best_payload[5],
                    best_payload[6],
                    best_payload[7],
                    best_payload[8],
                    best_payload[9],
                    best_payload[10],
                    best_payload[11],
                    best_payload[12],
                    best_payload[13],
                    best_payload[14],
                    best_payload[15],
                    best_payload[16],
                    best_payload[17],
                    nonlinear_iteration,
                    newton_accept_count,
                    damped_newton_count,
                    best_state_norm,
                )

            if nonlinear_iteration <= PICARD_ITERATIONS_BEFORE_NEWTON:
                picard_only_count += 1
                current_state = candidate_state
                continue

            picard_only_count += 1
            current_state = candidate_state

        assert best_payload is not None
        return (
            best_payload[0],
            best_payload[1],
            best_payload[2],
            best_payload[3],
            best_payload[4],
            best_payload[5],
            best_payload[6],
            best_payload[7],
            best_payload[8],
            best_payload[9],
            best_payload[10],
            best_payload[11],
            best_payload[12],
            best_payload[13],
            best_payload[14],
            best_payload[15],
            best_payload[16],
            best_payload[17],
            MAX_COUPLED_NONLINEAR_ITERATIONS,
            newton_accept_count,
            damped_newton_count,
            best_state_norm,
        )

    current_storage_mm = np.asarray(storage_mm, dtype=float).copy()
    current_soil_temp_c = np.asarray(soil_temp_c, dtype=float).copy()
    thickness_mm = np.asarray(thickness_mm, dtype=float)
    porosity = np.asarray(porosity, dtype=float)
    theta_residual = np.asarray(theta_residual, dtype=float)

    best_storage_mm = current_storage_mm.copy()
    best_soil_temp_c = current_soil_temp_c.copy()
    best_theta_liquid = current_storage_mm / thickness_mm
    best_theta_ice = np.zeros_like(best_theta_liquid)
    best_interlayer_flux_mm = np.zeros(max(len(storage_mm) - 1, 0), dtype=float)
    best_aet_by_layer_mm = np.zeros(len(storage_mm), dtype=float)
    best_deep_drainage_mm = 0.0
    best_thermal_diffusivity = np.full(len(storage_mm), 0.1e-6, dtype=float)
    best_latent_heat_capacity = np.zeros(len(storage_mm), dtype=float)
    best_swe_mm = float(swe_mm)
    best_snow_depth_mm = float(snow_depth_mm)
    current_snow_layer_swe_mm = np.asarray(snow_layer_swe_mm, dtype=float).copy()
    current_snow_layer_age_days = np.asarray(snow_layer_age_days, dtype=float).copy()
    best_snow_layer_swe_mm = current_snow_layer_swe_mm.copy()
    best_snow_layer_age_days = current_snow_layer_age_days.copy()
    best_snow_density_kg_m3 = 1000.0 * (best_swe_mm / max(best_snow_depth_mm, 1.0e-9)) if best_snow_depth_mm > 0.0 else 80.0
    accumulated_rain_mm = 0.0
    accumulated_snowfall_water_mm = 0.0
    accumulated_snowmelt_mm = 0.0
    accumulated_effective_precip_mm = 0.0
    accumulated_thaw_degree_c = 0.0
    accepted_substeps = 0
    nonlinear_iterations_used = 0
    newton_iterations_used = 0
    picard_fallback_iterations_used = 0
    damped_newton_iterations_used = 0
    worst_substep_residual_norm = 0.0
    accumulated_interlayer_flux_mm = np.zeros_like(best_interlayer_flux_mm)
    accumulated_aet_by_layer_mm = np.zeros_like(best_aet_by_layer_mm)
    accumulated_deep_drainage_mm = 0.0
    remaining_seconds = 86400.0
    dt_try_seconds = remaining_seconds

    while remaining_seconds > 1e-6 and accepted_substeps < max_iterations:
        remaining_substeps = max(max_iterations - accepted_substeps, 1)
        minimum_substep_seconds = remaining_seconds / float(remaining_substeps)
        dt_try_seconds = min(dt_try_seconds, remaining_seconds)
        precip_substep_mm = precip_mm * (dt_try_seconds / 86400.0)
        pet_substep_mm = pet_mm * (dt_try_seconds / 86400.0)

        (
            trial_storage_mm,
            trial_soil_temp_c,
            theta_liquid_trial,
            theta_ice_trial,
            interlayer_flux_mm,
            deep_drainage_mm,
            aet_by_layer_mm,
            thermal_diffusivity,
            latent_heat_capacity,
            trial_swe_mm,
            rain_substep_mm,
            snowfall_substep_mm,
            snowmelt_substep_mm,
            thaw_degree_substep_c,
            trial_snow_depth_mm,
            trial_snow_density_kg_m3,
            trial_snow_layer_swe_mm,
            trial_snow_layer_age_days,
            nonlinear_iterations_substep,
            newton_iterations_substep,
            damped_newton_substep,
            substep_residual_norm,
        ) = advance_one_substep(
            current_storage_mm=current_storage_mm,
            current_enthalpy_j_m3=current_enthalpy_j_m3,
            current_soil_temp_c=current_soil_temp_c,
            current_swe_mm=best_swe_mm,
            current_snow_depth_mm=best_snow_depth_mm,
            current_snow_layer_swe_mm=current_snow_layer_swe_mm,
            current_snow_layer_age_days=current_snow_layer_age_days,
            dt_seconds=dt_try_seconds,
            precip_mm=precip_substep_mm,
            pet_substep_mm=pet_substep_mm,
        )

        storage_jump = float(
            np.max(np.abs((trial_storage_mm - current_storage_mm) / np.maximum(thickness_mm, 1.0)))
        )
        temp_jump = float(np.max(np.abs(trial_soil_temp_c - current_soil_temp_c)))
        unstable_storage = storage_jump > 0.08
        unstable_temp = temp_jump > 6.0

        if (unstable_storage or unstable_temp) and dt_try_seconds > minimum_substep_seconds * 1.01:
            dt_try_seconds = max(dt_try_seconds * 0.5, minimum_substep_seconds)
            continue

        current_storage_mm = trial_storage_mm
        current_soil_temp_c = trial_soil_temp_c
        current_enthalpy_j_m3 = compute_volumetric_enthalpy_j_m3(
            soil_temp_c=trial_soil_temp_c,
            theta_total=np.clip(trial_storage_mm / thickness_mm, theta_residual, porosity),
            theta_liquid=theta_liquid_trial,
            theta_ice=theta_ice_trial,
            porosity=porosity,
        )
        best_swe_mm = float(trial_swe_mm)
        best_snow_depth_mm = float(trial_snow_depth_mm)
        best_snow_density_kg_m3 = float(trial_snow_density_kg_m3)
        current_snow_layer_swe_mm = np.asarray(trial_snow_layer_swe_mm, dtype=float).copy()
        current_snow_layer_age_days = np.asarray(trial_snow_layer_age_days, dtype=float).copy()
        best_snow_layer_swe_mm = current_snow_layer_swe_mm.copy()
        best_snow_layer_age_days = current_snow_layer_age_days.copy()
        accumulated_interlayer_flux_mm += interlayer_flux_mm
        accumulated_aet_by_layer_mm += aet_by_layer_mm
        accumulated_deep_drainage_mm += deep_drainage_mm
        accumulated_rain_mm += rain_substep_mm
        accumulated_snowfall_water_mm += snowfall_substep_mm
        accumulated_snowmelt_mm += snowmelt_substep_mm
        accumulated_effective_precip_mm += rain_substep_mm + snowmelt_substep_mm
        accumulated_thaw_degree_c += thaw_degree_substep_c * (dt_try_seconds / 86400.0)
        nonlinear_iterations_used += nonlinear_iterations_substep
        newton_iterations_used += newton_iterations_substep
        damped_newton_iterations_used += damped_newton_substep
        picard_fallback_iterations_used += max(
            nonlinear_iterations_substep - newton_iterations_substep,
            0,
        )
        worst_substep_residual_norm = max(worst_substep_residual_norm, float(substep_residual_norm))
        remaining_seconds -= dt_try_seconds
        accepted_substeps += 1

        best_storage_mm = current_storage_mm.copy()
        best_soil_temp_c = current_soil_temp_c.copy()
        best_theta_liquid = theta_liquid_trial
        best_theta_ice = theta_ice_trial
        best_interlayer_flux_mm = accumulated_interlayer_flux_mm.copy()
        best_aet_by_layer_mm = accumulated_aet_by_layer_mm.copy()
        best_deep_drainage_mm = float(accumulated_deep_drainage_mm)
        best_thermal_diffusivity = thermal_diffusivity
        best_latent_heat_capacity = latent_heat_capacity

        if remaining_seconds <= 1e-6:
            break

        if storage_jump < 0.005 and temp_jump < 0.3:
            dt_try_seconds = remaining_seconds
        elif storage_jump < 0.015 and temp_jump < 0.75:
            dt_try_seconds = min(remaining_seconds, dt_try_seconds * 1.5)
        else:
            dt_try_seconds = min(remaining_seconds, max(dt_try_seconds, minimum_substep_seconds))

    return (
        best_storage_mm,
        best_soil_temp_c,
        best_theta_liquid,
        best_theta_ice,
        best_interlayer_flux_mm,
        best_deep_drainage_mm,
        best_aet_by_layer_mm,
        best_thermal_diffusivity,
        best_latent_heat_capacity,
        accumulated_rain_mm,
        accumulated_snowfall_water_mm,
        accumulated_snowmelt_mm,
        accumulated_effective_precip_mm,
        accumulated_thaw_degree_c,
        best_swe_mm,
        best_snow_depth_mm,
        best_snow_density_kg_m3,
        best_snow_layer_swe_mm,
        best_snow_layer_age_days,
        nonlinear_iterations_used,
        accepted_substeps,
        newton_iterations_used,
        picard_fallback_iterations_used,
        damped_newton_iterations_used,
        worst_substep_residual_norm,
    )


def run_one_site_year_crop(
    group_df: pd.DataFrame,
    soil_df: pd.DataFrame,
    crop_row: dict[str, Any],
    initial_snowpack_depth_in_mm: float,
    deep_drainage_impedence_factor: float,
    previous_year_theta_by_depth: dict[str, float] | None = None,
    previous_input_theta_by_depth: dict[str, float] | None = None,
    previous_year_soil_temp_by_depth: dict[str, float] | None = None,
    previous_year_snow_depth_mm: float | None = None,
) -> pd.DataFrame:
    """Run the model for one site-year-crop combination."""
    group_df = group_df.sort_values("date").reset_index(drop=True)
    soil_df = soil_df.sort_values(["top_cm", "bottom_cm"]).reset_index(drop=True)
    initial_snowpack_depth_in_mm = float(initial_snowpack_depth_in_mm)
    deep_drainage_impedence_factor = float(np.clip(deep_drainage_impedence_factor, 0.0, 1.0))

    n_layers = len(soil_df)

    porosity = soil_df["porosity"].to_numpy(dtype=float)
    fc = soil_df["field_capacity"].to_numpy(dtype=float)
    wp = soil_df["wilting_point"].to_numpy(dtype=float)
    theta = soil_df["initial_soil_moisture"].to_numpy(dtype=float).copy()
    theta_residual = soil_df["theta_residual"].to_numpy(dtype=float)
    alpha_vg = soil_df["van_genuchten_alpha_per_mpa"].to_numpy(dtype=float)
    n_vg = soil_df["van_genuchten_n"].to_numpy(dtype=float)
    m_vg = soil_df["van_genuchten_m"].to_numpy(dtype=float)
    ksat_m_per_s = soil_df["ksat_m_per_s"].to_numpy(dtype=float)

    thickness_mm = soil_df["thickness_mm"].to_numpy(dtype=float)
    mid_cm = soil_df["mid_cm"].to_numpy(dtype=float)
    safe_depth_names = [
        depth_to_safe_name(float(depth_name))
        for depth_name in soil_df["depth"].tolist()
    ]
    if previous_year_theta_by_depth is None:
        previous_year_theta_by_depth = {}
    if previous_input_theta_by_depth is None:
        previous_input_theta_by_depth = {}
    if previous_year_soil_temp_by_depth is None:
        previous_year_soil_temp_by_depth = {}
    for j, safe_depth in enumerate(safe_depth_names):
        if theta[j] < 0.0:
            if safe_depth in previous_year_theta_by_depth:
                theta[j] = float(previous_year_theta_by_depth[safe_depth])
                continue
            if safe_depth in previous_input_theta_by_depth:
                theta[j] = float(previous_input_theta_by_depth[safe_depth])
                continue
            if safe_depth not in previous_year_theta_by_depth:
                raise ValueError(
                    "soil_data.txt has negative initial_soil_moisture but neither previous modeled soil "
                    "moisture nor a prior soil-data seed value is available for "
                    f"site={group_df.loc[0, 'site']}, "
                    f"year={group_df.loc[0, 'year']}, crop={group_df.loc[0, 'crop']}, depth={safe_depth}."
                )

    latitude_deg = float(crop_row["latitude_degrees"])
    alberta_township = float(crop_row["alberta_township"])
    if previous_year_snow_depth_mm is None:
        initial_snowpack_depth = initial_snowpack_depth_in_mm
        if initial_snowpack_depth < 0.0:
            raise ValueError(
                "boundary.txt must provide a non-negative initial_snowpack_depth_in_mm "
                f"for the first simulated year at site={group_df.loc[0, 'site']}, "
                f"year={group_df.loc[0, 'year']}, crop={group_df.loc[0, 'crop']}."
            )
    else:
        initial_snowpack_depth = float(previous_year_snow_depth_mm)

    theta = np.minimum(np.maximum(theta, theta_residual), porosity)

    storage_mm = theta * thickness_mm
    porosity_mm = porosity * thickness_mm
    fc_mm = fc * thickness_mm
    wp_mm = wp * thickness_mm
    residual_mm = theta_residual * thickness_mm

    snow_depth_mm = max(initial_snowpack_depth, 0.0)
    swe_mm = snow_depth_mm * 0.10
    snow_layer_swe_mm, snow_layer_age_days = initialize_snowpack_layers(
        swe_mm=swe_mm,
        snow_depth_mm=snow_depth_mm,
        air_temp_c=float(group_df.loc[0, 'daily_min_air_temp']),
    )
    annual_mean_air_temp_c = float(
        np.mean(0.5 * (group_df["daily_max_air_temp"].to_numpy(dtype=float) + group_df["daily_min_air_temp"].to_numpy(dtype=float)))
    )
    simulation_year = int(group_df.loc[0, "year"])
    days_in_simulation_year = 366 if pd.Timestamp(year=simulation_year, month=1, day=1).is_leap_year else 365
    daily_mean_air_temp_c = 0.5 * (
        group_df["daily_max_air_temp"].to_numpy(dtype=float)
        + group_df["daily_min_air_temp"].to_numpy(dtype=float)
    )
    harmonic_mean_temp_c, harmonic_cosine_coeff_c, harmonic_sine_coeff_c = fit_annual_air_temperature_harmonic(
        dates=group_df["date"],
        daily_mean_air_temp_c=daily_mean_air_temp_c,
        days_in_year=days_in_simulation_year,
    )

    initial_air_temp_mean, _, _ = compute_daily_air_temperature_metrics(
        daily_max_air_temp_c=group_df.loc[0, "daily_max_air_temp"],
        daily_min_air_temp_c=group_df.loc[0, "daily_min_air_temp"],
    )
    theta_liquid_initial, theta_ice_initial, _, _, _ = compute_freezing_state(
        theta_total=theta,
        theta_sat=porosity,
        theta_r=theta_residual,
        alpha_vg_mpa_inv=alpha_vg,
        n_vg=n_vg,
        m_vg=m_vg,
        soil_temp_c=np.full(n_layers, initial_air_temp_mean, dtype=float),
    )
    deep_thermal_diffusivity_m2_s = float(
        np.median(
            compute_soil_thermal_diffusivity_m2_s(
                theta_total=theta,
                theta_liquid=theta_liquid_initial,
                theta_ice=theta_ice_initial,
                porosity=porosity,
            )
        )
    )
    initial_day_of_year = int(pd.Timestamp(group_df.loc[0, "date"]).dayofyear)
    harmonic_initial_soil_temp = compute_harmonic_soil_temperature_c(
        mean_temp_c=harmonic_mean_temp_c,
        cosine_coeff_c=harmonic_cosine_coeff_c,
        sine_coeff_c=harmonic_sine_coeff_c,
        depth_cm=mid_cm,
        day_of_year=initial_day_of_year,
        days_in_year=days_in_simulation_year,
        thermal_diffusivity_m2_s=deep_thermal_diffusivity_m2_s,
    )
    soil_temp = np.asarray(harmonic_initial_soil_temp, dtype=float).copy()
    for j, safe_depth in enumerate(safe_depth_names):
        if safe_depth in previous_year_soil_temp_by_depth:
            soil_temp[j] = float(previous_year_soil_temp_by_depth[safe_depth])
    lower_boundary_reference_depth_cm = compute_lower_boundary_reference_depth_cm(
        column_bottom_depth_cm=float(soil_df["bottom_cm"].iloc[-1]),
        thermal_diffusivity_m2_s=deep_thermal_diffusivity_m2_s,
    )
    lower_boundary_distance_m = max(
        (lower_boundary_reference_depth_cm - float(mid_cm[-1])) / 100.0,
        float(0.5 * thickness_mm[-1] / 1000.0),
    )
    annual_daylength_total = sum(
        compute_daylength_hours(latitude_deg, doy)
        for doy in range(1, days_in_simulation_year + 1)
    )
    results = []

    for _, row in group_df.iterrows():
        date = row["date"]
        daily_air_temp_max_c = float(row["daily_max_air_temp"])
        daily_air_temp_min_c = float(row["daily_min_air_temp"])
        (
            tair,
            daily_air_temp_range_c,
            diurnal_thaw_bias_c,
        ) = compute_daily_air_temperature_metrics(
            daily_max_air_temp_c=daily_air_temp_max_c,
            daily_min_air_temp_c=daily_air_temp_min_c,
        )
        precip_mm = max(float(row["precp"]), 0.0)
        crop_coefficient, growth_stage = compute_stage_crop_coefficient(date, crop_row)
        season_progress = compute_season_progress(date, crop_row)
        rooting_depth_cm = compute_rooting_depth_cm(growth_stage)
        root_weights = compute_root_weights(mid_cm=mid_cm, rooting_depth_cm=rooting_depth_cm)
        bc_p_factor, daylength_hours = compute_blaney_criddle_p_factor(
            date, latitude_deg, annual_daylength_total
        )
        melt_factor_mm_per_degC_day = compute_degree_day_melt_factor_mm_per_degC_day(
            daylength_hours=daylength_hours,
        )
        bc_temperature_term = max((0.46 * tair) + 8.0, 0.0)
        rainfall_cloudiness_damping = 1.0
        reference_et_mm = (
            bc_p_factor
            * bc_temperature_term
            * rainfall_cloudiness_damping
        )

        # FAO training material gives the Blaney-Criddle form
        # ETo = p * (0.46 * Tmean + 8). Here it is combined with the stage-based
        # Kc framework (FAO-56) and evaluated on a daily timestep.
        pet_mm = crop_coefficient * reference_et_mm
        if growth_stage == "off_season":
            seasonal_et_mm = 0.0
            off_season_et_mm = pet_mm
        else:
            seasonal_et_mm = pet_mm
            off_season_et_mm = 0.0
        deep_boundary_temp_c = float(
            compute_damped_lower_boundary_temperature_c(
                mean_air_temp_c=harmonic_mean_temp_c,
                cosine_coeff_c=harmonic_cosine_coeff_c,
                sine_coeff_c=harmonic_sine_coeff_c,
                depth_cm=lower_boundary_reference_depth_cm,
                day_of_year=int(pd.Timestamp(date).dayofyear),
                days_in_year=days_in_simulation_year,
            )
        )
        (
            storage_mm,
            soil_temp,
            theta_liquid,
            theta_ice,
            interlayer_flux_mm,
            deep_drainage_mm,
            aet_by_layer_mm,
            thermal_diffusivity_layer,
            apparent_latent_heat_capacity_layer,
            rain_mm,
            snowfall_water_mm,
            snowmelt_mm,
            effective_precip_mm,
            thaw_degree_c,
            swe_mm,
            snow_depth_mm,
            snow_density_kg_m3,
            snow_layer_swe_mm,
            snow_layer_age_days,
            solver_iterations,
            coupled_substeps_used,
            newton_iterations_used,
            picard_fallback_iterations_used,
            damped_newton_iterations_used,
            worst_substep_residual_norm,
        ) = solve_coupled_day(
            storage_mm=storage_mm,
            soil_temp_c=soil_temp,
            thickness_mm=thickness_mm,
            mid_cm=mid_cm,
            porosity=porosity,
            fc=fc,
            wp=wp,
            theta_residual=theta_residual,
            alpha_vg=alpha_vg,
            n_vg=n_vg,
            m_vg=m_vg,
            ksat_m_per_s=ksat_m_per_s,
            precip_mm=precip_mm,
            pet_mm=pet_mm,
            season_progress=season_progress,
            root_weights=root_weights,
            swe_mm=swe_mm,
            snow_depth_mm=snow_depth_mm,
            snow_layer_swe_mm=snow_layer_swe_mm,
            snow_layer_age_days=snow_layer_age_days,
            air_temp_c=tair,
            deep_boundary_temp_c=deep_boundary_temp_c,
            lower_boundary_distance_m=lower_boundary_distance_m,
            deep_drainage_impedence_factor=deep_drainage_impedence_factor,
            melt_factor_mm_per_degC_day=melt_factor_mm_per_degC_day,
            max_iterations=HYDROLOGY_SUBSTEPS_PER_DAY,
        )

        aet_mm = float(np.sum(aet_by_layer_mm))

        theta = storage_mm / thickness_mm
        _, _, _, _, depressed_freezing_temp_c = compute_freezing_state(
            theta_total=theta,
            theta_sat=porosity,
            theta_r=theta_residual,
            alpha_vg_mpa_inv=alpha_vg,
            n_vg=n_vg,
            m_vg=m_vg,
            soil_temp_c=soil_temp,
        )
        freeze_thaw_damping = np.clip(apparent_latent_heat_capacity_layer / 3.0e8, 0.0, 1.0)

        liquid_storage_mm = theta_liquid * thickness_mm
        ice_storage_mm = theta_ice * thickness_mm

        paw_storage_mm = np.maximum(storage_mm - wp_mm, 0.0)
        out = {
            "site": row["site"],
            "year": row["year"],
            "crop": row["crop"],
            "date": date,
            "alberta_township": alberta_township,
            "latitude_deg": latitude_deg,
            "daily_max_air_temp_c": daily_air_temp_max_c,
            "daily_min_air_temp_c": daily_air_temp_min_c,
            "air_temp_c": tair,
            "daily_air_temp_range_c": daily_air_temp_range_c,
            "surface_diurnal_thaw_bias_c": diurnal_thaw_bias_c,
            "deep_boundary_temp_c": deep_boundary_temp_c,
            "lower_boundary_reference_depth_cm": lower_boundary_reference_depth_cm,
            "precip_mm": precip_mm,
            "rain_mm": rain_mm,
            "snowfall_water_mm": snowfall_water_mm,
            "snowmelt_mm": snowmelt_mm,
            "effective_precip_mm": effective_precip_mm,
            "snow_depth_mm": snow_depth_mm,
            "snow_water_equivalent_mm": swe_mm,
            "snow_density_kg_m3": snow_density_kg_m3,
            "snowmelt_degree_c": thaw_degree_c,
            "daylength_hours": daylength_hours,
            "bc_p_factor_pct": bc_p_factor,
            "bc_temperature_term": bc_temperature_term,
            "rainfall_cloudiness_damping": rainfall_cloudiness_damping,
            "growth_stage": growth_stage,
            "crop_coefficient": crop_coefficient,
            "off_season_kc": crop_row["kc_off_season"],
            "reference_et_mm": reference_et_mm,
            "seasonal_et_mm": seasonal_et_mm,
            "off_season_et_mm": off_season_et_mm,
            "pet_mm": pet_mm,
            "aet_mm": aet_mm,
            "deep_drainage_impedence_factor": deep_drainage_impedence_factor,
            "deep_drainage_mm": deep_drainage_mm,
            "coupled_solver_iterations": solver_iterations,
            "coupled_substeps_used": coupled_substeps_used,
            "newton_iterations_used": newton_iterations_used,
            "picard_fallback_iterations_used": picard_fallback_iterations_used,
            "damped_newton_iterations_used": damped_newton_iterations_used,
            "worst_substep_residual_norm": worst_substep_residual_norm,
        }

        for j, safe_depth in enumerate(safe_depth_names):
            out[f"soil_temp_{safe_depth}_c"] = soil_temp[j]

        for j, safe_depth in enumerate(safe_depth_names):
            out[f"theta_{safe_depth}"] = theta[j]
            out[f"theta_liquid_{safe_depth}"] = theta_liquid[j]
            out[f"theta_ice_{safe_depth}"] = theta_ice[j]
            out[f"theta_residual_{safe_depth}"] = theta_residual[j]

        for j, safe_depth in enumerate(safe_depth_names):
            out[f"storage_{safe_depth}_mm"] = storage_mm[j]
            out[f"liquid_storage_{safe_depth}_mm"] = liquid_storage_mm[j]
            out[f"ice_storage_{safe_depth}_mm"] = ice_storage_mm[j]
            out[f"storage_{safe_depth}_at_porosity_mm"] = porosity_mm[j]
            out[f"storage_{safe_depth}_at_field_capacity_mm"] = fc_mm[j]
            out[f"storage_{safe_depth}_at_wilting_point_mm"] = wp_mm[j]
            out[f"storage_{safe_depth}_at_residual_mm"] = residual_mm[j]
            out[f"plant_available_storage_{safe_depth}_mm"] = paw_storage_mm[j]
            out[f"et_uptake_{safe_depth}_mm"] = aet_by_layer_mm[j]
            out[f"freeze_thaw_damping_{safe_depth}"] = freeze_thaw_damping[j]
            out[f"thermal_diffusivity_{safe_depth}_m2_s"] = thermal_diffusivity_layer[j]
            out[f"apparent_latent_heat_capacity_{safe_depth}_j_m3_k"] = (
                apparent_latent_heat_capacity_layer[j]
            )
            out[f"depressed_freezing_temp_{safe_depth}_c"] = depressed_freezing_temp_c[j]
            out[f"ksat_{safe_depth}_m_per_s"] = ksat_m_per_s[j]
            out[f"van_genuchten_alpha_{safe_depth}_per_mpa"] = alpha_vg[j]
            out[f"van_genuchten_n_{safe_depth}"] = n_vg[j]
            out[f"van_genuchten_m_{safe_depth}"] = m_vg[j]
        for j, safe_depth in enumerate(safe_depth_names[:-1]):
            out[f"vertical_flux_{safe_depth}_to_next_mm"] = interlayer_flux_mm[j]

        results.append(out)

    return pd.DataFrame(results)


def run_model_many_site_years(
    weather_file: str,
    soil_file: str,
    crop_file: str,
    layer_file: str,
    boundary_file: str,
    output_file: str,
) -> pd.DataFrame:
    """Run the full model across all matched site-year-crop combinations."""
    dfw = read_delimited_text_table(weather_file)
    dfs = read_delimited_text_table(soil_file)
    dfc = read_delimited_text_table(crop_file)
    dfb = read_delimited_text_table(boundary_file)
    target_depths = load_output_layer_depths(layer_file)

    required_weather_cols = {
        "site",
        "year",
        "crop",
        "date",
        "precp",
        "daily_max_air_temp",
        "daily_min_air_temp",
    }
    missing_weather = required_weather_cols - set(dfw.columns)
    if missing_weather:
        raise ValueError(f"missing daily_weather.txt columns: {sorted(missing_weather)}")

    dfw = dfw.copy()
    dfw["date"] = parse_weather_dates(dfw["date"])
    bad_temp_order = dfw["daily_max_air_temp"] < dfw["daily_min_air_temp"]
    if bad_temp_order.any():
        bad_rows = dfw.loc[
            bad_temp_order,
            ["site", "year", "crop", "date", "daily_max_air_temp", "daily_min_air_temp"],
        ]
        raise ValueError(
            "daily_weather.txt has rows where daily_max_air_temp < daily_min_air_temp:\n"
            f"{bad_rows.to_string(index=False)}"
        )

    soil_lookup = build_soil_lookup(dfs, target_depths=target_depths)
    crop_lookup = build_crop_lookup(dfc)
    boundary_lookup = build_boundary_lookup(dfb)

    outputs = []
    missing_soil = []
    missing_crop = []
    previous_year_theta_lookup: dict[tuple[str, str], dict[str, float]] = {}
    previous_input_theta_lookup: dict[tuple[str, str, int], dict[str, float]] = {}
    previous_year_soil_temp_lookup: dict[tuple[str, str], dict[str, float]] = {}
    previous_year_snow_depth_lookup: dict[tuple[str, str], float] = {}

    for soil_key in sorted(soil_lookup.keys()):
        site_name, year, crop_name = soil_key
        carry_key = (str(site_name), str(crop_name))
        previous_input_theta_lookup[soil_key] = previous_year_theta_lookup.get(carry_key, {}).copy()
        soil_df = soil_lookup[soil_key].sort_values(["top_cm", "bottom_cm"]).reset_index(drop=True)
        current_seed_theta = {}
        for _, soil_row in soil_df.iterrows():
            safe_depth = depth_to_safe_name(float(soil_row["depth"]))
            theta_value = float(soil_row["initial_soil_moisture"])
            if theta_value >= 0.0:
                current_seed_theta[safe_depth] = theta_value
        if current_seed_theta:
            merged_seed_theta = previous_input_theta_lookup[soil_key].copy()
            merged_seed_theta.update(current_seed_theta)
            previous_year_theta_lookup[carry_key] = merged_seed_theta

    previous_year_theta_lookup = {}

    for key, group_df in dfw.groupby(["site", "year", "crop"], sort=True):
        if key not in soil_lookup:
            missing_soil.append(key)
            continue
        if key not in crop_lookup:
            missing_crop.append(key)
            continue

        site_name = str(key[0])
        crop_name = str(key[2])
        carry_key = (site_name, crop_name)
        if site_name not in boundary_lookup:
            raise ValueError(
                f"boundary.txt is missing a required site row for site={site_name}."
            )
        boundary_row = boundary_lookup[site_name]
        out = run_one_site_year_crop(
            group_df=group_df,
            soil_df=soil_lookup[key],
            crop_row=crop_lookup[key],
            initial_snowpack_depth_in_mm=boundary_row["initial_snowpack_depth_in_mm"],
            deep_drainage_impedence_factor=boundary_row["deep_drainage_impedence_factor"],
            previous_year_theta_by_depth=previous_year_theta_lookup.get(carry_key, {}),
            previous_input_theta_by_depth=previous_input_theta_lookup.get(key, {}),
            previous_year_soil_temp_by_depth=previous_year_soil_temp_lookup.get(carry_key, {}),
            previous_year_snow_depth_mm=previous_year_snow_depth_lookup.get(carry_key),
        )
        outputs.append(out)
        final_row = out.sort_values("date").iloc[-1]
        previous_year_theta_lookup[carry_key] = {
            depth_to_safe_name(float(depth_cm)): float(
                final_row[f"theta_{depth_to_safe_name(float(depth_cm))}"]
            )
            for depth_cm in target_depths
        }
        previous_year_soil_temp_lookup[carry_key] = {
            depth_to_safe_name(float(depth_cm)): float(
                final_row[f"soil_temp_{depth_to_safe_name(float(depth_cm))}_c"]
            )
            for depth_cm in target_depths
        }
        previous_year_snow_depth_lookup[carry_key] = float(final_row["snow_depth_mm"])

    if missing_soil:
        print("\nwarning: these site-year-crop combinations were in weather but not in soil_data.txt:")
        print(pd.DataFrame(missing_soil, columns=["site", "year", "crop"]).to_string(index=False))

    if missing_crop:
        print("\nwarning: these site-year-crop combinations were in weather but not in crop.txt:")
        print(pd.DataFrame(missing_crop, columns=["site", "year", "crop"]).to_string(index=False))

    if not outputs:
        raise ValueError("no matching site-year-crop combinations found across the input files.")

    final_output = pd.concat(outputs, ignore_index=True)
    final_output = final_output.sort_values(["site", "year", "crop", "date"]).reset_index(drop=True)
    final_output = reorder_output_columns(final_output, target_depths=target_depths)
    final_output.to_csv(output_file, index=False)

    return final_output


def main() -> None:
    """Run the default command-line workflow using repository input files."""
    output = run_model_many_site_years(
        weather_file=str(Path(__file__).resolve().parents[2] / "data" / "daily_weather.txt"),
        soil_file=str(Path(__file__).resolve().parents[2] / "data" / "soil_data.txt"),
        crop_file=str(Path(__file__).resolve().parents[2] / "data" / "crop.txt"),
        layer_file=str(Path(__file__).resolve().parents[2] / "data" / "layer.txt"),
        boundary_file=str(Path(__file__).resolve().parents[2] / "data" / "boundary.txt"),
        output_file=str(
            Path(__file__).resolve().parents[2]
            / "outputs"
            / "ab_ag_soil_hydro_temp_model_outputs.txt"
        ),
    )

    print(output.head())
    print("\nfinished. output written to ab_ag_soil_hydro_temp_model_outputs.txt")


if __name__ == "__main__":
    main()
