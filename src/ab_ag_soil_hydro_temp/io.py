from __future__ import annotations

"""Input/output helpers for validating, expanding, and ordering model tables."""

import numpy as np
import pandas as pd

from .common import (
    depth_contains,
    depth_to_safe_name,
    format_bottom_depth_cm,
    normalize_depth,
    parse_bottom_depth_to_cm,
    parse_depth_to_cm,
)
from .water import tune_inflection_point_neg_mpa


def read_delimited_text_table(file_path: str) -> pd.DataFrame:
    """Read a comma- or tab-delimited text table with trimmed column names."""
    candidates = []
    for kwargs in ({}, {"sep": "\t"}):
        try:
            df = pd.read_csv(file_path, **kwargs)
        except Exception:
            continue
        cleaned_columns = [str(col).strip() for col in df.columns]
        score = (len(cleaned_columns), sum(bool(col) for col in cleaned_columns))
        candidates.append((score, cleaned_columns, df))
    if not candidates:
        raise ValueError(f"Could not read delimited text table: {file_path}")
    _, cleaned_columns, best_df = max(candidates, key=lambda item: item[0])
    best_df = best_df.copy()
    best_df.columns = cleaned_columns
    return best_df



def _normalize_soil_input_columns(df_soil: pd.DataFrame) -> pd.DataFrame:
    """Normalize supported legacy soil-input column names."""
    df_soil = df_soil.copy()
    rename_map = {}
    if "Ksat_mm_h" in df_soil.columns and "ksat_mm_h" not in df_soil.columns:
        rename_map["Ksat_mm_h"] = "ksat_mm_h"
    if "theta_residual" in df_soil.columns and "inflection_point_negMPa" not in df_soil.columns:
        rename_map["theta_residual"] = "theta_residual_legacy"
    if "inflection_point_neg_mpa" in df_soil.columns and "inflection_point_negMPa" not in df_soil.columns:
        rename_map["inflection_point_neg_mpa"] = "inflection_point_negMPa"
    return df_soil.rename(columns=rename_map)



def _build_layer_geometry_from_bottom_depths(bottom_depths: list[float]) -> pd.DataFrame:
    """Build the internal layer geometry table from layer bottom depths."""
    bottom_depths = [float(depth_cm) for depth_cm in bottom_depths]
    if not bottom_depths:
        raise ValueError("At least one depth is required.")

    bottom_depth_series = pd.Series(bottom_depths, dtype=float)
    if (bottom_depth_series <= 0.0).any():
        raise ValueError("All bottom depths must be > 0 cm.")
    if not bottom_depth_series.is_monotonic_increasing:
        raise ValueError("Bottom depths must be strictly increasing from the surface downward.")
    if bottom_depth_series.duplicated(keep=False).any():
        duplicates = sorted(bottom_depth_series[bottom_depth_series.duplicated(keep=False)].unique())
        raise ValueError(f"Duplicate bottom depths are not allowed: {duplicates}")

    top_depths = [0.0] + bottom_depths[:-1]
    thicknesses_cm = [bottom - top for top, bottom in zip(top_depths, bottom_depths)]
    if any(thickness_cm <= 0.0 for thickness_cm in thicknesses_cm):
        raise ValueError("Each layer must have positive thickness.")

    layer_df = pd.DataFrame(
        {
            "depth": [format_bottom_depth_cm(depth_cm) for depth_cm in bottom_depths],
            "top_cm": top_depths,
            "bottom_cm": bottom_depths,
            "thickness_cm": thicknesses_cm,
        }
    )
    layer_df["thickness_mm"] = layer_df["thickness_cm"] * 10.0
    layer_df["mid_cm"] = 0.5 * (layer_df["top_cm"] + layer_df["bottom_cm"])
    return layer_df


def load_output_layer_depths(layer_file: str) -> list[float]:
    """Read the requested output layer depths from ``layer.txt``."""
    df_layers = read_delimited_text_table(layer_file)
    if "depth" not in df_layers.columns:
        raise ValueError("layer.txt must contain a 'depth' column.")

    depth_values = [value for value in df_layers["depth"].tolist() if str(value).strip()]
    if not depth_values:
        raise ValueError("layer.txt must contain at least one non-empty comma-separated depth.")

    bottom_depths = [parse_bottom_depth_to_cm(value) for value in depth_values]
    _build_layer_geometry_from_bottom_depths(bottom_depths)
    return bottom_depths



def _carry_forward_soil_inputs(df_soil: pd.DataFrame) -> pd.DataFrame:
    """
    Resolve year-to-year hydraulic inputs, treating negative values as sentinels.

    Negative hydraulic inputs carry forward the previous year's value for the same
    site-depth sequence. A special case is applied when every `ksat_mm_h` entry in
    a site-depth sequence is negative: the first year's `ksat_mm_h` is estimated
    from porosity, field capacity, wilting point, and the required first-year
    `inflection_point_negMPa` using a Mualem-van Genuchten inversion anchored at
    `K(field capacity) = 0.1 mm day-1`. The van Genuchten hydraulic parameters are
    solved once from the first resolved year for each site-depth and then reused
    for all subsequent years at that site-depth.
    """
    df_soil = _normalize_soil_input_columns(df_soil)
    df_soil["depth"] = df_soil["depth"].map(normalize_depth)

    def estimate_ksat_mm_h_from_fc_porosity(
        theta_sat: float,
        theta_fc: float,
        theta_wp: float,
        initial_inflection_point_neg_mpa: float,
        conductivity_at_fc_mm_day: float = 0.1,
    ) -> float:
        """
        Estimate saturated hydraulic conductivity from porosity and field capacity.

        A constrained van Genuchten retention curve is tuned from porosity, field
        capacity, wilting point, and a required positive inflection-point suction.
        The Mualem relative conductivity at field capacity is then inverted from an
        assumed conductivity at field capacity. The result is returned in mm h-1.
        """
        theta_sat = float(theta_sat)
        theta_fc = float(theta_fc)
        theta_wp = float(theta_wp)
        initial_inflection_point_neg_mpa = float(initial_inflection_point_neg_mpa)
        if (
            not np.isfinite(theta_sat)
            or not np.isfinite(theta_fc)
            or not np.isfinite(theta_wp)
            or not np.isfinite(initial_inflection_point_neg_mpa)
        ):
            raise ValueError(
                "Cannot estimate ksat_mm_h from non-finite porosity, field_capacity, "
                "wilting_point, or inflection_point_negMPa."
            )
        if (
            theta_sat <= 0.0
            or theta_fc <= 0.0
            or theta_fc >= theta_sat
            or theta_wp <= 0.0
            or theta_wp >= theta_fc
            or initial_inflection_point_neg_mpa <= 0.0
        ):
            raise ValueError(
                "Cannot estimate ksat_mm_h unless porosity, field_capacity, wilting_point, "
                "and inflection_point_negMPa satisfy 0 < wilting_point < field_capacity < "
                "porosity and inflection_point_negMPa > 0. "
                f"Got porosity={theta_sat}, field_capacity={theta_fc}, "
                f"wilting_point={theta_wp}, inflection_point_negMPa={initial_inflection_point_neg_mpa}."
            )

        _, theta_r, _, n_vg = tune_inflection_point_neg_mpa(
            theta_sat=theta_sat,
            theta_fc=theta_fc,
            theta_wp=theta_wp,
            initial_inflection_point_neg_mpa=initial_inflection_point_neg_mpa,
            n_upper=10.0,
        )
        m_vg = float(1.0 - (1.0 / max(n_vg, 1.000001)))
        effective_saturation_at_fc = np.clip(
            (theta_fc - theta_r) / max(theta_sat - theta_r, 1.0e-9),
            1.0e-9,
            0.999999,
        )
        relative_conductivity_at_fc = (
            effective_saturation_at_fc ** 0.5
        ) * (1.0 - (1.0 - effective_saturation_at_fc ** (1.0 / max(m_vg, 1.0e-9))) ** m_vg) ** 2
        ksat_mm_day = conductivity_at_fc_mm_day / max(relative_conductivity_at_fc, 1.0e-12)
        return float(max(ksat_mm_day / 24.0, 0.0))

    hydraulic_input_cols = [
        "porosity",
        "field_capacity",
        "wilting_point",
        "inflection_point_negMPa",
        "ksat_mm_h",
    ]
    out_groups = []
    for _, group in df_soil.groupby(["site", "depth"], sort=False):
        group = group.sort_values("year").reset_index(drop=True).copy()
        all_years_negative_ksat = bool((pd.to_numeric(group["ksat_mm_h"], errors="coerce") < 0.0).all())
        previous_values: dict[str, float] = {}
        fixed_hydraulic_parameters: dict[str, float] | None = None
        fixed_inflection_point: float | None = None

        resolved_rows = []
        for row_index, row in group.iterrows():
            row = row.copy()

            raw_negative = False
            for col in hydraulic_input_cols:
                value = float(row[col])
                if value < 0.0:
                    if col == "ksat_mm_h" and all_years_negative_ksat and col not in previous_values:
                        value = estimate_ksat_mm_h_from_fc_porosity(
                            theta_sat=float(row["porosity"]),
                            theta_fc=float(row["field_capacity"]),
                            theta_wp=float(row["wilting_point"]),
                            initial_inflection_point_neg_mpa=float(row["inflection_point_negMPa"]),
                        )
                    else:
                        raw_negative = True
                    if value < 0.0 and col not in previous_values:
                        raise ValueError(
                            f"soil_data.txt has negative {col} without a previous year's value for "
                            f"site={row['site']}, crop={row['crop']}, depth={row['depth']}, year={row['year']}."
                        )
                    if value < 0.0:
                        value = previous_values[col]
                previous_values[col] = value
                row[col] = value

            row["used_previous_year_hydraulics"] = bool(raw_negative)

            inflection_value = float(row["inflection_point_negMPa"])
            if not np.isfinite(inflection_value) or inflection_value <= 0.0:
                raise ValueError(
                    "soil_data.txt requires a positive inflection_point_negMPa for the first resolved "
                    f"year of each site-depth. Got inflection_point_negMPa={row['inflection_point_negMPa']} "
                    f"for site={row['site']}, crop={row['crop']}, depth={row['depth']}, year={row['year']}."
                )

            if row_index == 0:
                theta_sat = float(row["porosity"])
                theta_fc = float(row["field_capacity"])
                theta_wp = float(row["wilting_point"])

                tuned_inflection, theta_r, alpha_vg, n_vg = tune_inflection_point_neg_mpa(
                    theta_sat=theta_sat,
                    theta_fc=theta_fc,
                    theta_wp=theta_wp,
                    initial_inflection_point_neg_mpa=inflection_value,
                    n_upper=10.0,
                )
                m_vg = float(1.0 - (1.0 / max(n_vg, 1.000001)))

                row["inflection_point_negMPa"] = tuned_inflection
                previous_values["inflection_point_negMPa"] = float(tuned_inflection)
                row["theta_residual"] = float(theta_r)
                row["van_genuchten_alpha_per_mpa"] = float(alpha_vg)
                row["van_genuchten_n"] = float(n_vg)
                row["van_genuchten_m"] = float(m_vg)
                fixed_inflection_point = float(tuned_inflection)
                fixed_hydraulic_parameters = {
                    "theta_residual": float(theta_r),
                    "van_genuchten_alpha_per_mpa": float(alpha_vg),
                    "van_genuchten_n": float(n_vg),
                    "van_genuchten_m": float(m_vg),
                }
            else:
                if fixed_hydraulic_parameters is None or fixed_inflection_point is None:
                    raise ValueError(
                        "Missing first-year site-depth hydraulic solution before processing later years "
                        f"for site={row['site']}, crop={row['crop']}, depth={row['depth']}, year={row['year']}."
                    )
                row["inflection_point_negMPa"] = fixed_inflection_point
                previous_values["inflection_point_negMPa"] = float(fixed_inflection_point)
                for key, value in fixed_hydraulic_parameters.items():
                    row[key] = value

            resolved_rows.append(row)

        out_groups.append(pd.DataFrame(resolved_rows))

    return pd.concat(out_groups, ignore_index=True)



def build_soil_lookup(
    df_soil: pd.DataFrame,
    target_depths: list[float],
) -> dict[tuple[str, int, str], pd.DataFrame]:
    """Expand source soil horizons onto the model's target layers."""
    df_soil = _normalize_soil_input_columns(df_soil)
    required_cols = {
        "site",
        "year",
        "crop",
        "depth",
        "porosity",
        "field_capacity",
        "wilting_point",
        "inflection_point_negMPa",
        "initial_soil_moisture",
        "ksat_mm_h",
    }
    missing = required_cols - set(df_soil.columns)
    if missing:
        raise ValueError(f"missing soil_data.txt columns: {sorted(missing)}")

    df_soil = _carry_forward_soil_inputs(df_soil)
    target_depth_df = _build_layer_geometry_from_bottom_depths(target_depths)
    soil_lookup = {}

    for key, g in df_soil.groupby(["site", "year", "crop"], sort=True):
        g = g.copy().reset_index(drop=True)
        depth_has_range = g["depth"].str.contains("-", regex=False)
        if depth_has_range.any() and not depth_has_range.all():
            raise ValueError(
                f"soil_data.txt mixes depth ranges and bottom-depth-only values for site-year-crop={key}."
            )

        if depth_has_range.all():
            g[["top_cm", "bottom_cm"]] = g["depth"].apply(
                lambda value: pd.Series(parse_depth_to_cm(value))
            )
            g["depth"] = g["bottom_cm"].map(format_bottom_depth_cm)
            dupes = g.duplicated(subset=["bottom_cm"], keep=False)
        else:
            g["bottom_cm"] = g["depth"].map(parse_bottom_depth_to_cm)
            dupes = g.duplicated(subset=["bottom_cm"], keep=False)
            g = g.sort_values("bottom_cm").reset_index(drop=True)
            g["top_cm"] = [0.0] + g["bottom_cm"].tolist()[:-1]
            if ((g["bottom_cm"] - g["top_cm"]) <= 0.0).any():
                raise ValueError(
                    f"soil_data.txt has non-positive layer thickness for site-year-crop={key}."
                )
            g["depth"] = g["bottom_cm"].map(format_bottom_depth_cm)

        if dupes.any():
            bad = g.loc[dupes, ["site", "year", "crop", "depth"]]
            raise ValueError(
                "soil_data.txt has duplicate site-year-crop-depth rows:\n"
                f"{bad.to_string(index=False)}"
            )

        expanded_rows = []
        for _, thin in target_depth_df.iterrows():
            matches = g[
                g.apply(
                    lambda r: depth_contains(
                        r["top_cm"], r["bottom_cm"],
                        thin["top_cm"], thin["bottom_cm"], 1e-9
                    ),
                    axis=1,
                )
            ]
            if len(matches) == 0:
                raise ValueError(
                    f"soil_data.txt does not contain a source layer covering depth "
                    f"{thin['depth']} for site-year-crop={key}"
                )
            if len(matches) > 1:
                matches = matches.copy()
                matches["source_thickness"] = matches["bottom_cm"] - matches["top_cm"]
                src = matches.sort_values("source_thickness").iloc[0]
            else:
                src = matches.iloc[0]

            theta_fc = float(src["field_capacity"])
            theta_wp = float(src["wilting_point"])
            theta_sat = float(src["porosity"])
            theta_r = float(src["theta_residual"])
            alpha_vg = float(src["van_genuchten_alpha_per_mpa"])
            n_vg = float(src["van_genuchten_n"])
            m_vg = float(src["van_genuchten_m"])
            ksat_mm_h = max(float(src["ksat_mm_h"]), 0.0)
            ksat_m_per_s = ksat_mm_h * 1.0e-3 / 3600.0

            expanded_rows.append(
                {
                    "depth": thin["depth"],
                    "top_cm": thin["top_cm"],
                    "bottom_cm": thin["bottom_cm"],
                    "thickness_cm": thin["thickness_cm"],
                    "thickness_mm": thin["thickness_mm"],
                    "mid_cm": thin["mid_cm"],
                    "porosity": theta_sat,
                    "field_capacity": theta_fc,
                    "wilting_point": theta_wp,
                    "inflection_point_negMPa": float(src["inflection_point_negMPa"]),
                    "initial_soil_moisture": float(src["initial_soil_moisture"]),
                    "theta_residual": theta_r,
                    "van_genuchten_alpha_per_mpa": alpha_vg,
                    "van_genuchten_n": n_vg,
                    "van_genuchten_m": m_vg,
                    "ksat_mm_h": ksat_mm_h,
                    "ksat_m_per_s": ksat_m_per_s,
                    "used_previous_year_hydraulics": bool(src.get("used_previous_year_hydraulics", False)),
                }
            )

        soil_lookup[key] = pd.DataFrame(expanded_rows)

    return soil_lookup


def build_boundary_lookup(df_boundary: pd.DataFrame) -> dict[str, dict[str, float]]:
    """Validate site boundary conditions and return a site lookup."""
    required_cols = {
        "site",
        "drainage",
        "initial_snowpack_depth_in_mm",
    }
    missing = required_cols - set(df_boundary.columns)
    if missing:
        raise ValueError(f"missing boundary.txt columns: {sorted(missing)}")

    df_boundary = df_boundary.copy()
    df_boundary["site"] = df_boundary["site"].astype(str).str.strip()
    numeric_cols = [
        "drainage",
        "initial_snowpack_depth_in_mm",
    ]
    for col in numeric_cols:
        df_boundary[col] = pd.to_numeric(df_boundary[col], errors="coerce")

    if (df_boundary["site"] == "").any():
        raise ValueError("boundary.txt has blank site names.")
    dupes = df_boundary.duplicated(subset=["site"], keep=False)
    if dupes.any():
        bad = df_boundary.loc[
            dupes,
            ["site", "drainage", "initial_snowpack_depth_in_mm"],
        ]
        raise ValueError(
            "boundary.txt has duplicate site rows:\n"
            f"{bad.to_string(index=False)}"
        )
    for col in numeric_cols:
        if df_boundary[col].isna().any():
            bad = df_boundary.loc[
                df_boundary[col].isna(),
                ["site", "drainage", "initial_snowpack_depth_in_mm"],
            ]
            raise ValueError(
                f"boundary.txt has non-numeric {col} values:\n"
                f"{bad.to_string(index=False)}"
            )
    bad_drainage = (df_boundary["drainage"] < 0.0) | (df_boundary["drainage"] > 1.0)
    if bad_drainage.any():
        bad = df_boundary.loc[
            bad_drainage,
            ["site", "drainage", "initial_snowpack_depth_in_mm"],
        ]
        raise ValueError(
            "boundary.txt drainage values must be between 0 and 1:\n"
            f"{bad.to_string(index=False)}"
        )
    bad_snow = df_boundary["initial_snowpack_depth_in_mm"] < 0.0
    if bad_snow.any():
        bad = df_boundary.loc[
            bad_snow,
            ["site", "drainage", "initial_snowpack_depth_in_mm"],
        ]
        raise ValueError(
            "boundary.txt initial_snowpack_depth_in_mm values must be >= 0:\n"
            f"{bad.to_string(index=False)}"
        )

    return {
        str(row["site"]): {
            "deep_drainage_impedence_factor": float(row["drainage"]),
            "initial_snowpack_depth_in_mm": float(row["initial_snowpack_depth_in_mm"]),
        }
        for _, row in df_boundary.iterrows()
    }


def reorder_output_columns(
    df: pd.DataFrame,
    target_depths: list[float],
) -> pd.DataFrame:
    """Arrange output columns into a stable publication-friendly order."""
    id_cols = [
        "site",
        "year",
        "crop",
        "date",
        "alberta_township",
        "latitude_deg",
        "daily_max_air_temp_c",
        "daily_min_air_temp_c",
        "air_temp_c",
        "daily_air_temp_range_c",
        "surface_diurnal_thaw_bias_c",
        "precip_mm",
        "rain_mm",
        "snowfall_water_mm",
        "snowmelt_mm",
        "effective_precip_mm",
        "snow_depth_mm",
        "snow_water_equivalent_mm",
        "snow_density_kg_m3",
        "snowmelt_degree_c",
        "daylength_hours",
        "bc_p_factor_pct",
        "bc_temperature_term",
        "rainfall_cloudiness_damping",
        "growth_stage",
        "crop_coefficient",
        "off_season_kc",
        "reference_et_mm",
        "seasonal_et_mm",
        "off_season_et_mm",
        "pet_mm",
        "aet_mm",
        "deep_drainage_impedence_factor",
        "deep_drainage_mm",
        "coupled_solver_iterations",
        "coupled_substeps_used",
        "newton_iterations_used",
        "picard_fallback_iterations_used",
        "damped_newton_iterations_used",
        "worst_substep_residual_norm",
    ]

    depth_order = [depth_to_safe_name(depth_cm) for depth_cm in target_depths]

    soil_temp_cols = [f"soil_temp_{d}_c" for d in depth_order if f"soil_temp_{d}_c" in df.columns]

    theta_cols = [f"theta_{d}" for d in depth_order if f"theta_{d}" in df.columns]
    theta_liquid_cols = [f"theta_liquid_{d}" for d in depth_order if f"theta_liquid_{d}" in df.columns]
    theta_ice_cols = [f"theta_ice_{d}" for d in depth_order if f"theta_ice_{d}" in df.columns]
    theta_residual_cols = [f"theta_residual_{d}" for d in depth_order if f"theta_residual_{d}" in df.columns]

    storage_cols = [f"storage_{d}_mm" for d in depth_order if f"storage_{d}_mm" in df.columns]
    liquid_storage_cols = [
        f"liquid_storage_{d}_mm" for d in depth_order if f"liquid_storage_{d}_mm" in df.columns
    ]
    ice_storage_cols = [
        f"ice_storage_{d}_mm" for d in depth_order if f"ice_storage_{d}_mm" in df.columns
    ]
    storage_fc_cols = [
        f"storage_{d}_at_field_capacity_mm" for d in depth_order
        if f"storage_{d}_at_field_capacity_mm" in df.columns
    ]
    storage_porosity_cols = [
        f"storage_{d}_at_porosity_mm" for d in depth_order
        if f"storage_{d}_at_porosity_mm" in df.columns
    ]
    storage_wp_cols = [
        f"storage_{d}_at_wilting_point_mm" for d in depth_order
        if f"storage_{d}_at_wilting_point_mm" in df.columns
    ]
    storage_residual_cols = [
        f"storage_{d}_at_residual_mm" for d in depth_order
        if f"storage_{d}_at_residual_mm" in df.columns
    ]
    vg_alpha_cols = [
        f"van_genuchten_alpha_{d}_per_mpa" for d in depth_order
        if f"van_genuchten_alpha_{d}_per_mpa" in df.columns
    ]
    vg_n_cols = [f"van_genuchten_n_{d}" for d in depth_order if f"van_genuchten_n_{d}" in df.columns]
    vg_m_cols = [f"van_genuchten_m_{d}" for d in depth_order if f"van_genuchten_m_{d}" in df.columns]
    psi_inflection_cols = [
        f"psi_inflection_{d}_mpa" for d in depth_order
        if f"psi_inflection_{d}_mpa" in df.columns
    ]
    paw_storage_cols = [
        f"plant_available_storage_{d}_mm" for d in depth_order
        if f"plant_available_storage_{d}_mm" in df.columns
    ]
    et_uptake_cols = [
        f"et_uptake_{d}_mm" for d in depth_order if f"et_uptake_{d}_mm" in df.columns
    ]
    thermal_diffusivity_cols = [
        f"thermal_diffusivity_{d}_m2_s" for d in depth_order
        if f"thermal_diffusivity_{d}_m2_s" in df.columns
    ]
    apparent_latent_heat_capacity_cols = [
        f"apparent_latent_heat_capacity_{d}_j_m3_k" for d in depth_order
        if f"apparent_latent_heat_capacity_{d}_j_m3_k" in df.columns
    ]
    depressed_freezing_temp_cols = [
        f"depressed_freezing_temp_{d}_c" for d in depth_order
        if f"depressed_freezing_temp_{d}_c" in df.columns
    ]
    ksat_cols = [
        f"ksat_{d}_m_per_s" for d in depth_order
        if f"ksat_{d}_m_per_s" in df.columns
    ]
    vertical_flux_cols = [
        f"vertical_flux_{depth_order[i]}_to_next_mm"
        for i in range(max(len(depth_order) - 1, 0))
        if f"vertical_flux_{depth_order[i]}_to_next_mm" in df.columns
    ]
    freeze_thaw_damping_cols = [
        f"freeze_thaw_damping_{d}" for d in depth_order
        if f"freeze_thaw_damping_{d}" in df.columns
    ]
    ordered = (
        [c for c in id_cols if c in df.columns]
        + soil_temp_cols
        + theta_cols
        + theta_liquid_cols
        + theta_ice_cols
        + theta_residual_cols
        + storage_cols
        + liquid_storage_cols
        + ice_storage_cols
        + storage_porosity_cols
        + storage_fc_cols
        + storage_wp_cols
        + storage_residual_cols
        + paw_storage_cols
        + et_uptake_cols
        + thermal_diffusivity_cols
        + apparent_latent_heat_capacity_cols
        + depressed_freezing_temp_cols
        + ksat_cols
        + vg_alpha_cols
        + vg_n_cols
        + vg_m_cols
        + psi_inflection_cols
        + vertical_flux_cols
        + freeze_thaw_damping_cols
    )

    remaining = [c for c in df.columns if c not in ordered]
    return df[ordered + remaining]
