from __future__ import annotations

"""Crop and astronomical helper functions used by the daily site simulations.

References:
- Allen et al. (1998), FAO-56, for daylength and reference ET support terms.
- FAO Blaney-Criddle training material for the daylight percentage term.
"""

from typing import Any

import numpy as np
import pandas as pd

from .common import (
    ALBERTA_SOUTH_BORDER_LATITUDE_DEG,
    MILES_PER_LATITUDE_DEGREE,
    TOWNSHIP_HEIGHT_MILES,
)


def township_to_latitude_degrees(alberta_township: float) -> float:
    """Approximate latitude from Alberta township number."""
    township = float(alberta_township)
    if township < 1:
        raise ValueError("alberta_township must be >= 1.")

    # Alberta DLS townships are 6 miles tall and counted north from 49 N.
    # Here the latitude is approximated as the township midpoint.
    # Source: Alberta Energy Regulator DLS help.
    township_height_deg = TOWNSHIP_HEIGHT_MILES / MILES_PER_LATITUDE_DEGREE
    return ALBERTA_SOUTH_BORDER_LATITUDE_DEG + (township - 0.5) * township_height_deg


def compute_daylength_hours(latitude_deg: float, day_of_year: int) -> float:
    """Compute astronomical daylength from latitude and day of year."""
    # FAO-56 solar declination, sunset hour angle, and daylight hours
    # relationships (Allen et al., 1998; Ch. 3-4, Eq. 24-25 and N = 24/pi * ws).
    latitude_rad = np.deg2rad(latitude_deg)
    declination = 0.409 * np.sin((2.0 * np.pi * day_of_year / 365.0) - 1.39)

    sunset_hour_angle_arg = -np.tan(latitude_rad) * np.tan(declination)
    sunset_hour_angle_arg = np.clip(sunset_hour_angle_arg, -1.0, 1.0)
    sunset_hour_angle = np.arccos(sunset_hour_angle_arg)

    return float(24.0 * sunset_hour_angle / np.pi)


def compute_blaney_criddle_p_factor(
    date_value: object,
    latitude_deg: float,
    annual_daylength_total: float | None,
) -> tuple[float, float]:
    """Return the daily Blaney-Criddle daylight percentage and daylength."""
    # The daily p term is implemented as the day's fraction of annual daylight
    # hours so that p sums to 100 over the year. This is a daily-scale extension
    # of the Blaney-Criddle monthly "mean daily percentage of annual daytime
    # hours" term described in FAO training material.
    timestamp = pd.Timestamp(date_value)
    day_of_year = int(timestamp.dayofyear)
    days_in_year = 366 if timestamp.is_leap_year else 365

    daylength_today = compute_daylength_hours(latitude_deg, day_of_year)
    if annual_daylength_total is None:
        annual_daylength_total = sum(
            compute_daylength_hours(latitude_deg, doy)
            for doy in range(1, days_in_year + 1)
        )

    return 100.0 * daylength_today / annual_daylength_total, daylength_today


def parse_crop_stage_start_date(raw_value: object, year_value: int, field_name: str) -> pd.Timestamp:
    """Parse crop-stage dates supplied without an explicit year."""
    timestamp = pd.to_datetime(raw_value, errors="coerce")
    if pd.isna(timestamp):
        for fmt in ("%d-%b", "%d-%B", "%b-%d", "%B-%d"):
            try:
                timestamp = pd.to_datetime(str(raw_value), format=fmt, errors="raise")
                break
            except Exception:
                continue
    if pd.isna(timestamp):
        raise ValueError(
            f"could not parse {field_name}={raw_value!r} as a date for crop stage timing."
        )

    return pd.Timestamp(year=int(year_value), month=int(timestamp.month), day=int(timestamp.day))


def compute_stage_crop_coefficient(
    date_value: object,
    crop_row: dict[str, Any],
) -> tuple[float, str]:
    """Interpolate the crop coefficient across the configured growth stages."""
    # The crop curve below keeps the same CSV inputs but maps them to
    # phenology-based phases:
    #   - seeding to emergence: low Kc
    #   - emergence to anthesis: gradual increase to peak Kc
    #   - around anthesis: peak Kc
    #   - anthesis to early grain fill: hold high Kc
    #   - grain fill to maturity: begin decreasing Kc
    #   - maturity to harvest: low Kc again
    #
    current_date = pd.Timestamp(date_value).normalize()

    season_start = crop_row["season_start_date"]
    emergence_date = crop_row["emergence_date"]
    anthesis_date = crop_row["anthesis_date"]
    early_grain_fill_date = crop_row["early_grain_fill_date"]
    maturity_date = crop_row["maturity_date"]
    harvest_date = crop_row["harvest_date"]

    seeding_emergence_end = emergence_date
    emergence_anthesis_end = anthesis_date
    anthesis_peak_end = min(
        anthesis_date + pd.Timedelta(days=7),
        early_grain_fill_date,
    )
    early_grain_fill_end = early_grain_fill_date
    grain_fill_maturity_end = maturity_date
    harvest_end = harvest_date

    kc_initial = float(crop_row["kc_initial"])
    kc_mid = float(crop_row["kc_mid"])
    kc_end = float(crop_row["kc_end"])
    kc_off_season = float(crop_row["kc_off_season"])

    if current_date < season_start or current_date >= harvest_end:
        return kc_off_season, "off_season"

    if current_date < seeding_emergence_end:
        return kc_initial, "seeding_to_emergence"

    if current_date < emergence_anthesis_end:
        development_days = max((emergence_anthesis_end - seeding_emergence_end).days, 1)
        progress = (current_date - seeding_emergence_end).days / development_days
        kc_value = kc_initial + progress * (kc_mid - kc_initial)
        return kc_value, "emergence_to_anthesis"

    if current_date < anthesis_peak_end:
        return kc_mid, "anthesis_peak"

    if current_date < early_grain_fill_end:
        return kc_mid, "anthesis_to_early_grain_fill"

    if current_date < grain_fill_maturity_end:
        grain_fill_to_maturity_days = max(
            (grain_fill_maturity_end - early_grain_fill_end).days,
            1,
        )
        progress = (current_date - early_grain_fill_end).days / grain_fill_to_maturity_days
        kc_value = kc_mid + progress * (kc_end - kc_mid)
        return kc_value, "grain_fill_to_maturity"

    return kc_end, "maturity_to_harvest"


def compute_season_progress(
    date_value: object,
    crop_row: dict[str, Any],
) -> float:
    """Map a date onto a 0 to 1 seasonal progress scale."""
    current_date = pd.Timestamp(date_value).normalize()
    season_start = crop_row["season_start_date"]
    harvest_date = crop_row["harvest_date"]

    if current_date <= season_start:
        return 0.0
    if current_date >= harvest_date:
        return 1.0

    season_duration_days = max((harvest_date - season_start).days, 1)
    elapsed_days = (current_date - season_start).days
    return float(np.clip(elapsed_days / season_duration_days, 0.0, 1.0))


def build_crop_lookup(df_crop: pd.DataFrame) -> dict[tuple[str, int, str], dict[str, Any]]:
    """Validate crop metadata and build a site-year-crop lookup table."""
    # Stage timing and Kc values are read from crop.txt so crop-specific
    # assumptions remain data-driven rather than hard-coded in the model.
    required_cols = {
        "site",
        "year",
        "crop",
        "alberta_township",
        "season_start_date",
        "emergence_date",
        "anthesis_date",
        "early_grain_fill_date",
        "maturity_date",
        "harvest_date",
        "kc_initial",
        "kc_mid",
        "kc_end",
        "kc_off_season",
    }
    missing = required_cols - set(df_crop.columns)
    if missing:
        raise ValueError(f"missing crop.txt columns: {sorted(missing)}")

    df_crop = df_crop.dropna(how="all").copy()

    dupes = df_crop.duplicated(subset=["site", "year", "crop"], keep=False)
    if dupes.any():
        bad = df_crop.loc[dupes, ["site", "year", "crop"]]
        raise ValueError(
            "crop.txt has duplicate site-year-crop rows:\n"
            f"{bad.to_string(index=False)}"
        )

    crop_lookup = {}
    for _, row in df_crop.iterrows():
        season_start_date = parse_crop_stage_start_date(
            row["season_start_date"], row["year"], "season_start_date"
        )
        parsed_required_dates = {
            field_name: parse_crop_stage_start_date(
                row[field_name],
                row["year"],
                field_name,
            )
            for field_name in [
                "emergence_date",
                "anthesis_date",
                "early_grain_fill_date",
                "maturity_date",
                "harvest_date",
            ]
        }
        stage_dates = [
            season_start_date,
            parsed_required_dates["emergence_date"],
            parsed_required_dates["anthesis_date"],
            parsed_required_dates["early_grain_fill_date"],
            parsed_required_dates["maturity_date"],
            parsed_required_dates["harvest_date"],
        ]
        if stage_dates != sorted(stage_dates):
            raise ValueError(
                "crop stage dates must be chronological for "
                f"site-year-crop={(row['site'], row['year'], row['crop'])}."
            )

        crop_lookup[(row["site"], row["year"], row["crop"])] = {
            "alberta_township": float(row["alberta_township"]),
            "latitude_degrees": float(
                township_to_latitude_degrees(row["alberta_township"])
            ),
            "season_start_date": season_start_date,
            "kc_initial": float(row["kc_initial"]),
            "kc_mid": float(row["kc_mid"]),
            "kc_end": float(row["kc_end"]),
            "kc_off_season": float(row["kc_off_season"]),
            **parsed_required_dates,
        }

    return crop_lookup


def compute_rooting_depth_cm(growth_stage: str) -> float:
    """Return rooting depth assigned to the current growth stage."""
    stage_to_depth_cm = {
        "off_season": 35.0,
        "seeding_to_emergence": 45.0,
        "emergence_to_anthesis": 80.0,
        "anthesis_peak": 90.0,
        "anthesis_to_early_grain_fill": 90.0,
        "grain_fill_to_maturity": 85.0,
        "maturity_to_harvest": 75.0,
    }
    return stage_to_depth_cm.get(growth_stage, 75.0)


def compute_root_weights(mid_cm: np.ndarray, rooting_depth_cm: float) -> np.ndarray:
    """Distribute root uptake weights across the soil profile."""
    mid_cm = np.asarray(mid_cm, dtype=float)
    rooting_depth_cm = max(float(rooting_depth_cm), 1.0)

    normalized_depth = np.clip(mid_cm / rooting_depth_cm, 0.0, None)
    active_weights = np.exp(-0.55 * normalized_depth)
    active_weights = np.where(mid_cm <= rooting_depth_cm, active_weights, 0.0)
    background_weights = np.exp(-0.03 * np.clip(mid_cm / 90.0, 0.0, None))
    weights = 0.40 * active_weights + 0.60 * background_weights
    minimum_access = np.select(
        [
            mid_cm <= 15.0,
            mid_cm <= 30.0,
            mid_cm <= 45.0,
            mid_cm <= 60.0,
            mid_cm <= 75.0,
            mid_cm <= 90.0,
        ],
        [
            0.10,
            0.13,
            0.16,
            0.20,
            0.24,
            0.28,
        ],
        default=0.10,
    )
    weights = np.maximum(weights, minimum_access)

    if float(np.sum(weights)) <= 1e-12:
        weights = np.where(mid_cm == np.min(mid_cm), 1.0, 0.0)

    return weights / np.sum(weights)
