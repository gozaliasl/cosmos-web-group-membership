"""
Optimise membership search radius per group.

Strategies:
- xray: Use R500 or R200 from X-ray catalog when available (with optional scaling/clamp).
- redshift: Redshift-dependent default (typical R200 scale for group-mass halos).
- fixed: Single fixed radius (no optimisation).
- hybrid: X-ray when available, else redshift default, else fixed.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Union

# Project root (membership_determination/radius_optimizer.py -> parent.parent)
BASE_DIR = Path(__file__).resolve().parent.parent


def _load_xray_r200_r500(
    xray_path: Path,
    catalog_name: str,
) -> pd.DataFrame:
    """Load Group_ID, R200_kpc, R500_kpc for a catalog from X-ray results."""
    if not xray_path.exists():
        return pd.DataFrame()
    try:
        if xray_path.suffix.lower() == ".csv":
            df = pd.read_csv(xray_path)
        else:
            from astropy.table import Table
            t = Table.read(xray_path)
            df = t.to_pandas()
    except Exception:
        return pd.DataFrame()
    if "Group_ID" not in df.columns:
        return pd.DataFrame()
    if "Catalog_Name" in df.columns:
        df = df[df["Catalog_Name"].astype(str).str.strip().str.upper() == catalog_name.upper()].copy()
    if df.empty:
        return df
    out = df[["Group_ID"]].copy()
    for key in ["R200_kpc", "R500_kpc", "R200_Luminosity_kpc", "R500_Luminosity_kpc"]:
        if key in df.columns:
            out[key] = pd.to_numeric(df[key], errors="coerce")
    return out


def get_radius_from_xray(
    group_id: Union[int, str],
    xray_df: pd.DataFrame,
    use_r200: bool = True,
    use_r500_fallback: bool = True,
    scale: float = 1.0,
    r_min_kpc: float = 300.0,
    r_max_kpc: float = 1500.0,
) -> Optional[float]:
    """
    Get optimal radius for one group from X-ray R200/R500.
    Returns None if not found or invalid.
    When use_r200=True and use_r500_fallback=False, only R200 is used (no R500 fallback).
    """
    if xray_df is None or xray_df.empty:
        return None
    row = xray_df[xray_df["Group_ID"].astype(str) == str(group_id)]
    if row.empty:
        return None
    row = row.iloc[0]
    r = None
    if use_r200:
        for col in ["R200_kpc", "R200_Luminosity_kpc"]:
            if col in row and np.isfinite(row[col]) and row[col] > 0:
                r = float(row[col])
                break
    if r is None and use_r500_fallback:
        for col in ["R500_kpc", "R500_Luminosity_kpc"]:
            if col in row and np.isfinite(row[col]) and row[col] > 0:
                r = float(row[col])
                break
    if r is None:
        return None
    r = r * scale
    r = max(r_min_kpc, min(r_max_kpc, r))
    return r


def redshift_default_radius(
    z: float,
    r_min_kpc: float = 300.0,
    r_max_kpc: float = 1200.0,
) -> float:
    """
    Redshift-dependent default radius (typical R200 scale for ~10^13 Msun at z).
    Simple scaling: R ~ 400 + 80*z kpc (approximate).
    """
    r = 400.0 + 80.0 * z
    return float(max(r_min_kpc, min(r_max_kpc, r)))


def get_per_group_radius(
    groups: pd.DataFrame,
    catalog_name: str,
    id_col: str,
    z_col: str,
    radius_mode: str = "fixed",
    fixed_radius_kpc: float = 500.0,
    xray_path: Optional[Path] = None,
    xray_scale: float = 1.0,
    xray_r_min: float = 300.0,
    xray_r_max: float = 1500.0,
    use_r200: bool = True,
    use_r500_fallback: bool = True,
) -> Tuple[np.ndarray, pd.Series]:
    """
    Compute per-group membership radius.

    Parameters
    ----------
    groups : DataFrame
        Group catalog with id and redshift columns.
    catalog_name : str
        'CW-All' or 'CW-HCG'.
    id_col : str
        Group ID column name.
    z_col : str
        Redshift column name.
    radius_mode : str
        'fixed' | 'xray' | 'redshift' | 'hybrid'
    fixed_radius_kpc : float
        Used for 'fixed' and as fallback in 'hybrid'.
    xray_path : Path, optional
        Path to X-ray catalog (FITS/CSV). If None, uses default cw_all path for CW-All.
    xray_scale : float
        Multiply X-ray R200/R500 by this (e.g. 1.2 to include outskirts).
    xray_r_min, xray_r_max : float
        Clamp X-ray-derived radius to [r_min, r_max].
    use_r200 : bool
        Prefer R200 over R500 when both available.
    use_r500_fallback : bool
        If False, use only R200 (no R500 fallback) when radius_mode is xray/hybrid.

    Returns
    -------
    radius_arr : np.ndarray
        One radius per row in groups (same order).
    radius_series : pd.Series
        Index = group index, value = radius (for logging).
    """
    n = len(groups)
    radius_arr = np.full(n, fixed_radius_kpc, dtype=float)
    radius_series = pd.Series(index=groups.index, data=fixed_radius_kpc, dtype=float)

    xray_df = None
    if radius_mode in ("xray", "hybrid"):
        path = xray_path
        if path is None:
            path = BASE_DIR / "outputs" / "results" / "cw_all" / "xray_catalog.fits"
            if not path.exists():
                path = BASE_DIR / "outputs" / "results" / "cw_all" / "xray_catalog.csv"
        path = Path(path) if path is not None else None
        if path is not None and path.exists():
            if path.suffix.lower() == ".csv":
                xray_df = pd.read_csv(path)
            else:
                xray_df = _load_xray_r200_r500(path, catalog_name)
            if xray_df.empty and path.suffix.lower() == ".fits":
                from astropy.table import Table
                t = Table.read(path)
                xray_df = t.to_pandas()
                if "Catalog_Name" in xray_df.columns:
                    xray_df = xray_df[
                        xray_df["Catalog_Name"].astype(str).str.strip().str.upper()
                        == catalog_name.upper()
                    ]

    for i, (idx, row) in enumerate(groups.iterrows()):
        gid = row[id_col]
        z = row[z_col]
        if not np.isfinite(z) or z < 0:
            radius_arr[i] = fixed_radius_kpc
            continue
        r = None
        if radius_mode in ("xray", "hybrid") and xray_df is not None and not xray_df.empty:
            r = get_radius_from_xray(
                gid,
                xray_df,
                use_r200=use_r200,
                use_r500_fallback=use_r500_fallback,
                scale=xray_scale,
                r_min_kpc=xray_r_min,
                r_max_kpc=xray_r_max,
            )
        if r is None and radius_mode == "redshift":
            r = redshift_default_radius(z)
        if r is None and radius_mode == "hybrid":
            r = redshift_default_radius(z)
        if r is None:
            r = fixed_radius_kpc
        radius_arr[i] = r
        radius_series.loc[idx] = r

    return radius_arr, radius_series
