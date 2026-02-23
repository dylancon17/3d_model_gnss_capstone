import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pyproj import Transformer
import geopandas as gpd
import contextily as ctx
from matplotlib.colors import LogNorm

# ============================================================
# Small utilities
# ============================================================
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def set_pandas_wide() -> None:
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    pd.set_option("display.max_colwidth", None)


def to_numeric_inplace(df: pd.DataFrame, cols) -> None:
    df[cols] = df[cols].apply(pd.to_numeric, errors="coerce")


def ecdf(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.array([]), np.array([])
    xs = np.sort(x)
    ys = np.arange(1, xs.size + 1) / xs.size
    return xs, ys


def p_xlim(x, pct):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (0.0, 1.0)
    return (0.0, np.percentile(x, pct))


def savefig_close(out_dir: str, filename: str, dpi=300) -> None:
    plt.savefig(os.path.join(out_dir, filename), dpi=dpi, bbox_inches="tight")
    plt.close()


def plot_line(out_dir, filename, x, ys, labels, xlabel, ylabel, title, grid=True):
    plt.figure()
    for y, lab in zip(ys, labels):
        plt.plot(x, y, label=lab)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if labels:
        plt.legend()
    if grid:
        plt.grid()
    savefig_close(out_dir, filename)


def plot_single(out_dir, filename, x, y, xlabel, ylabel, title, grid=True):
    plt.figure()
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if grid:
        plt.grid()
    savefig_close(out_dir, filename)


def plot_hist(out_dir, filename, data, bins, xlim, xlabel, ylabel, title, grid=True):
    plt.figure()
    plt.hist(data, bins=bins)
    if xlim is not None:
        plt.xlim(*xlim)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if grid:
        plt.grid()
    savefig_close(out_dir, filename)


# ============================================================
# Loaders
# ============================================================
def load_rtklib_pos(rtklib_file: str) -> pd.DataFrame:
    rtk = pd.read_csv(
        rtklib_file,
        comment="%",
        header=None,
        sep=",",
        engine="python",
    )
    if rtk.shape[1] != 15:
        raise ValueError(f"Unexpected RTKLIB column count: {rtk.shape[1]} (expected 15)")

    rtk.columns = [
        "Week", "GPSTime", "Lat_deg", "Lon_deg", "Height_m",
        "Q", "ns", "sdn", "sde", "sdu",
        "sdne", "sdeu", "sdun", "age", "ratio"
    ]

    num_cols = ["GPSTime", "Lat_deg", "Lon_deg", "Height_m", "Q", "ns"]
    to_numeric_inplace(rtk, num_cols)

    rtk["t"] = rtk["GPSTime"] - rtk["GPSTime"].iloc[0]
    return rtk


def read_sat_rows(rtklib_stats: str):
    sat_rows = []
    with open(rtklib_stats, "r") as f:
        for line in f:
            if line.startswith("$SAT"):
                parts = line.strip().split(",")
                sat_rows.append(parts[1:])
    return sat_rows


def load_rtklib_sat_stats(rtklib_stats: str) -> pd.DataFrame:
    sat_rows = read_sat_rows(rtklib_stats)

    sat_cols = [
        "week", "tow", "sat", "frq", "az", "el", "resp", "resc",
        "vsat", "snr", "fix", "slip", "lock", "outc", "slipc", "rejc",
        "scal", "prob", "idx"
    ]

    sol_stats = pd.DataFrame(sat_rows, columns=sat_cols)

    numeric_cols = [
        "week", "tow", "frq", "az", "el", "resp", "resc",
        "vsat", "snr", "fix", "slip", "lock", "outc", "slipc", "rejc",
        "scal", "prob", "idx"
    ]
    to_numeric_inplace(sol_stats, numeric_cols)

    # Drop unused columns
    sol_stats = sol_stats.drop(
        columns=["vsat", "snr", "fix", "slip", "lock", "outc", "slipc", "rejc", "scal"]
    )

    # Relative time in seconds from first epoch
    sol_stats["t"] = sol_stats["tow"] - sol_stats["tow"].iloc[0]

    sol_stats["resp"] = sol_stats["resp"].abs()
    sol_stats["resc"] = sol_stats["resc"].abs()

    sol_stats = sol_stats[sol_stats["prob"] >= 0.0]

    sol_stats = replace_zero_residuals(sol_stats, col="resp")
    sol_stats = replace_zero_residuals(sol_stats, col="resc")

    return sol_stats


def find_truth_header_line(truth_file: str, startswith="UTCDate") -> int:
    with open(truth_file, "r") as f:
        for i, line in enumerate(f):
            if line.strip().startswith(startswith):
                return i
    return -1


def load_truth_file(truth_file: str) -> pd.DataFrame:
    start_idx = find_truth_header_line(truth_file, "UTCDate")
    if start_idx < 0:
        raise RuntimeError("Could not find CSV header in truth file (line starting with 'UTCDate')")

    truth = pd.read_csv(truth_file, skiprows=start_idx, sep=",", engine="python")

    numeric_cols = truth.columns[1:]
    to_numeric_inplace(truth, numeric_cols)

    truth = truth.dropna(subset=["GPSTime", "Latitude", "Longitude", "H-Ell"])
    truth["t"] = truth["GPSTime"] - truth["GPSTime"].iloc[0]
    return truth


# ============================================================
# Processing
# ============================================================
def replace_zero_residuals(sol_stats: pd.DataFrame, col: str) -> pd.DataFrame:
    sol_stats = sol_stats.copy()
    sat_sys = sol_stats["sat"].astype(str).str[0]

    group_min = (
        sol_stats[col]
        .where(sol_stats[col] != 0)
        .groupby([sol_stats["tow"], sat_sys])
        .transform("min")
    )

    mask = sol_stats[col] == 0
    sol_stats.loc[mask, col] = group_min[mask].fillna(1e-6)
    return sol_stats


def match_truth_to_rtk_time(rtk: pd.DataFrame, truth: pd.DataFrame, tol_s=0.001) -> pd.DataFrame:
    # Use absolute GPSTime as numeric
    rtk = rtk.copy()
    truth = truth.copy()

    rtk["GPSTime_abs"] = pd.to_numeric(rtk["GPSTime"], errors="coerce")
    truth["GPSTime_abs"] = pd.to_numeric(truth["GPSTime"], errors="coerce")

    rtk = rtk.dropna(subset=["GPSTime_abs", "Lat_deg", "Lon_deg", "Height_m"]).copy()
    truth = truth.dropna(subset=["GPSTime_abs", "Latitude", "Longitude", "H-Ell"]).copy()

    rtk_sorted = rtk.sort_values("GPSTime_abs").reset_index(drop=True)
    truth_sorted = truth.sort_values("GPSTime_abs").reset_index(drop=True)

    matched = pd.merge_asof(
        rtk_sorted,
        truth_sorted[["GPSTime_abs", "Latitude", "Longitude", "H-Ell"]],
        on="GPSTime_abs",
        direction="nearest",
        tolerance=tol_s,
    )

    matched = matched.dropna(subset=["Latitude", "Longitude", "H-Ell"]).copy()

    print("\nMatch-only alignment:")
    print(f"  tolerance: {tol_s} s")
    print(f"  matched epochs: {len(matched)} / {len(rtk_sorted)}")

    if len(matched) == 0:
        raise RuntimeError("No matched epochs found. Increase TOL_S or fix time sync between logs.")

    # Redefine RTK stream to matched-only epochs
    rtk = matched.reset_index(drop=True)

    # Define analysis time axis as seconds since first matched epoch
    rtk["t"] = rtk["GPSTime_abs"] - rtk["GPSTime_abs"].iloc[0]
    return rtk


def compute_availability_metrics(rtk: pd.DataFrame):
    # Compute time differences between consecutive epochs
    dt = np.diff(rtk["GPSTime_abs"].values)
    dt = dt[dt > 0]

    typical_dt = np.median(dt)  # seconds per epoch
    data_rate_hz = 1.0 / typical_dt

    logging_duration = rtk["GPSTime_abs"].iloc[-1] - rtk["GPSTime_abs"].iloc[0]
    expected_epochs = int(round(logging_duration / typical_dt)) + 1

    solution_epochs_raw = len(rtk)
    missing_epochs = expected_epochs - solution_epochs_raw
    true_availability = solution_epochs_raw / expected_epochs if expected_epochs > 0 else np.nan

    print("\n==== SOLUTION AVAILABILITY ====")
    print(f"Logging duration (s): {logging_duration:.1f}")
    print(f"Expected epochs     : {expected_epochs}")
    print(f"Solution epochs     : {solution_epochs_raw}")
    print(f"Missing epochs      : {missing_epochs}")
    print(f"Solution availability   : {true_availability:.3%}")

    return typical_dt, data_rate_hz, logging_duration, expected_epochs, solution_epochs_raw, missing_epochs, true_availability


def build_truth_interp_from_matched(rtk: pd.DataFrame) -> pd.DataFrame:
    truth_interp = pd.DataFrame()
    truth_interp["t"] = rtk["t"]
    truth_interp["Lat_deg"] = rtk["Latitude"]
    truth_interp["Lon_deg"] = rtk["Longitude"]
    truth_interp["Height_m"] = rtk["H-Ell"]
    return truth_interp


def llt_to_ecef_arrays(lat_deg, lon_deg, h_m):
    llt_to_ecef = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
    arr = np.array(llt_to_ecef.transform(lon_deg, lat_deg, h_m)).T
    return arr


def ecef_to_enu_errors(truth_interp: pd.DataFrame, rtk: pd.DataFrame):
    truth_ecef = llt_to_ecef_arrays(
        truth_interp["Lat_deg"].values,
        truth_interp["Lon_deg"].values,
        truth_interp["Height_m"].values,
    )

    rtk_ecef = llt_to_ecef_arrays(
        rtk["Lat_deg"].values,
        rtk["Lon_deg"].values,
        rtk["Height_m"].values,
    )

    ecef_err = rtk_ecef - truth_ecef

    ref_lat = truth_interp["Lat_deg"].iloc[0]
    ref_lon = truth_interp["Lon_deg"].iloc[0]

    lat0 = np.deg2rad(ref_lat)
    lon0 = np.deg2rad(ref_lon)

    R = np.array([
        [-np.sin(lon0),               np.cos(lon0),               0],
        [-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0),  np.cos(lat0)],
        [ np.cos(lat0)*np.cos(lon0),  np.cos(lat0)*np.sin(lon0),  np.sin(lat0)],
    ])

    enu_err = (R @ ecef_err.T).T
    return enu_err


def add_error_columns(rtk: pd.DataFrame, enu_err: np.ndarray) -> pd.DataFrame:
    rtk = rtk.copy()
    rtk["E_err"] = enu_err[:, 0]
    rtk["N_err"] = enu_err[:, 1]
    rtk["U_err"] = enu_err[:, 2]
    rtk["Horz_err"] = np.sqrt(rtk["E_err"] ** 2 + rtk["N_err"] ** 2)
    return rtk


def compute_rmse(rtk: pd.DataFrame):
    E = rtk["E_err"].values
    N = rtk["N_err"].values
    U = rtk["U_err"].values

    rmse_h = np.sqrt(np.mean(E**2 + N**2))
    rmse_v = np.sqrt(np.mean(U**2))
    rmse_3d = np.sqrt(np.mean(E**2 + N**2 + U**2))

    print("\n==== RMSE Positioning Errors ====")
    print(f"Horizontal RMSE : {rmse_h:.4f} m")
    print(f"Vertical RMSE   : {rmse_v:.4f} m")
    print(f"3D RMSE         : {rmse_3d:.4f} m")
    return rmse_h, rmse_v, rmse_3d


# ============================================================
# Plotting blocks
# ============================================================
def make_basic_timeseries_plots(out_dir: str, rtk: pd.DataFrame, truth_interp: pd.DataFrame):
    plot_line(
        out_dir,
        "Latitude Over Time.png",
        rtk["t"],
        [rtk["Lat_deg"], truth_interp["Lat_deg"]],
        ["RTK", "Truth"],
        "Time (s)",
        "Latitude (deg)",
        "Latitude Over Time",
    )

    plot_line(
        out_dir,
        "Longitude Over Time.png",
        rtk["t"],
        [rtk["Lon_deg"], truth_interp["Lon_deg"]],
        ["RTK", "Truth"],
        "Time (s)",
        "Longitude (deg)",
        "Longitude Over Time",
    )

    plot_line(
        out_dir,
        "Height Over Time.png",
        rtk["t"],
        [rtk["Height_m"], truth_interp["Height_m"]],
        ["RTK", "Truth"],
        "Time (s)",
        "Height (m)",
        "Height Over Time",
    )

    plot_line(
        out_dir,
        "ENU Error vs Time.png",
        rtk["t"],
        [rtk["N_err"], rtk["E_err"], rtk["U_err"]],
        ["North Err", "East Err", "Up Err"],
        "Time (s)",
        "Error (m)",
        "ENU Error vs Time",
    )

    plot_single(
        out_dir,
        "Horizontal Error Over Time.png",
        rtk["t"],
        rtk["Horz_err"],
        "Time (s)",
        "Horizontal Error (m)",
        "Horizontal Error Over Time",
    )

    plot_single(
        out_dir,
        "Vertical Error Over Time.png",
        rtk["t"],
        rtk["U_err"],
        "Time (s)",
        "Vertical Error (m)",
        "Vertical Error Over Time",
    )

    xmax = np.percentile(rtk["Horz_err"], 50)
    plot_hist(
        out_dir,
        "Horizontal Error Histogram.png",
        rtk["Horz_err"],
        bins=10000,
        xlim=(0, xmax),
        xlabel="Horizontal Error (m)",
        ylabel="Count",
        title="Horizontal Error Histogram",
    )

    xmax = np.percentile(rtk["Horz_err"], 50)
    plot_hist(
        out_dir,
        "Vertical Error Histogram.png",
        rtk["U_err"],
        bins=10000,
        xlim=(0, xmax),
        xlabel="Vertical Error (m)",
        ylabel="Count",
        title="Vertical Error Histogram",
    )

    plot_single(
        out_dir,
        "Satellite Count vs Time.png",
        rtk["t"],
        rtk["ns"],
        "Time (s)",
        "# of Satellites",
        "Satellite Count vs Time",
    )


def make_residual_probability_plots(out_dir: str, sol_stats: pd.DataFrame):
    sol_stats_primary = sol_stats[sol_stats["frq"] == 1]
    sol_stats_secondary = sol_stats[sol_stats["frq"] == 2]

    plt.figure()
    plt.scatter(sol_stats_primary["prob"], sol_stats_primary["resp"].abs(), s=8, alpha=0.6, label="Primary Frequency")
    plt.scatter(sol_stats_secondary["prob"], sol_stats_secondary["resp"].abs(), s=8, alpha=0.6, label="Secondary Frequency")
    plt.xlabel("Probability of Obstruction")
    plt.ylabel("True Double Differenced Pseudorange Error (m)")
    plt.title("Pseudorange Errors at Estimated Probability Levels")
    plt.grid()
    plt.legend()
    savefig_close(out_dir, "Primary_vs_Secondary_Pseudorange_Errors_vs_Obstruction_Probability.png")

    plt.figure()
    plt.scatter(sol_stats_primary["prob"], sol_stats_primary["resc"].abs(), s=8, alpha=0.6, label="Primary Frequency")
    plt.scatter(sol_stats_secondary["prob"], sol_stats_secondary["resc"].abs(), s=8, alpha=0.6, label="Secondary Frequency")
    plt.xlabel("Probability of Obstruction")
    plt.ylabel("True Double Differenced Carrier Phase Error (m)")
    plt.title("Carrier Phase Errors at Estimated Probability Levels")
    plt.grid()
    plt.legend()
    savefig_close(out_dir, "Primary_vs_Secondary_Carrier_Phase_Errors_vs_Obstruction_Probability.png")

    # 2D heatmap with text counts
    plt.figure()
    counts, xedges, yedges, im = plt.hist2d(
        sol_stats["prob"],
        sol_stats["resp"].abs(),
        bins=[20, int(300 / 10)],
        range=[[0.0, 1.0], [0.0, 300.0]],
        norm=LogNorm(),
    )
    plt.colorbar(label="Count (log scale)")

    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])

    for i, x in enumerate(xcenters):
        for j, y in enumerate(ycenters):
            count = counts[i, j]
            if count > 0:
                plt.text(x, y, f"{int(count)}", color="white", ha="center", va="center", fontsize=5)

    plt.xlabel("Probability of Obstruction")
    plt.ylabel("True Double Differenced Pseudorange Error (m)")
    plt.title("Primary Pseudorange Errors at Estimated Probability Levels")
    plt.grid()
    savefig_close(out_dir, "Pseudorange Errors vs Obstruction Probability (Heatmap).png")

    plt.figure()
    plt.hist(sol_stats["prob"], bins=20, range=(0.0, 1.0))
    plt.xlabel("Probability of Obstruction")
    plt.ylabel("Count")
    plt.title("Distribution of Estimated Obstruction Probability")
    plt.grid()
    savefig_close(out_dir, "Obstruction Probability Histogram.png")


def make_skyplots(out_dir: str, sol_stats: pd.DataFrame):
    sol_stats_primary = sol_stats[sol_stats["frq"] == 1]

    theta = np.deg2rad(sol_stats_primary["az"])
    r = 90 - sol_stats_primary["el"]

    tn = (sol_stats_primary["prob"] < 0.5) & (sol_stats_primary["resp"] < 10)
    tp = (sol_stats_primary["prob"] >= 0.5) & (sol_stats_primary["resp"] >= 10)
    fp = (sol_stats_primary["prob"] >= 0.5) & (sol_stats_primary["resp"] < 10)
    fn = (sol_stats_primary["prob"] < 0.5) & (sol_stats_primary["resp"] >= 10)

    fig = plt.figure(figsize=(7, 7))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlim(0, 90)

    ax.scatter(theta[tn], r[tn], s=12, label="True Negative")
    ax.scatter(theta[tp], r[tp], s=12, label="True Positive")
    ax.scatter(theta[fp], r[fp], s=12, label="False Positive")
    ax.scatter(theta[fn], r[fn], s=12, label="False Negative")

    ax.set_rgrids([0, 30, 60, 90], labels=["0", "30", "60°", "90°"])
    ax.set_title("Skyplot: Probability vs Residual Classification", pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))

    plt.savefig(os.path.join(out_dir, "Pseudorange Residuals.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # Skyplot by satellite
    theta = np.deg2rad(sol_stats_primary["az"])
    r = 90 - sol_stats_primary["el"]

    fig = plt.figure(figsize=(7, 7))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlim(0, 90)

    for sat in sorted(sol_stats_primary["sat"].unique()):
        mask = sol_stats_primary["sat"] == sat
        ax.scatter(theta[mask], r[mask], s=12, label=sat)

    ax.set_rgrids([0, 30, 60, 90], labels=["90°", "60°", "30°", "0°"])
    ax.set_title("Skyplot by Satellite", pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), title="Satellite", fontsize=8)

    plt.savefig(os.path.join(out_dir, "Skyplot_by_Satellite.png"), dpi=300, bbox_inches="tight")
    plt.close()


def make_trajectory_map(out_dir: str, rtk: pd.DataFrame, truth_interp: pd.DataFrame):
    gdf_rtk = gpd.GeoDataFrame(
        rtk,
        geometry=gpd.points_from_xy(rtk["Lon_deg"], rtk["Lat_deg"]),
        crs="EPSG:4326",
    )

    gdf_truth = gpd.GeoDataFrame(
        truth_interp,
        geometry=gpd.points_from_xy(truth_interp["Lon_deg"], truth_interp["Lat_deg"]),
        crs="EPSG:4326",
    )

    gdf_rtk = gdf_rtk.to_crs(epsg=3857)
    gdf_truth = gdf_truth.to_crs(epsg=3857)

    fig, ax = plt.subplots(figsize=(10, 10))

    gdf_truth.plot(ax=ax, color="green", linewidth=5, label="Truth")
    gdf_rtk.plot(ax=ax, color="red", linewidth=1, alpha=0.8, label="Solution")

    ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)

    ax.set_title("Solution vs Truth (Map View)")
    ax.set_axis_off()
    ax.legend()

    plt.savefig(os.path.join(out_dir, "Map View - Solution vs Truth.png"), dpi=300, bbox_inches="tight")
    plt.close()


# ============================================================
# CDF / percentile outputs
# ============================================================
def percentiles_by_q_and_write_csv(out_dir: str, rtk: pd.DataFrame, expected_epochs: int, solution_epochs_raw: int, solution_rate: float):
    print("\nSolution quality (Q) breakdown:")
    q_counts = rtk["Q"].value_counts(dropna=False).sort_index()
    print(q_counts)

    q_rows = []
    for q in sorted(rtk["Q"].dropna().unique()):
        vals = rtk.loc[rtk["Q"] == q, "Horz_err"].astype(float).values
        vals = vals[np.isfinite(vals)]
        if vals.size < 10:
            continue

        q_rows.append({
            "Q": int(q),
            "n": int(vals.size),
            "P50": float(np.nanpercentile(vals, 50)),
            "P95": float(np.nanpercentile(vals, 95)),
            "P99": float(np.nanpercentile(vals, 99)),
            "Max": float(np.nanmax(vals)),
            "epochs_with_solution": int(solution_epochs_raw),
            "total_epochs": int(expected_epochs),
            "solution_rate": float(solution_rate),
        })

    q_table = pd.DataFrame(q_rows)
    q_table.to_csv(os.path.join(out_dir, "HorzError_percentiles_by_Q.csv"), index=False)

    q_table_print = q_table.sort_values("Q").copy()
    if not q_table_print.empty:
        q_table_print[["P50", "P95", "P99", "Max", "solution_rate"]] = q_table_print[
            ["P50", "P95", "P99", "Max", "solution_rate"]
        ].round(4)

    print("\n==== Horz_err percentiles by Q ====")
    print(q_table_print.to_string(index=False))

    return q_table


def make_cdf_plots(out_dir: str, rtk: pd.DataFrame, cdf_pct=99):
    E_abs = np.abs(rtk["E_err"].astype(float).values)
    N_abs = np.abs(rtk["N_err"].astype(float).values)
    U_abs = np.abs(rtk["U_err"].astype(float).values)
    H = np.abs(rtk["Horz_err"].astype(float).values)

    plt.figure()
    for data, lab in zip([E_abs, N_abs, U_abs], ["|E|", "|N|", "|U|"]):
        xs, ys = ecdf(data)
        if xs.size:
            plt.plot(xs, ys, label=lab)
    plt.xlim(p_xlim(np.r_[E_abs, N_abs, U_abs], cdf_pct))
    plt.xlabel("Absolute error (m)")
    plt.ylabel("CDF")
    plt.title(f"CDF Absolute ENU Component Errors (0–P{cdf_pct})")
    plt.grid(True)
    plt.legend()
    savefig_close(out_dir, f"CDF_ENU_abs_P{cdf_pct}.png")

    plt.figure()
    xs, ys = ecdf(H)
    plt.plot(xs, ys, label="Horizontal |EN|")
    plt.xlim(p_xlim(H, cdf_pct))
    plt.xlabel("Horizontal error (m)")
    plt.ylabel("CDF")
    plt.title(f"CDF Horizontal Error (0–P{cdf_pct})")
    plt.grid(True)
    plt.legend()
    savefig_close(out_dir, f"CDF_Horizontal_P{cdf_pct}.png")

    plt.figure()
    for q in sorted(rtk["Q"].dropna().unique()):
        vals = rtk.loc[rtk["Q"] == q, "Horz_err"].values
        vals = vals[np.isfinite(vals)]
        if vals.size < 20:
            continue
        xs, ys = ecdf(vals)
        plt.plot(xs, ys, label=f"Q={int(q)} (n={vals.size})")
    plt.xlim(p_xlim(H, 99))
    plt.xlabel("Horizontal error (m)")
    plt.ylabel("CDF")
    plt.title("CDF Horizontal Error by Solution Quality (Q)")
    plt.grid(True)
    plt.legend()
    savefig_close(out_dir, "CDF_Horizontal_by_Q_P99.png")


# ============================================================
# Solution availability over time
# ============================================================
def add_has_solution_column(rtk: pd.DataFrame) -> pd.DataFrame:
    rtk = rtk.copy()
    has_solution = (
        rtk["Q"].notna() &
        (rtk["Q"] > 0) &
        np.isfinite(rtk["Lat_deg"]) &
        np.isfinite(rtk["Lon_deg"]) &
        np.isfinite(rtk["Height_m"])
    )
    rtk["has_solution"] = has_solution.astype(int)
    return rtk


def plot_solution_availability(out_dir: str, rtk: pd.DataFrame, expected_epochs: int, typical_dt: float, solution_rate: float):
    # Build full expected timeline
    t_full = np.arange(0, expected_epochs * typical_dt, typical_dt)
    has_solution_full = np.zeros_like(t_full, dtype=int)

    indices = np.round(rtk["t"].values / typical_dt).astype(int)
    indices = indices[indices < len(t_full)]
    has_solution_full[indices] = rtk["has_solution"].values

    plt.figure(figsize=(12, 4))
    plt.step(t_full, has_solution_full, where="post", color="blue", linewidth=1.5, label="Solution Available")
    plt.fill_between(t_full, 0, 1, where=(has_solution_full == 0), color="red", alpha=0.2, label="No Solution")

    plt.ylim(-0.1, 1.1)
    plt.xlabel("Time (s)")
    plt.ylabel("Solution Availability")
    plt.title(f"Solution Availability Over Time ({solution_rate:.2%} available)")
    plt.grid(True)
    plt.legend(loc="upper right")
    plt.tight_layout()

    savefig_close(out_dir, "Solution Availability vs Time.png")


# ============================================================
# Main
# ============================================================
def parse_args(argv):
    if len(argv) != 6:
        print("\nUsage:")
        print("  python plots.py <rtklib_output.pos> <rtklib_output.pos.stat> <truth_output.txt> <truth_residuals.pos.stat><output_dir>\n")
        sys.exit(1)

    return argv[1], argv[2], argv[3], argv[4], argv[5]


def main(argv):
    rtklib_file, rtklib_stats, truth_file, truth_stats, out_dir = parse_args(argv)

    ensure_dir(out_dir)
    set_pandas_wide()

    print(f"\nSolution: {rtklib_file}")
    print(f"Truth file : {truth_file}\n")

    # Load data
    rtk = load_rtklib_pos(rtklib_file)
    sol_stats = load_rtklib_sat_stats(rtklib_stats)
    truth = load_truth_file(truth_file)

    # Match truth to RTK time
    TOL_S = 0.001
    rtk = match_truth_to_rtk_time(rtk, truth, tol_s=TOL_S)

    # Availability metrics
    typical_dt, data_rate_hz, logging_duration, expected_epochs, solution_epochs_raw, missing_epochs, true_availability = (
        compute_availability_metrics(rtk)
    )

    # Build truth_interp for downstream compatibility
    truth_interp = build_truth_interp_from_matched(rtk)

    # ENU errors
    enu_err = ecef_to_enu_errors(truth_interp, rtk)
    rtk = add_error_columns(rtk, enu_err)

    # RMSE
    compute_rmse(rtk)

    # Plots
    make_basic_timeseries_plots(out_dir, rtk, truth_interp)
    make_residual_probability_plots(out_dir, sol_stats)
    make_skyplots(out_dir, sol_stats)
    make_trajectory_map(out_dir, rtk, truth_interp)

    # CDF metrics and plots
    percentiles_by_q_and_write_csv(out_dir, rtk, expected_epochs, solution_epochs_raw, true_availability)
    make_cdf_plots(out_dir, rtk, cdf_pct=99)

    # Solution availability over time
    rtk = add_has_solution_column(rtk)
    plot_solution_availability(out_dir, rtk, expected_epochs, typical_dt, true_availability)

    # plt.show()

if __name__ == "__main__":
    main(sys.argv)
