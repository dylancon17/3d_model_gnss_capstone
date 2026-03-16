import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pyproj import Transformer
import geopandas as gpd
import contextily as ctx
from matplotlib.colors import LogNorm

import traceback

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

# ============================
# COMMAND LINE ARGUMENTS
# ============================
if len(sys.argv) != 6:
    print("\nUsage:")
    print("  python plots.py <rtklib_output.pos> <rtklib_output.pos.stat> <truth_output.txt> <truth_residuals.pos.stat><output_dir>\n")
    sys.exit(1)

rtklib_file = sys.argv[1]
rtklib_stats = sys.argv[2]
truth_file = sys.argv[3]
truth_stats = sys.argv[4]
out_dir = sys.argv[5]

os.makedirs(out_dir, exist_ok=True)

print(f"\nSolution: {rtklib_file}")
print(f"Truth file : {truth_file}\n")

# ============================
# LOAD RTKLIB OUTPUT (.pos)
# ============================
rtk = pd.read_csv(
    rtklib_file,
    comment='%',
    header=None,
    sep=",",
    engine="python"
)

if rtk.shape[1] != 15:
    raise ValueError(f"Unexpected RTKLIB column count: {rtk.shape[1]} (expected 15)")

rtk.columns = [
    "Week", "GPSTime", "Lat_deg", "Lon_deg", "Height_m",
    "Q", "ns", "sdn", "sde", "sdu",
    "sdne", "sdeu", "sdun", "age", "ratio"
]

# Make sure key numeric columns are numeric (helps merge_asof)
for c in ["GPSTime", "Lat_deg", "Lon_deg", "Height_m", "Q", "ns"]:
    rtk[c] = pd.to_numeric(rtk[c], errors="coerce")

# Relative time (seconds from start)
rtk["t"] = rtk["GPSTime"] - rtk["GPSTime"].iloc[0]

# Time-based solution availability
start_time = rtk["GPSTime"].iloc[0]
end_time = rtk["GPSTime"].iloc[-1]
duration = end_time - start_time

epochs = len(rtk)

if duration > 0:
    solution_rate_time = epochs / duration
else:
    solution_rate_time = float("nan")

print(f"\nTime-based solution availability: {solution_rate_time:.3f} epochs/sec")


# ============================
# LOAD RTKLIB $SAT STATS FILE
# ============================
sat_rows = []

with open(rtklib_stats, "r") as f:
    for line in f:
        if not line.startswith("$SAT"):
            continue

        parts = line.strip().split(",")

        sat_rows.append(parts[1:])  # drop "$SAT"

sat_cols = [
    "week",   # GPS week
    "tow",    # GPS time of week (s)
    "sat",    # satellite ID (e.g., G06, R21)
    "frq",    # frequency index
    "az",     # azimuth (deg)
    "el",     # elevation (deg)
    "resp",   # pseudorange residual (m)
    "resc",   # carrier-phase residual (m)
    "vsat",   # valid data flag
    "snr",    # signal strength (dB-Hz)
    "fix",    # ambiguity flag
    "slip",   # cycle slip flag
    "lock",   # carrier lock count
    "outc",   # outage count
    "slipc",  # cycle slip count
    "rejc",   # reject count
    "scal",   # observation weight scaling
    "prob",   # obstruction probability
    "idx",     # Index
    "dHeight",
    "dprob", #
    "dsearchmax"
]

try:
    sol_stats = pd.DataFrame(sat_rows, columns=sat_cols)

    numeric_cols = [
        "week", "tow", "frq",
        "az", "el",
        "resp", "resc",
        "vsat", "snr", "fix", "slip",
        "lock", "outc", "slipc", "rejc",
        "scal", "prob", "idx", "dHeight", "dprob", "dsearchmax"
    ]

    sol_stats[numeric_cols] = sol_stats[numeric_cols].apply(
        pd.to_numeric, errors="coerce"
    )

    sol_stats = sol_stats.drop(
        columns=["vsat", "snr", "fix", "slip", "lock", "outc", "slipc", "rejc", "scal"]
    )

    # Relative time in seconds from first epoch
    sol_stats["t"] = sol_stats["tow"] - sol_stats["tow"].iloc[0]

    sol_stats["resp"] = abs(sol_stats["resp"])
    sol_stats["resc"] = abs(sol_stats["resc"])

    # Remove -1-probability entries (indicates always out of bounds)
    sol_stats = sol_stats[sol_stats["prob"] >= 0.0]

    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    pd.set_option("display.max_colwidth", None)

    # Assign the base the minimum residual across the other satellites, as it can't be better than that
    sat_sys = sol_stats["sat"].astype(str).str[0]

    # Compute per-(tow, sat_sys) minimum non-zero resp
    group_min = (
        sol_stats["resp"]
        .where(sol_stats["resp"] != 0)
        .groupby([sol_stats["tow"], sat_sys])
        .transform("min")
    )

    # Replace zeros with group minimum (fallback epsilon if group_min is NaN)
    mask = sol_stats["resp"] == 0
    sol_stats.loc[mask, "resp"] = group_min[mask].fillna(1e-6)

    # Compute per-(tow, sat_sys) minimum non-zero resc
    group_min = (
        sol_stats["resc"]
        .where(sol_stats["resc"] != 0)
        .groupby([sol_stats["tow"], sat_sys])
        .transform("min")
    )

    # Replace zeros with group minimum (fallback epsilon if group_min is NaN)
    mask = sol_stats["resc"] == 0
    sol_stats.loc[mask, "resc"] = group_min[mask].fillna(1e-6)
except:
    traceback.print_exc()


# ============================
# LOAD TRUTH FILE (NovAtel-style CSV)
# ============================
with open(truth_file, "r") as f:
    lines = f.readlines()

start_idx = None
for i, line in enumerate(lines):
    if line.strip().startswith("UTCDate"):
        start_idx = i
        break

if start_idx is None:
    raise RuntimeError("Could not find CSV header in truth file (line starting with 'UTCDate')")

truth = pd.read_csv(
    truth_file,
    skiprows=start_idx,
    sep=",",
    engine="python"
)

numeric_cols = truth.columns[1:]
truth[numeric_cols] = truth[numeric_cols].apply(pd.to_numeric, errors="coerce")

truth = truth.dropna(subset=["GPSTime", "Latitude", "Longitude", "H-Ell"])

truth["t"] = truth["GPSTime"] - truth["GPSTime"].iloc[0]

# ============================
# MATCH TRUTH TO RTK TIME
# ============================

# Use absolute GPSTime as numeric
rtk["GPSTime_abs"] = pd.to_numeric(rtk["GPSTime"], errors="coerce")
truth["GPSTime_abs"] = pd.to_numeric(truth["GPSTime"], errors="coerce")

# Drop rows missing the needed fields
rtk = rtk.dropna(subset=["GPSTime_abs", "Lat_deg", "Lon_deg", "Height_m"]).copy()
truth = truth.dropna(subset=["GPSTime_abs", "Latitude", "Longitude", "H-Ell"]).copy()

# Sort for merge_asof
rtk_sorted = rtk.sort_values("GPSTime_abs").reset_index(drop=True)
truth_sorted = truth.sort_values("GPSTime_abs").reset_index(drop=True)

# Tolerance for "aligned" epochs (seconds)
TOL_S = 0.001

matched = pd.merge_asof(
    rtk_sorted,
    truth_sorted[["GPSTime_abs", "Latitude", "Longitude", "H-Ell"]],
    on="GPSTime_abs",
    direction="nearest",
    tolerance=TOL_S
)

# Drop RTK epochs that had no truth match within tolerance
matched = matched.dropna(subset=["Latitude", "Longitude", "H-Ell"]).copy()

print("\nMatch-only alignment:")
print(f"  tolerance: {TOL_S} s")
print(f"  matched epochs: {len(matched)} / {len(rtk_sorted)}")

if len(matched) == 0:
    raise RuntimeError("No matched epochs found. Increase TOL_S or fix time sync between logs.")

# Redefine RTK stream to matched-only epochs (so everything below uses matched epochs)
rtk = matched.reset_index(drop=True)

# Define analysis time axis as seconds since first matched epoch
rtk["t"] = rtk["GPSTime_abs"] - rtk["GPSTime_abs"].iloc[0]

# Build truth_interp with the same shape as rtk (keeps your downstream code unchanged)
truth_interp = pd.DataFrame()
truth_interp["t"] = rtk["t"]
truth_interp["Lat_deg"] = rtk["Latitude"]
truth_interp["Lon_deg"] = rtk["Longitude"]
truth_interp["Height_m"] = rtk["H-Ell"]


# ============================
# COORDINATE TRANSFORM SETUP
# ============================
llt_to_ecef = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)

truth_ecef = np.array(llt_to_ecef.transform(
    truth_interp["Lon_deg"].values,
    truth_interp["Lat_deg"].values,
    truth_interp["Height_m"].values
)).T

rtk_ecef = np.array(llt_to_ecef.transform(
    rtk["Lon_deg"].values,
    rtk["Lat_deg"].values,
    rtk["Height_m"].values
)).T

ecef_err = rtk_ecef - truth_ecef

# ============================
# ENU ERROR COMPUTATION
# ============================
ref_lat = truth_interp["Lat_deg"].iloc[0]
ref_lon = truth_interp["Lon_deg"].iloc[0]
ref_h   = truth_interp["Height_m"].iloc[0]

lat0 = np.deg2rad(ref_lat)
lon0 = np.deg2rad(ref_lon)

R = np.array([
    [-np.sin(lon0),              np.cos(lon0),              0],
    [-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0), np.cos(lat0)],
    [ np.cos(lat0)*np.cos(lon0),  np.cos(lat0)*np.sin(lon0), np.sin(lat0)]
])

enu_err = (R @ ecef_err.T).T

rtk["E_err"] = enu_err[:, 0]
rtk["N_err"] = enu_err[:, 1]
rtk["U_err"] = enu_err[:, 2]
rtk["Horz_err"] = np.sqrt(rtk["E_err"]**2 + rtk["N_err"]**2)

# ============================
# RMSE COMPUTATION
# ============================
E = rtk["E_err"].values
N = rtk["N_err"].values
U = rtk["U_err"].values

rmse_h = np.sqrt(np.mean(E**2 + N**2))
rmse_v = np.sqrt(np.mean(U**2))
rmse_3d = np.sqrt(np.mean(E**2 + N**2 + U**2))

rmse_text = (
    "\n==== RMSE Positioning Errors ====\n"
    f"Horizontal RMSE : {rmse_h:.4f} m\n"
    f"Vertical RMSE   : {rmse_v:.4f} m\n"
    f"3D RMSE         : {rmse_3d:.4f} m\n"
)

# Print to console
print(rmse_text, end="")

# Save to file
with open(os.path.join(out_dir, "RMSE_Positioning_Errors.txt"), "w") as f:
    f.write(rmse_text)

# ============================
# PLOTTING
# ============================
plt.figure()
plt.plot(rtk["t"], rtk["Lat_deg"], label="RTK")
plt.plot(rtk["t"], truth_interp["Lat_deg"], label="Truth")
plt.xlabel("Time (s)")
plt.ylabel("Latitude (deg)")
plt.title("Latitude Over Time")
plt.legend()
plt.grid()
plt.savefig(os.path.join(out_dir, "Latitude Over Time.png"), dpi=300, bbox_inches="tight")
plt.close()

plt.figure()
plt.plot(rtk["t"], rtk["Lon_deg"], label="RTK")
plt.plot(rtk["t"], truth_interp["Lon_deg"], label="Truth")
plt.xlabel("Time (s)")
plt.ylabel("Longitude (deg)")
plt.title("Longitude Over Time")
plt.legend()
plt.grid()
plt.savefig(os.path.join(out_dir, "Longitude Over Time.png"), dpi=300, bbox_inches="tight")
plt.close()

plt.figure()
plt.plot(rtk["t"], rtk["Height_m"], label="RTK")
plt.plot(rtk["t"], truth_interp["Height_m"], label="Truth")
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.title("Height Over Time")
plt.legend()
plt.grid()
plt.savefig(os.path.join(out_dir, "Height Over Time.png"), dpi=300, bbox_inches="tight")
plt.close()

plt.figure()
plt.plot(rtk["t"], rtk["N_err"], label="North Err")
plt.plot(rtk["t"], rtk["E_err"], label="East Err")
plt.plot(rtk["t"], rtk["U_err"], label="Up Err")
plt.xlabel("Time (s)")
plt.ylabel("Error (m)")
plt.title("ENU Error vs Time")
plt.legend()
plt.grid()
plt.savefig(os.path.join(out_dir, "ENU Error vs Time.png"), dpi=300, bbox_inches="tight")
plt.close()

plt.figure()
plt.plot(rtk["t"], rtk["Horz_err"])
plt.xlabel("Time (s)")
plt.ylabel("Horizontal Error (m)")
plt.title("Horizontal Error Over Time")
plt.grid()
plt.savefig(os.path.join(out_dir, "Horizontal Error Over Time.png"), dpi=300, bbox_inches="tight")
plt.close()

plt.figure()
plt.plot(rtk["t"], rtk["U_err"])
plt.xlabel("Time (s)")
plt.ylabel("Vertical Error (m)")
plt.title("Vertical Error Over Time")
plt.grid()
plt.savefig(os.path.join(out_dir, "Vertical Error Over Time.png"), dpi=300, bbox_inches="tight")
plt.close()

xmax = np.percentile(rtk["Horz_err"], 99)
plt.figure()
plt.hist(rtk["Horz_err"], bins=10000)
plt.xlim(0, xmax)
plt.xlabel("Horizontal Error (m)")
plt.ylabel("Count")
plt.title("Horizontal Error Histogram")
plt.grid()
plt.savefig(os.path.join(out_dir, "Horizontal Error Histogram.png"), dpi=300, bbox_inches="tight")
plt.close()

xmax = np.percentile(rtk["Horz_err"], 99)
plt.figure()
plt.hist(rtk["U_err"], bins=10000)
plt.xlim(0, xmax)
plt.xlabel("Vertical Error (m)")
plt.ylabel("Count")
plt.title("Vertical Error Histogram")
plt.grid()
plt.savefig(os.path.join(out_dir, "Vertical Error Histogram.png"), dpi=300, bbox_inches="tight")
plt.close()

plt.figure()
plt.plot(rtk["t"], rtk["ns"])
plt.xlabel("Time (s)")
plt.ylabel("# of Satellites")
plt.title("Satellite Count vs Time")
plt.grid()
plt.savefig(os.path.join(out_dir, "Satellite Count vs Time.png"), dpi=300, bbox_inches="tight")
plt.close()

try:
    sol_stats_primary = sol_stats[sol_stats["frq"] == 1]
    sol_stats_secondary = sol_stats[sol_stats["frq"] == 2]

    plt.figure()
    plt.scatter(sol_stats_primary["prob"], abs(sol_stats_primary["resp"]), s=8, alpha=0.6, label="Primary Frequency")
    plt.scatter(sol_stats_secondary["prob"], abs(sol_stats_secondary["resp"]), s=8, alpha=0.6, label="Secondary Frequency")
    plt.xlabel("Probability of Obstruction")
    plt.ylabel("True Double Differenced Pseudorange Error (m)")
    plt.title("Pseudorange Errors at Estimated Probability Levels")
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(out_dir, "Primary_vs_Secondary_Pseudorange_Errors_vs_Obstruction_Probability.png"), dpi=300, bbox_inches="tight")
    plt.close()

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)

    plt.figure()
    plt.scatter(sol_stats_primary["prob"], abs(sol_stats_primary["resc"]), s=8, alpha=0.6, label="Primary Frequency")
    plt.scatter(sol_stats_secondary["prob"], abs(sol_stats_secondary["resc"]), s=8, alpha=0.6, label="Secondary Frequency")
    plt.xlabel("Probability of Obstruction")
    plt.ylabel("True Double Differenced Carrier Phase Error (m)")
    plt.title("Carrier Phase Errors at Estimated Probability Levels")
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(out_dir, "Primary_vs_Secondary_Carrier_Phase_Errors_vs_Obstruction_Probability.png"), dpi=300, bbox_inches="tight")
    plt.close()
except:
    traceback.print_exc()

try:
    if sol_stats.empty:
        pass
    else:
        plt.figure()
        
        counts, xedges, yedges, im = plt.hist2d(
            sol_stats["prob"],
            abs(sol_stats["resp"]),
            bins=[20, int(300 / 10)],
            range=[[0.0, 1.0], [0.0, 300.0]],
            norm=LogNorm()
        )

        plt.colorbar(label="Count (log scale)")

        xcenters = 0.5 * (xedges[:-1] + xedges[1:])
        ycenters = 0.5 * (yedges[:-1] + yedges[1:])

        for i, x in enumerate(xcenters):
            for j, y in enumerate(ycenters):
                count = counts[i, j]
                if count > 0:
                    plt.text(
                        x, y,
                        f"{int(count)}",
                        color="white",
                        ha="center",
                        va="center",
                        fontsize=5
                    )

        plt.xlabel("Probability of Obstruction")
        plt.ylabel("True Double Differenced Pseudorange Error (m)")
        plt.title("Primary Pseudorange Errors at Estimated Probability Levels")
        plt.grid()
        plt.savefig(os.path.join(out_dir, "Pseudorange Errors vs Obstruction Probability (Heatmap).png"), dpi=300, bbox_inches="tight")
        plt.close()
except:
    traceback.print_exc()

try:
    plt.figure()
    if sol_stats_primary.empty:
        pass
    else:
        counts, xedges, yedges, im = plt.hist2d(
            sol_stats_primary["prob"],
            abs(sol_stats_primary["resp"]),
            bins=[
                np.arange(0, 1.025, 0.025),   # x‑bins
                np.arange(0, 51, 1)           # y‑bins
            ],
            range=[[0.0, 1.0], [0.0, 50.0]],
            norm=LogNorm()
        )


        plt.colorbar(label="Count (log scale)")

        # Compute bin centers
        xcenters = 0.5 * (xedges[:-1] + xedges[1:])
        ycenters = 0.5 * (yedges[:-1] + yedges[1:])

        # Annotate each bin with count
        for i, x in enumerate(xcenters):
            for j, y in enumerate(ycenters):
                count = counts[i, j]
                if count > 0:  # avoid cluttering empty bins
                    plt.text(
                        x, y,
                        f"{int(count)}",
                        color="white",
                        ha="center",
                        va="center",
                        fontsize=5
                    )

        plt.xlabel("Probability of Obstruction")
        plt.ylabel("True Double Differenced Pseudorange Error (m)")
        plt.title("Primary Pseudorange Errors at Estimated Probability Levels")
        plt.grid()
        plt.savefig(os.path.join(out_dir, "Pseudorange Errors vs Obstruction Probability Zoomed (Heatmap).png"),dpi=300,bbox_inches="tight")
        plt.close()
except:
    traceback.print_exc()


try:
    # Probability Histogram
    plt.figure()
    plt.hist(sol_stats["prob"], bins=20, range=(0.0, 1.0))
    plt.xlabel("Probability of Obstruction")
    plt.ylabel("Count")
    plt.title("Distribution of Estimated Obstruction Probability")
    plt.grid()
    plt.savefig(os.path.join(out_dir, "Obstruction Probability Histogram.png"), dpi=300, bbox_inches="tight")
    plt.close()

    plt.show()
except:
    traceback.print_exc()


try:
    theta = np.deg2rad(sol_stats_primary["az"])
    r = 90 - sol_stats_primary["el"]

    # Classification masks
    tn = (sol_stats_primary["prob"] < 0.95) & (sol_stats_primary["resp"] < 3)
    tp = (sol_stats_primary["prob"] >= 0.95) & (sol_stats_primary["resp"] >= 3)
    fp = (sol_stats_primary["prob"] >= 0.95) & (sol_stats_primary["resp"] < 3)
    fn = (sol_stats_primary["prob"] < 0.95) & (sol_stats_primary["resp"] >= 3)

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
    plt.show()

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

    ax.legend(
        loc="upper right",
        bbox_to_anchor=(1.35, 1.1),
        title="Satellite",
        fontsize=8
    )

    plt.savefig(
        os.path.join(out_dir, "Skyplot_by_Satellite.png"),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()
    # ============================================================
    # Shared dprob histogram bounds and bins
    # ============================================================
    try:
        all_dprob = sol_stats_primary["dprob"].dropna()

        xmin = all_dprob.min()
        xmax = all_dprob.max()

        bins = np.linspace(xmin, xmax, 101)  # 100 equal-width bins
    except:
        traceback.print_exc()


    # ============================================================
    # dprob Histogram — True Negative
    # ============================================================
    try:
        data = sol_stats_primary.loc[tn, "dprob"].dropna()
        count = len(data)

        plt.figure(figsize=(10, 6))
        plt.hist(data, bins=bins, alpha=0.8)

        plt.xlim(xmin, xmax)
        plt.xlabel("Largest Obstruction Distance (dprob)")
        plt.ylabel("Count")
        plt.title(f"Largest Obstruction Distance Distribution — True Negative (n={count})")
        plt.grid()

        plt.savefig(
            os.path.join(out_dir, "dprob_histogram_true_negative.png"),
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()
    except:
        traceback.print_exc()


    # ============================================================
    # dprob Histogram — True Positive
    # ============================================================
    try:
        data = sol_stats_primary.loc[tp, "dprob"].dropna()
        count = len(data)

        plt.figure(figsize=(10, 6))
        plt.hist(data, bins=bins, alpha=0.8)

        plt.xlim(xmin, xmax)
        plt.xlabel("Largest Obstruction Distance (dprob)")
        plt.ylabel("Count")
        plt.title(f"Largest Obstruction Distance Distribution — True Positive (n={count})")
        plt.grid()

        plt.savefig(
            os.path.join(out_dir, "dprob_histogram_true_positive.png"),
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()
    except:
        traceback.print_exc()


    # ============================================================
    # dprob Histogram — False Negative
    # ============================================================
    try:
        data = sol_stats_primary.loc[fn, "dprob"].dropna()
        count = len(data)

        plt.figure(figsize=(10, 6))
        plt.hist(data, bins=bins, alpha=0.8)

        plt.xlim(xmin, xmax)
        plt.xlabel("Largest Obstruction Distance (dprob)")
        plt.ylabel("Count")
        plt.title(f"Largest Obstruction Distance Distribution — False Negative (n={count})")
        plt.grid()

        plt.savefig(
            os.path.join(out_dir, "dprob_histogram_false_negative.png"),
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()
    except:
        traceback.print_exc()


    # ============================================================
    # dprob Histogram — False Positive
    # ============================================================
    try:
        data = sol_stats_primary.loc[fp, "dprob"].dropna()
        count = len(data)

        plt.figure(figsize=(10, 6))
        plt.hist(data, bins=bins, alpha=0.8)

        plt.xlim(xmin, xmax)
        plt.xlabel("Largest Obstruction Distance (dprob)")
        plt.ylabel("Count")
        plt.title(f"Largest Obstruction Distance Distribution — False Positive (n={count})")
        plt.grid()

        plt.savefig(
            os.path.join(out_dir, "dprob_histogram_false_positive.png"),
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()
    except:
        traceback.print_exc()

    plt.figure(figsize=(10, 6))

    plt.scatter(
        sol_stats_primary.loc[tn, "dprob"],
        sol_stats_primary.loc[tn, "el"],
        s=8,
        alpha=0.6,
        label="True Negative"
    )

    plt.scatter(
        sol_stats_primary.loc[tp, "dprob"],
        sol_stats_primary.loc[tp, "el"],
        s=8,
        alpha=0.6,
        label="True Positive"
    )

    plt.scatter(
        sol_stats_primary.loc[fp, "dprob"],
        sol_stats_primary.loc[fp, "el"],
        s=8,
        alpha=0.6,
        label="False Positive"
    )

    plt.scatter(
        sol_stats_primary.loc[fn, "dprob"],
        sol_stats_primary.loc[fn, "el"],
        s=8,
        alpha=0.6,
        label="False Negative"
    )

    plt.scatter(
        sol_stats_primary.loc[tn | tp | fp | fn, "dsearchmax"],
        sol_stats_primary.loc[tn | tp | fp | fn, "el"],
        s=8,
        alpha=0.6,
        label="Max Distance Search"
    )

    plt.xlabel("dprob")
    plt.ylabel("Satellite Elevation (deg)")
    plt.title("dprob vs Satellite Elevation (Primary Frequency)")
    plt.grid()
    plt.legend()

    plt.savefig(
        os.path.join(out_dir, "dprob_vs_satellite_elevation.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

except:
    traceback.print_exc()

try:
    with open(os.path.join(out_dir, "Pseudorange Error Summary.txt"), "w") as file:
        total = sol_stats_primary.shape[0]

        file.write(f"Total Observations: {total}\n")
        file.write(f"TN: {sol_stats_primary[tn].shape[0]} : {100 * sol_stats_primary[tn].shape[0] / total:.2f}%\n")
        file.write(f"TP: {sol_stats_primary[tp].shape[0]} : {100 * sol_stats_primary[tp].shape[0] / total:.2f}%\n")
        file.write(f"FN: {sol_stats_primary[fn].shape[0]} : {100 * sol_stats_primary[fn].shape[0] / total:.2f}%\n")
        file.write(f"FP: {sol_stats_primary[fp].shape[0]} : {100 * sol_stats_primary[fp].shape[0] / total:.2f}%\n")

        # One entry per time (tow): take satellite with max residual magnitude
        sol_stats_primary_1pt = (
            sol_stats_primary
            .assign(abs_resp=lambda df: np.abs(df["resp"]))
            .sort_values("abs_resp", ascending=False)
            .drop_duplicates(subset="tow", keep="first")
            .drop(columns="abs_resp")
        )
        # Recompute masks on 1-entry-per-time dataframe
        tn = (sol_stats_primary_1pt["prob"] < 0.95) & (sol_stats_primary_1pt["resp"] < 3)
        tp = (sol_stats_primary_1pt["prob"] >= 0.95) & (sol_stats_primary_1pt["resp"] >= 3)
        fp = (sol_stats_primary_1pt["prob"] >= 0.95) & (sol_stats_primary_1pt["resp"] < 3)
        fn = (sol_stats_primary_1pt["prob"] < 0.95) & (sol_stats_primary_1pt["resp"] >= 3)

        dH = sol_stats_primary_1pt["dHeight"].astype(float)

        dH.to_csv(
            os.path.join(out_dir, "dHeight.csv"),
            index=False,
            header=["dHeight"]
        )

        dH_tn = dH[tn]
        dH_tp = dH[tp]
        dH_fn = dH[fn]
        dH_fp = dH[fp]

        
except:
    traceback.print_exc()

try:
    bins_full = np.linspace(-2000, 2000, 81)  # 50 m bins

    plt.figure(figsize=(10, 6))
    plt.hist(dH_tn, bins=bins_full, alpha=0.6, label="TN")
    plt.hist(dH_tp, bins=bins_full, alpha=0.6, label="TP")
    plt.hist(dH_fn, bins=bins_full, alpha=0.6, label="FN")
    plt.hist(dH_fp, bins=bins_full, alpha=0.6, label="FP")

    plt.xlabel("dHeight (m)")
    plt.ylabel("Count")
    plt.title("Classification Counts vs Height (-2000 to 2000 m)")
    plt.legend()
    plt.grid()

    plt.savefig(
        os.path.join(out_dir, "dHeight_Classification_Hist_Full.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
except:
    traceback.print_exc()

try:
    bins_zoom = np.linspace(-50, 50, 100)  # 5 m bins

    plt.figure(figsize=(10, 6))
    plt.hist(dH_tn, bins=bins_zoom, alpha=0.6, label="TN")
    plt.hist(dH_tp, bins=bins_zoom, alpha=0.6, label="TP")
    plt.hist(dH_fn, bins=bins_zoom, alpha=0.6, label="FN")
    plt.hist(dH_fp, bins=bins_zoom, alpha=0.6, label="FP")

    plt.xlabel("dHeight (m)")
    plt.ylabel("Count")
    plt.title("Classification Counts vs Height (-200 to 200 m)")
    plt.legend()
    plt.grid()

    plt.savefig(
        os.path.join(out_dir, "dHeight_Classification_Hist_Zoom.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
except:
    traceback.print_exc()

#Trajectory map
# Create GeoDataFrame for RTK
gdf_rtk = gpd.GeoDataFrame(
    rtk,
    geometry=gpd.points_from_xy(rtk["Lon_deg"], rtk["Lat_deg"]),
    crs="EPSG:4326"
)

gdf_truth = gpd.GeoDataFrame(
    truth_interp,
    geometry=gpd.points_from_xy(truth_interp["Lon_deg"], truth_interp["Lat_deg"]),
    crs="EPSG:4326"
)

gdf_rtk = gdf_rtk.to_crs(epsg=3857)
gdf_truth = gdf_truth.to_crs(epsg=3857)

fig, ax = plt.subplots(figsize=(10, 10))

gdf_truth.plot(
    ax=ax,
    color="green",
    linewidth=5,
    label="Truth"
)

gdf_rtk.plot(
    ax=ax,
    color="red",
    linewidth=1,
    alpha=0.8,
    label="Solution"
)

ctx.add_basemap(
    ax,
    source=ctx.providers.OpenStreetMap.Mapnik
)

ax.set_title("Solution vs Truth (Map View)")
ax.set_axis_off()
ax.legend()

plt.savefig(
    os.path.join(out_dir, "Map View - Solution vs Truth.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ============================
# CDF Metrics and Plots
# ============================
E_abs = np.abs(rtk["E_err"].astype(float).values)
N_abs = np.abs(rtk["N_err"].astype(float).values)
U_abs = np.abs(rtk["U_err"].astype(float).values)
H = np.abs(rtk["Horz_err"].astype(float).values)

print("\nSolution quality (Q) breakdown:")
q_counts = rtk["Q"].value_counts(dropna=False).sort_index()
print(q_counts)

has_solution = (
    rtk["Q"].notna() &
    (rtk["Q"] > 0) &
    np.isfinite(rtk["Lat_deg"]) &
    np.isfinite(rtk["Lon_deg"]) &
    np.isfinite(rtk["Height_m"])
)
epochs_with_solution = int(has_solution.sum())
total_epochs = int(len(rtk))
solution_rate = epochs_with_solution / total_epochs if total_epochs else np.nan

# Percentiles by Q using matched errors
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

        # --- NEW columns in the CSV ---
        "epochs_with_solution": epochs_with_solution,
        "total_epochs": total_epochs,
        "solution_rate": float(solution_rate),
    })

vals = rtk["Horz_err"].astype(float).values
vals = vals[np.isfinite(vals)]


q_rows.append({
    "Q": "All",
    "n": int(vals.size),
    "P50": float(np.nanpercentile(vals, 50)),
    "P95": float(np.nanpercentile(vals, 95)),
    "P99": float(np.nanpercentile(vals, 99)),
    "Max": float(np.nanmax(vals)),

    # --- NEW columns in the CSV ---
    "epochs_with_solution": epochs_with_solution,
    "total_epochs": total_epochs,
    "solution_rate": float(solution_rate),
})

q_table = pd.DataFrame(q_rows)
q_table.to_csv(os.path.join(out_dir, "HorzError_percentiles_by_Q.csv"), index=False)


# Build table
q_table = pd.DataFrame(q_rows)

# Nice formatting for console
q_table[["P50", "P95", "P99", "Max", "solution_rate"]] = q_table[
    ["P50", "P95", "P99", "Max", "solution_rate"]
].round(4)

print("\n==== Horz_err percentiles by Q ====")
print(q_table.to_string(index=False))

print(f"\nSolution availability: {epochs_with_solution} / {total_epochs} = {solution_rate:.3%}")



#CDF Plotting
CDF_PCT = 99

plt.figure()
for data, lab in zip([E_abs, N_abs, U_abs], ["|E|", "|N|", "|U|"]):
    xs, ys = ecdf(data)
    if xs.size:
        plt.plot(xs, ys, label=lab)
plt.xlim(p_xlim(np.r_[E_abs, N_abs, U_abs], CDF_PCT))
plt.xlabel("Absolute error (m)")
plt.ylabel("CDF")
plt.title(f"CDF Absolute ENU Component Errors (0–P{CDF_PCT})")
plt.grid(True)
plt.legend()
plt.savefig(os.path.join(out_dir, f"CDF_ENU_abs_P{CDF_PCT}.png"), dpi=300, bbox_inches="tight")
plt.close()


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
