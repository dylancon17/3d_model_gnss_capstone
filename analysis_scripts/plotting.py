import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Transformer
import sys
import os
import geopandas as gpd
import contextily as ctx
from shapely.geometry import Point
from matplotlib.colors import LogNorm

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

# Relative time (seconds from start)
rtk["t"] = rtk["GPSTime"] - rtk["GPSTime"].iloc[0]


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
    "prob"    # obstruction probability
]

sol_stats = pd.DataFrame(sat_rows, columns=sat_cols)

# Columns that should be numeric
numeric_cols = [
    "week", "tow", "frq",
    "az", "el",
    "resp", "resc",
    "vsat", "snr", "fix", "slip",
    "lock", "outc", "slipc", "rejc",
    "scal", "prob"
]

sol_stats[numeric_cols] = sol_stats[numeric_cols].apply(
    pd.to_numeric, errors="coerce"
)

# Relative time in seconds from first epoch
sol_stats["t"] = sol_stats["tow"] - sol_stats["tow"].iloc[0]



# ============================
# LOAD RTKLIB $SAT STATS FILE FOR TRUTH
# ============================
truth_stat_rows = []

with open(truth_stats, "r") as f:
    for line in f:
        if not line.startswith("$SAT"):
            continue

        parts = line.strip().split(",")

        truth_stat_rows.append(parts[1:])  # drop "$SAT"

truth_stats = pd.DataFrame(truth_stat_rows, columns=sat_cols)

truth_stats[numeric_cols] = truth_stats[numeric_cols].apply(
    pd.to_numeric, errors="coerce"
)

# Relative time in seconds from first epoch
truth_stats["t"] = truth_stats["tow"] - truth_stats["tow"].iloc[0]



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

# Drop rows where GPSTime is NaN (e.g., the units row)
truth = truth.dropna(subset=["GPSTime", "Latitude", "Longitude", "H-Ell"])

# Truth relative time
truth["t"] = truth["GPSTime"] - truth["GPSTime"].iloc[0]


# ============================
# INTERPOLATE TRUTH TO RTK TIME
# ============================
truth_interp = pd.DataFrame()
truth_interp["t"] = rtk["t"]
truth_interp["Lat_deg"]  = np.interp(rtk["t"], truth["t"], truth["Latitude"])
truth_interp["Lon_deg"]  = np.interp(rtk["t"], truth["t"], truth["Longitude"])
truth_interp["Height_m"] = np.interp(rtk["t"], truth["t"], truth["H-Ell"])


# ============================
# ALIGN SOLUTION STATS with TRUTH STATS
# ============================
sat_error_prob = sol_stats.merge(
    truth_stats,
    on=["tow", "sat"],
    how="inner",
    suffixes=("_sol", "_truth")
)

# Remove zero-probability entries (temporary until DEM of full city)
sat_error_prob = sat_error_prob[sat_error_prob["prob_sol"] > 0.0]


# ============================
# COORDINATE TRANSFORM SETUP
# ============================
llt_to_ecef = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)

# Convert truth LLH → ECEF
truth_ecef = np.array(llt_to_ecef.transform(
    truth_interp["Lon_deg"].values,
    truth_interp["Lat_deg"].values,
    truth_interp["Height_m"].values
)).T

# Convert RTK LLH → ECEF
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
# Reference LLH in radians
lat0 = np.deg2rad(ref_lat)
lon0 = np.deg2rad(ref_lon)

# Rotation matrix ECEF → ENU
R = np.array([
    [-np.sin(lon0),              np.cos(lon0),              0],
    [-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0), np.cos(lat0)],
    [ np.cos(lat0)*np.cos(lon0),  np.cos(lat0)*np.sin(lon0), np.sin(lat0)]
])

# Compute ENU
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

# Print to console
print("\n==== RMSE Positioning Errors ====")
print(f"Horizontal RMSE : {rmse_h:.4f} m")
print(f"Vertical RMSE   : {rmse_v:.4f} m")
print(f"3D RMSE         : {rmse_3d:.4f} m")

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

# ENU Errors
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

# Horizontal error
plt.figure()
plt.plot(rtk["t"], rtk["Horz_err"])
plt.xlabel("Time (s)")
plt.ylabel("Horizontal Error (m)")
plt.title("Horizontal Error Over Time")
plt.grid()
plt.savefig(os.path.join(out_dir, "Horizontal Error Over Time.png"), dpi=300, bbox_inches="tight")
plt.close()

# Vertical error
plt.figure()
plt.plot(rtk["t"], rtk["U_err"])
plt.xlabel("Time (s)")
plt.ylabel("Vertical Error (m)")
plt.title("Vertical Error Over Time")
plt.grid()
plt.savefig(os.path.join(out_dir, "Vertical Error Over Time.png"), dpi=300, bbox_inches="tight")
plt.close()

# Histograms
xmax = np.percentile(rtk["Horz_err"], 50)
plt.figure()
plt.hist(rtk["Horz_err"], bins=10000)
plt.xlim(0, xmax)
plt.xlabel("Horizontal Error (m)")
plt.ylabel("Count")
plt.title("Horizontal Error Histogram")
plt.grid()
plt.savefig(os.path.join(out_dir, "Horizontal Error Histogram.png"), dpi=300, bbox_inches="tight")
plt.close()

xmax = np.percentile(rtk["Horz_err"], 50)
plt.figure()
plt.hist(rtk["U_err"], bins=10000)
plt.xlim(0, xmax)
plt.xlabel("Vertical Error (m)")
plt.ylabel("Count")
plt.title("Vertical Error Histogram")
plt.grid()
plt.savefig(os.path.join(out_dir, "Vertical Error Histogram.png"), dpi=300, bbox_inches="tight")
plt.close()

# Satellite count
plt.figure()
plt.plot(rtk["t"], rtk["ns"])
plt.xlabel("Time (s)")
plt.ylabel("# of Satellites")
plt.title("Satellite Count vs Time")
plt.grid()
plt.savefig(os.path.join(out_dir, "Satellite Count vs Time.png"), dpi=300, bbox_inches="tight")
plt.close()

# Sat prob vs pseudorange error
plt.figure()
plt.scatter(sat_error_prob["prob_sol"], abs(sat_error_prob["resp_truth"]), s=8, alpha=0.6)
plt.xlabel("Probability of Obstruction")
plt.ylabel("True Double Differenced Pseudorange Error (m)")
plt.title("Pseudorange Errors at Estimated Probability Levels")
plt.grid()
plt.savefig(os.path.join(out_dir, "Pseudorange Errors vs Obstruction Probability.png"),dpi=300,bbox_inches="tight")
plt.close()

# Sat prob vs carrier phase error heatmap
plt.figure()
plt.scatter(sat_error_prob["prob_sol"], abs(sat_error_prob["resc_truth"]), s=8, alpha=0.6)
plt.xlabel("Probability of Obstruction")
plt.ylabel("True Double Differenced Carrier Phase Error (m)")
plt.title("Carrier Phase Errors at Estimated Probability Levels")
plt.grid()
plt.savefig(os.path.join(out_dir, "Carrier Phase Errors vs Obstruction Probability Heatmap.png"),dpi=300,bbox_inches="tight")
plt.close()
plt.figure()

plt.hist2d(
    sat_error_prob["prob_sol"],
    abs(sat_error_prob["resp_truth"]),
    bins=[
        20,                    # probability bins
        int(300 / 10)           # 10 m bins up to 300 m
    ],
    range=[
        [0.0, 1.0],
        [0.0, 300.0]
    ],
    norm=LogNorm()              # <-- log-scaled color density
)

plt.colorbar(label="Count (log scale)")

plt.xlabel("Probability of Obstruction")
plt.ylabel("True Double Differenced Pseudorange Error (m)")
plt.title("Pseudorange Errors at Estimated Probability Levels")
plt.grid()

plt.savefig(
    os.path.join(out_dir, "Pseudorange Errors vs Obstruction Probability Heatmap (Log).png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Probability Histogram
plt.figure()
plt.hist(sat_error_prob["prob_sol"],bins=20,range=(0.0, 1.0))
plt.xlabel("Probability of Obstruction")
plt.ylabel("Count")
plt.title("Distribution of Estimated Obstruction Probability")
plt.grid()
plt.savefig(os.path.join(out_dir, "Obstruction Probability Histogram.png"),dpi=300,bbox_inches="tight")
plt.close()


plt.show()

#Trajectory map
# Create GeoDataFrame for RTK
gdf_rtk = gpd.GeoDataFrame(
    rtk,
    geometry=gpd.points_from_xy(rtk["Lon_deg"], rtk["Lat_deg"]),
    crs="EPSG:4326"  # WGS84 lat/lon
)

# Create GeoDataFrame for Truth
gdf_truth = gpd.GeoDataFrame(
    truth_interp,
    geometry=gpd.points_from_xy(truth_interp["Lon_deg"], truth_interp["Lat_deg"]),
    crs="EPSG:4326"
)

gdf_rtk = gdf_rtk.to_crs(epsg=3857)
gdf_truth = gdf_truth.to_crs(epsg=3857)

fig, ax = plt.subplots(figsize=(10, 10))

# Plot truth and RTK
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

# Add basemap
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


