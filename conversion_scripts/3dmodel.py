import numpy as np
import trimesh
import os
import re

# ============================================================
# CONFIGURATION
# ============================================================

BIN_FILE = r"C:\capstone\dsm_tiles\DSM_CGY_5x5km_Res1m_SnappedProperly\DSM_CGY_5x5km_1mResolution_SnappedProperly_-5989.47E_5657949.28N.bin"

DSM_SIZE = 5000          # 5 km × 5 km DSM (1 m resolution)
TILE_SIZE = 1000         # 1 km × 1 km tiles

FOOTPRINT_MM = 181.8      # Printed footprint (X/Y)
MIN_BASE_MM = 5.0       # REQUIRED minimum base thickness
MAX_PRINT_HEIGHT_MM = 49.0

# DSM vertical limits (UNCHANGED)
MIN_Z = 5675
MAX_Z = 32608
Z_INPUT_MAX_CM = MAX_Z - MIN_Z
Z_SCALE = MAX_PRINT_HEIGHT_MM / Z_INPUT_MAX_CM

# FDM parameters
LAYER_HEIGHT_MM = 0.2

# XY scale
CELL_MM = FOOTPRINT_MM / (TILE_SIZE - 1)

# Decimation (XY only, Z will be re-snapped)
DECIMATE_RATIO = 0.4   # 40% of faces retained

# ============================================================
# PARSE DSM ORIGIN
# ============================================================

fname = os.path.basename(BIN_FILE)
match = re.search(r'_([-0-9.]+)E_([-0-9.]+)N\.bin$', fname)
if not match:
    raise ValueError("Could not parse coordinates from DSM filename")

TOPLEFT_E = float(match.group(1))
TOPLEFT_N = float(match.group(2))

print(f"DSM origin (top-left): E={TOPLEFT_E}, N={TOPLEFT_N}")

# ============================================================
# HEIGHTFIELD GENERATION
# ============================================================

def build_heightfield(Z_mm):
    rows, cols = Z_mm.shape
    vertices = np.zeros((rows * cols, 3), dtype=np.float32)
    faces = []

    for y in range(rows):
        for x in range(cols):
            idx = y * cols + x
            vertices[idx] = [
                x * CELL_MM,
                y * CELL_MM,
                Z_mm[y, x]
            ]

    def vid(x, y):
        return y * cols + x

    for y in range(rows - 1):
        for x in range(cols - 1):
            faces.append([vid(x, y), vid(x + 1, y), vid(x + 1, y + 1)])
            faces.append([vid(x, y), vid(x + 1, y + 1), vid(x, y + 1)])

    return trimesh.Trimesh(
        vertices=vertices,
        faces=np.array(faces),
        process=False
    )

# ============================================================
# FDM-SAFE SOLID EXTRUSION
# ============================================================

def make_fdm_solid(surface):
    """
    Convert a terrain surface mesh into a watertight FDM-safe solid
    with a guaranteed minimum base thickness.
    """

    vertices = surface.vertices.copy()
    faces = surface.faces.copy()

    # --- Determine base thickness ---
    min_z = vertices[:, 2].min()
    base_thickness = max(MIN_BASE_MM, min_z)

    # Shift terrain upward
    vertices[:, 2] -= min_z
    vertices[:, 2] += base_thickness

    n_verts = len(vertices)

    # --- Create bottom vertices (flat at Z=0) ---
    bottom_vertices = vertices.copy()
    bottom_vertices[:, 2] = 0.0

    # --- Combine vertices ---
    all_vertices = np.vstack([vertices, bottom_vertices])

    # --- Bottom faces (reverse winding!) ---
    bottom_faces = faces[:, ::-1] + n_verts

    # --- Build side walls ---
    edge_map = {}
    for tri in faces:
        for i in range(3):
            a = tri[i]
            b = tri[(i + 1) % 3]
            edge = tuple(sorted((a, b)))
            edge_map[edge] = edge_map.get(edge, 0) + 1

    boundary_edges = [e for e, c in edge_map.items() if c == 1]

    side_faces = []
    for a, b in boundary_edges:
        a_top = a
        b_top = b
        a_bot = a + n_verts
        b_bot = b + n_verts

        side_faces.append([a_top, b_top, b_bot])
        side_faces.append([a_top, b_bot, a_bot])

    # --- Combine all faces ---
    all_faces = np.vstack([
        faces,
        bottom_faces,
        np.array(side_faces)
    ])

    solid = trimesh.Trimesh(
        vertices=all_vertices,
        faces=all_faces,
        process=False
    )

    return solid
# ============================================================
# LOAD DSM
# ============================================================

print("Loading DSM...")
data = np.fromfile(BIN_FILE, dtype=np.uint16)
print(f"DSM min: {data.min()}, max: {data.max()}")

Z_cm_full = data.reshape((DSM_SIZE, DSM_SIZE)).astype(np.float32)

# ============================================================
# TILE PROCESSING
# ============================================================

tile_count = 0

for ty in range(0, DSM_SIZE, TILE_SIZE):
    for tx in range(0, DSM_SIZE, TILE_SIZE):

        tile_count += 1
        print(f"Generating tile {tile_count} at pixel ({tx}, {ty})")

        tile_E = TOPLEFT_E + tx
        tile_N = TOPLEFT_N - ty

        # Extract DSM tile
        Z_cm = Z_cm_full[ty:ty + TILE_SIZE, tx:tx + TILE_SIZE] - MIN_Z
        Z_mm = Z_cm * Z_SCALE

        # Snap Z to FDM layers
        Z_mm = np.round(Z_mm / LAYER_HEIGHT_MM) * LAYER_HEIGHT_MM
        Z_mm[Z_mm < 0] = 0

        # Build heightfield
        terrain = build_heightfield(Z_mm)

        # XY decimation (Z error fixed afterward)
        target_faces = int(len(terrain.faces) * DECIMATE_RATIO)
        # Decimate in XY (Z will be re-snapped afterward)
        terrain = terrain.simplify_quadric_decimation(
            percent=DECIMATE_RATIO,
            aggression=5
        )
        # Re-snap Z after decimation (CRITICAL)
        terrain.vertices[:, 2] = (
            np.round(terrain.vertices[:, 2] / LAYER_HEIGHT_MM)
            * LAYER_HEIGHT_MM
        )

        # Build solid with guaranteed base
        mesh = make_fdm_solid(terrain)

        # Mesh cleanup (slicer safety)
        # mesh.remove_duplicate_faces()
        # mesh.remove_degenerate_faces()
        # mesh.merge_vertices()
        # mesh.remove_unreferenced_vertices()
        mesh.process(validate=True)

        out_file = os.path.join(
            os.path.dirname(BIN_FILE),
            f"tile_{tile_E:.2f}E_{tile_N:.2f}N_FDM_100mm.stl"
        )

        mesh.export(out_file)
        print(f"Saved {out_file}")

print("All tiles generated successfully.")