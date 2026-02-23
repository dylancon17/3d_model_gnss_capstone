import numpy as np
import trimesh
import shapely.geometry as geom
import os
import glob

print("Using trimesh version:", trimesh.__version__)

# -----------------------------
# CONFIGURATION
# -----------------------------
directory = r"C:\capstone\dsm_tiles\DSM_CGY_5x5km_Res1m_SnappedProperly"

# Original DSM size
WIDTH_ORIG  = 5000
HEIGHT_ORIG = 5000

# Downsampling
DOWNSAMPLE = 5

# Downsampled grid size
WIDTH  = WIDTH_ORIG  // DOWNSAMPLE
HEIGHT = HEIGHT_ORIG // DOWNSAMPLE

# Target print footprint (mm)
FOOTPRINT_MM = 120.0   # 12 cm × 12 cm

# Derived spacing
HORIZONTAL_SPACING_MM = FOOTPRINT_MM / (WIDTH - 1)

# Vertical scaling
BASE_THICKNESS_MM   = 10.0
MAX_PRINT_HEIGHT_MM = 150.0
Z_INPUT_MAX_CM = 40000.0
Z_SCALE = MAX_PRINT_HEIGHT_MM / Z_INPUT_MAX_CM

# Wedge geometry
WEDGE_WIDTH_MM  = 12.0
WEDGE_HEIGHT_MM = 5.0
WEDGE_DEPTH_MM  = 20.0
CLEARANCE_MM    = 0.25

# -----------------------------
# TRIANGULAR PRISM
# -----------------------------
def triangular_prism(width, height, depth):
    verts = np.array([
        [0, 0, 0],
        [width, 0, 0],
        [0, height, 0],
        [0, 0, depth],
        [width, 0, depth],
        [0, height, depth],
    ])

    faces = np.array([
        [0, 1, 2],
        [3, 5, 4],
        [0, 3, 1],
        [1, 3, 4],
        [1, 4, 2],
        [2, 4, 5],
        [2, 5, 0],
        [0, 5, 3],
    ])

    return trimesh.Trimesh(vertices=verts, faces=faces, process=True)

# -----------------------------
# PROCESS EACH TILE
# -----------------------------
bin_files = glob.glob(os.path.join(directory, "*.bin"))

for file in bin_files:
    print(f"Processing {file}...")

    # -----------------------------
    # LOAD + DOWNSAMPLE HEIGHTFIELD
    # -----------------------------
    data = np.fromfile(file, dtype=np.uint16)
    Z_cm = data.reshape((HEIGHT_ORIG, WIDTH_ORIG)).astype(np.float32)

    Z_cm = Z_cm[::DOWNSAMPLE, ::DOWNSAMPLE]
    Z_mm = Z_cm * Z_SCALE

    # -----------------------------
    # BUILD TERRAIN MESH
    # -----------------------------
    x = np.arange(WIDTH)  * HORIZONTAL_SPACING_MM
    y = np.arange(HEIGHT) * HORIZONTAL_SPACING_MM
    xx, yy = np.meshgrid(x, y)

    vertices = np.column_stack((
        xx.ravel(),
        yy.ravel(),
        Z_mm.ravel()
    ))

    faces = []
    for j in range(HEIGHT - 1):
        for i in range(WIDTH - 1):
            idx = j * WIDTH + i
            faces.append([idx, idx + 1, idx + WIDTH])
            faces.append([idx + 1, idx + WIDTH + 1, idx + WIDTH])

    terrain = trimesh.Trimesh(
        vertices=vertices,
        faces=np.asarray(faces),
        process=False
    )

    print("terrain size (mm):", FOOTPRINT_MM, "×", FOOTPRINT_MM)

    # -----------------------------
    # BASE WITH FEMALE SLOTS
    # -----------------------------
    outer = geom.Polygon([
        (0, 0),
        (FOOTPRINT_MM, 0),
        (FOOTPRINT_MM, FOOTPRINT_MM),
        (0, FOOTPRINT_MM)
    ])

    slot_w = WEDGE_WIDTH_MM + 2 * CLEARANCE_MM
    slot_d = WEDGE_DEPTH_MM + CLEARANCE_MM

    south_slot = geom.Polygon([
        (FOOTPRINT_MM/2 - slot_w/2, 0),
        (FOOTPRINT_MM/2 + slot_w/2, 0),
        (FOOTPRINT_MM/2 + slot_w/2, slot_d),
        (FOOTPRINT_MM/2 - slot_w/2, slot_d),
    ])

    west_slot = geom.Polygon([
        (0, FOOTPRINT_MM/2 - slot_w/2),
        (slot_d, FOOTPRINT_MM/2 - slot_w/2),
        (slot_d, FOOTPRINT_MM/2 + slot_w/2),
        (0, FOOTPRINT_MM/2 + slot_w/2),
    ])

    footprint = outer.difference(south_slot.union(west_slot))

    base = trimesh.creation.extrude_polygon(
        footprint,
        height=BASE_THICKNESS_MM
    )
    base.apply_translation([0, 0, -BASE_THICKNESS_MM])

    # -----------------------------
    # MALE WEDGES
    # -----------------------------
    wedge = triangular_prism(
        WEDGE_WIDTH_MM,
        WEDGE_HEIGHT_MM,
        WEDGE_DEPTH_MM
    )

    north = wedge.copy()
    north.apply_transform(
        trimesh.transformations.rotation_matrix(np.pi / 2, [1, 0, 0])
    )
    north.apply_translation([
        FOOTPRINT_MM/2 - WEDGE_WIDTH_MM/2,
        FOOTPRINT_MM,
        -BASE_THICKNESS_MM
    ])

    east = wedge.copy()
    east.apply_transform(
        trimesh.transformations.rotation_matrix(-np.pi / 2, [0, 1, 0])
    )
    east.apply_translation([
        FOOTPRINT_MM,
        FOOTPRINT_MM/2 - WEDGE_WIDTH_MM/2,
        -BASE_THICKNESS_MM
    ])

    # -----------------------------
    # FINAL ASSEMBLY
    # -----------------------------
    mesh = trimesh.util.concatenate([
        terrain,
        base,
        north,
        east
    ])

    out_file = file.replace(".bin", "_120mm_ds5.stl")
    mesh.export(out_file)
    print(f"Saved {out_file}")

print("All tiles processed.")