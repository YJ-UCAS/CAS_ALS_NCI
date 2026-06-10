"""
ALS NCI estimation for a single LAS/LAZ file.

This script estimates NCI from discrete return ALS data using three methods:
LX, CC, and CLX.

Outputs are GeoTIFF files.
"""

import os
import numpy as np
import laspy
import rasterio
import warnings
from rasterio.transform import from_origin
from rasterio.crs import CRS


# =============================================================================
# Basic utilities
# =============================================================================

def safe_divide(a, b):
    """Divide arrays and return NaN where the denominator is zero."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    out = np.full_like(a, np.nan, dtype=float)
    mask = b != 0
    out[mask] = a[mask] / b[mask]
    return out

def read_las_file(las_path):
    """Read a LAS or LAZ file."""
    return laspy.read(las_path)

def normalize_height(las_data, dem_path=None, already_normalized=False):
    """
    Normalize point heights using a DEM.

    If already_normalized is True, the z values in the LAS file are assumed to
    be normalized heights.
    """
    if already_normalized:
        crs = las_data.header.parse_crs()
        return las_data, crs

    if dem_path is None:
        raise ValueError("A DEM is required unless the LAS file is already height normalized.")

    x = np.asarray(las_data.x, dtype=float)
    y = np.asarray(las_data.y, dtype=float)
    z = np.asarray(las_data.z, dtype=float)

    with rasterio.open(dem_path) as dem:
        dem_array = dem.read(1).astype(float)
        if dem.nodata is not None:
            dem_array[dem_array == dem.nodata] = np.nan

        transform = dem.transform
        crs = dem.crs

        # Convert point coordinates to DEM row and column indices
        col = ((x - transform.c) / transform.a).astype(int)
        row = ((y - transform.f) / transform.e).astype(int)

        inside = ((row >= 0) & (row < dem.height) & (col >= 0) & (col < dem.width))

        dem_values = np.full(len(z), np.nan, dtype=float)
        dem_values[inside] = dem_array[row[inside], col[inside]]

    valid = np.isfinite(dem_values)
    las_data = las_data[valid]
    las_data.z = z[valid] - dem_values[valid]

    return las_data, crs


def get_scan_angle(las_data):
    """Get scan angle in degrees from LAS 1.0 to 1.4 files."""
    dimension_names = set(las_data.point_format.dimension_names)

    if "scan_angle" in dimension_names:
        angle = np.asarray(las_data.scan_angle, dtype=float)
        if np.nanmax(np.abs(angle)) > 90:
            angle = angle * 0.006
        return angle

    if "scan_angle_rank" in dimension_names:
        return np.asarray(las_data.scan_angle_rank, dtype=float)

    return None


def filter_scan_angle(las_data, scan_angle_threshold=None):
    """Remove points with absolute scan angles larger than the threshold."""
    if scan_angle_threshold is None:
        return las_data

    angle = get_scan_angle(las_data)
    if angle is None:
        print("Warning: scan angle was not found. Scan angle filtering was skipped.")
        return las_data

    mask = np.abs(angle) <= scan_angle_threshold
    return las_data[mask]


def get_las_extent(las_data, resolution):
    """Create a raster extent from LAS coordinates and snap it to the resolution."""
    x_min = np.floor(np.nanmin(las_data.x) / resolution) * resolution
    x_max = np.ceil(np.nanmax(las_data.x) / resolution) * resolution
    y_min = np.floor(np.nanmin(las_data.y) / resolution) * resolution
    y_max = np.ceil(np.nanmax(las_data.y) / resolution) * resolution
    return [x_min, x_max, y_min, y_max]


def subset_las(las_data, mask):
    """Return a subset of LAS points."""
    return las_data[mask]


# =============================================================================
# Raster aggregation
# =============================================================================

def aggregate_points_to_grid(las_data, attribute, extent, resolution):
    """
    Aggregate LAS point attributes to a raster grid using sum.

    attribute options: count, intensity, inverse_number_of_returns
    """
    x_min, x_max, y_min, y_max = extent
    n_cols = int(np.ceil((x_max - x_min) / resolution))
    n_rows = int(np.ceil((y_max - y_min) / resolution))

    grid = np.zeros((n_rows, n_cols), dtype=float)

    if len(las_data.x) == 0:
        return grid

    if attribute == "count":
        values = np.ones(len(las_data.x), dtype=float)
    elif attribute == "intensity":
        values = np.asarray(las_data.intensity, dtype=float)
    elif attribute == "inverse_number_of_returns":
        values = 1.0 / np.asarray(las_data.number_of_returns, dtype=float)
    else:
        raise ValueError("Unsupported attribute: {}".format(attribute))

    x = np.asarray(las_data.x, dtype=float)
    y = np.asarray(las_data.y, dtype=float)

    col = np.floor((x - x_min) / resolution).astype(int)
    row_from_bottom = np.floor((y - y_min) / resolution).astype(int)
    row = n_rows - 1 - row_from_bottom

    valid = (row >= 0) & (row < n_rows) & (col >= 0) & (col < n_cols)
    np.add.at(grid, (row[valid], col[valid]), values[valid])

    return grid


def select_layer(las_data, layer, canopy_height):
    """Select all, canopy, or ground points."""
    if layer == "all":
        return las_data
    if layer == "canopy":
        return las_data[las_data.z > canopy_height]
    if layer == "ground":
        return las_data[las_data.z <= canopy_height]
    raise ValueError("Unsupported layer: {}".format(layer))


def select_return_type(las_data, return_type):
    """Select points by return type."""
    return_number = np.asarray(las_data.return_number)
    number_of_returns = np.asarray(las_data.number_of_returns)

    if return_type == "all":
        mask = np.ones(len(return_number), dtype=bool)
    elif return_type == "single":
        mask = number_of_returns == 1
    elif return_type == "first":
        mask = (number_of_returns != 1) & (return_number == 1)
    elif return_type == "last":
        mask = (number_of_returns != 1) & (return_number == number_of_returns)
    elif return_type == "intermediate":
        mask = (number_of_returns != 1) & (return_number != 1) & (return_number != number_of_returns)
    elif return_type == "pulse_first":
        mask = return_number == 1
    else:
        raise ValueError("Unsupported return type: {}".format(return_type))

    return las_data[mask]


def aggregate_lpm_component(las_data, layer, return_type, attribute, extent, resolution, canopy_height):
    """Aggregate one component needed for laser penetration index calculation."""
    data = select_layer(las_data, layer, canopy_height)
    data = select_return_type(data, return_type)
    return aggregate_points_to_grid(data, attribute, extent, resolution)


# =============================================================================
# Gap fraction from laser penetration indices
# =============================================================================

def calculate_gap_fraction(las_data, extent, resolution=2, canopy_height=2, lpi="BL", min_returns_per_cell=4):
    """
    Calculate gap fraction using a selected laser penetration index.

    Supported LPI options:
    D2, ACI, FCI, LCI, SCI, RI, FCI_RI, BL

    Note:
    This function calculates all intermediate components required by the supported
    LPIs for clarity and easy comparison among indices. If only one LPI is needed,
    users can comment out the components unrelated to that LPI to improve efficiency.
    """
    lpi = lpi.upper()

    all_all_count = aggregate_lpm_component(las_data, "all", "all", "count", extent, resolution, canopy_height)
    all_pulse_count = aggregate_lpm_component(las_data, "all", "pulse_first", "count", extent, resolution, canopy_height)

    canopy_all_count = aggregate_lpm_component(las_data, "canopy", "all", "count", extent, resolution, canopy_height)
    canopy_pulse_count = aggregate_lpm_component(las_data, "canopy", "pulse_first", "count", extent, resolution, canopy_height)
    canopy_pulse_inverse_returns = aggregate_lpm_component(las_data, "canopy", "pulse_first", "inverse_number_of_returns", extent, resolution, canopy_height)

    all_single_count = aggregate_lpm_component(las_data, "all", "single", "count", extent, resolution, canopy_height)
    all_first_count = aggregate_lpm_component(las_data, "all", "first", "count", extent, resolution, canopy_height)
    all_last_count = aggregate_lpm_component(las_data, "all", "last", "count", extent, resolution, canopy_height)

    canopy_single_count = aggregate_lpm_component(las_data, "canopy", "single", "count", extent, resolution, canopy_height)
    canopy_first_count = aggregate_lpm_component(las_data, "canopy", "first", "count", extent, resolution, canopy_height)
    canopy_last_count = aggregate_lpm_component(las_data, "canopy", "last", "count", extent, resolution, canopy_height)

    all_all_intensity = aggregate_lpm_component(las_data, "all", "all", "intensity", extent, resolution, canopy_height)
    canopy_all_intensity = aggregate_lpm_component(las_data, "canopy", "all", "intensity", extent, resolution, canopy_height)
    ground_single_intensity = aggregate_lpm_component(las_data, "ground", "single", "intensity", extent, resolution, canopy_height)
    ground_last_intensity = aggregate_lpm_component(las_data, "ground", "last", "intensity", extent, resolution, canopy_height)
    all_single_intensity = aggregate_lpm_component(las_data, "all", "single", "intensity", extent, resolution, canopy_height)
    all_first_intensity = aggregate_lpm_component(las_data, "all", "first", "intensity", extent, resolution, canopy_height)
    all_intermediate_intensity = aggregate_lpm_component(las_data, "all", "intermediate", "intensity", extent, resolution, canopy_height)
    all_last_intensity = aggregate_lpm_component(las_data, "all", "last", "intensity", extent, resolution, canopy_height)

    gf_dict = {}
    gf_dict["D2"] = 1 - safe_divide(canopy_pulse_inverse_returns, all_pulse_count)
    gf_dict["ACI"] = 1 - safe_divide(canopy_all_count, all_all_count)
    gf_dict["FCI"] = 1 - safe_divide(canopy_single_count + canopy_first_count, all_single_count + all_first_count)
    gf_dict["LCI"] = 1 - safe_divide(canopy_single_count + canopy_last_count, all_single_count + all_last_count)
    gf_dict["SCI"] = 1 - safe_divide(
        canopy_single_count + canopy_first_count + canopy_last_count,
        all_single_count + all_first_count + all_last_count,
    )
    gf_dict["RI"] = 1 - safe_divide(canopy_all_intensity, all_all_intensity)
    gf_dict["FCI_RI"] = 0.5 * (gf_dict["ACI"] + gf_dict["RI"])

    bl_num = safe_divide(ground_single_intensity, all_all_intensity) + np.sqrt(
        safe_divide(ground_last_intensity, all_all_intensity)
    )
    bl_den = safe_divide(all_single_intensity + all_first_intensity, all_all_intensity) + np.sqrt(
        safe_divide(all_intermediate_intensity + all_last_intensity, all_all_intensity)
    )
    gf_dict["BL"] = safe_divide(bl_num, bl_den)

    if lpi not in gf_dict:
        raise ValueError("Unsupported LPI: {}. Choose from {}".format(lpi, list(gf_dict.keys())))

    gf = gf_dict[lpi]
    gf[all_all_count < min_returns_per_cell] = np.nan
    gf = np.clip(gf, 0, 1)

    return gf


def calculate_large_gap_fraction(las_data, extent, resolution=30, canopy_height=2):
    """Calculate the fraction of first returns reaching below the canopy height."""
    first_returns = select_return_type(las_data, "pulse_first")
    all_first = aggregate_points_to_grid(first_returns, "count", extent, resolution)
    ground_first = aggregate_points_to_grid(first_returns[first_returns.z <= canopy_height], "count", extent, resolution)

    large_gap_fraction = safe_divide(ground_first, all_first)
    large_gap_fraction[all_first == 0] = np.nan
    return large_gap_fraction


# =============================================================================
# NCI methods
# =============================================================================

def pad_to_block_size(array, block_size):
    """Pad a raster so that rows and columns are divisible by block_size."""
    rows, cols = array.shape
    pad_rows = int(np.ceil(rows / block_size) * block_size - rows)
    pad_cols = int(np.ceil(cols / block_size) * block_size - cols)
    return np.pad(array, ((0, pad_rows), (0, pad_cols)), constant_values=np.nan)

def block_nanmean(array, block_size):
    """Calculate block mean while ignoring NaN values."""
    rows, cols = array.shape
    new_shape = (rows // block_size, block_size, cols // block_size, block_size)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(array.reshape(new_shape), axis=(1, 3))

def block_nanstd(array, block_size):
    """Calculate block standard deviation while ignoring NaN values."""
    array = pad_to_block_size(array, block_size)
    rows, cols = array.shape
    reshaped = array.reshape(rows // block_size, block_size, cols // block_size, block_size)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanstd(reshaped, axis=(1, 3))

def calculate_lx_nci(gf_segment, nci_resolution=30, segment_size=2, min_gap_fraction=0.001):
    """Calculate NCI using the LX method."""
    block_size = int(round(nci_resolution / segment_size))
    if not np.isclose(block_size * segment_size, nci_resolution):
        raise ValueError("nci_resolution must be an integer multiple of segment_size.")

    gf = np.array(gf_segment, dtype=float)
    gf[gf < min_gap_fraction] = min_gap_fraction

    numerator = np.log(block_nanmean(gf, block_size))
    denominator = block_nanmean(np.log(gf), block_size)
    nci_lx = safe_divide(numerator, denominator)

    gf_mean = block_nanmean(gf, block_size)
    gf_std = block_nanstd(gf, block_size)
    cv_gf = safe_divide(gf_std, gf_mean)

    return nci_lx, cv_gf


def calculate_cc_nci(gf_grid, large_gap_fraction, min_gap_fraction=0.001):
    """Calculate NCI using the CC method."""
    fm = np.array(gf_grid, dtype=float)
    fm[fm < min_gap_fraction] = min_gap_fraction

    fmr = safe_divide(fm - large_gap_fraction, 1 - large_gap_fraction)
    fmr[fmr < min_gap_fraction] = min_gap_fraction
    fmr[fmr > fm] = fm[fmr > fm]

    nci_cc = safe_divide(np.log(fm), np.log(fmr)) * safe_divide(1 - fmr, 1 - fm)
    return nci_cc, fmr


def calculate_clx_nci(gf_segment, cc_segment, nci_resolution=30, segment_size=2, min_gap_fraction=0.001):
    """Calculate NCI using the CLX method."""
    block_size = int(round(nci_resolution / segment_size))
    if not np.isclose(block_size * segment_size, nci_resolution):
        raise ValueError("nci_resolution must be an integer multiple of segment_size.")

    gf = np.array(gf_segment, dtype=float)
    gf[gf < min_gap_fraction] = min_gap_fraction

    cc = np.array(cc_segment, dtype=float)
    cc[~np.isfinite(cc)] = 1.0
    cc[cc <= 0] = 1.0

    numerator = np.log(block_nanmean(gf, block_size))
    denominator = block_nanmean(np.log(gf) / cc, block_size)
    nci_clx = safe_divide(numerator, denominator)

    return nci_clx


def calculate_pulse_density(las_data, extent, resolution=30):
    """Calculate pulse density from first returns."""
    first_returns = select_return_type(las_data, "pulse_first")
    pulse_count = aggregate_points_to_grid(first_returns, "count", extent, resolution)
    density = pulse_count / (resolution * resolution)
    density[density == 0] = np.nan
    return density


# =============================================================================
# Output
# =============================================================================

def save_geotiff(array, output_path, extent, resolution, crs=None, nodata=-9999):
    """Save a 2D array as a GeoTIFF."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    x_min, x_max, y_min, y_max = extent
    output = np.array(array, dtype=np.float32)
    output[~np.isfinite(output)] = nodata

    transform = from_origin(x_min, y_max, resolution, resolution)

    with rasterio.open(
        output_path,
        "w",
        driver="GTiff",
        height=output.shape[0],
        width=output.shape[1],
        count=1,
        dtype="float32",
        crs=crs,
        transform=transform,
        nodata=nodata,
        compress="lzw",
    ) as dst:
        dst.write(output, 1)


def calculate_nci_from_las(
    las_path,
    output_dir,
    dem_path=None,
    already_normalized=False,
    epsg_code=None,
    lpi="BL",
    methods=("LX", "CC", "CLX"),
    nci_resolution=30,
    segment_size=2,
    canopy_height=2,
    scan_angle_threshold=14,
    min_gap_fraction=0.001,
    min_returns_per_cell=4,
    save_intermediate=True,
):
    """Main workflow for estimating NCI from one LAS/LAZ file."""
    os.makedirs(output_dir, exist_ok=True)

    las_data = read_las_file(las_path)
    las_data, crs = normalize_height(las_data, dem_path, already_normalized)
    las_data = filter_scan_angle(las_data, scan_angle_threshold)

    if crs is None and epsg_code is not None:
        crs = CRS.from_epsg(epsg_code)

    base_name = os.path.splitext(os.path.basename(las_path))[0]
    extent = get_las_extent(las_data, nci_resolution)

    gf_segment = calculate_gap_fraction(
        las_data,
        extent,
        resolution=segment_size,
        canopy_height=canopy_height,
        lpi=lpi,
        min_returns_per_cell=min_returns_per_cell,
    )
    gf_nci = calculate_gap_fraction(
        las_data,
        extent,
        resolution=nci_resolution,
        canopy_height=canopy_height,
        lpi=lpi,
        min_returns_per_cell=min_returns_per_cell,
    )

    if save_intermediate:
        save_geotiff(gf_nci, os.path.join(output_dir, base_name + "_GF.tif"), extent, nci_resolution, crs)
        pulse_density = calculate_pulse_density(las_data, extent, resolution=nci_resolution)
        save_geotiff(pulse_density, os.path.join(output_dir, base_name + "_PulseDensity.tif"), extent, nci_resolution, crs)

    if "LX" in methods:
        nci_lx, cv_gf = calculate_lx_nci(gf_segment, nci_resolution, segment_size, min_gap_fraction)
        save_geotiff(nci_lx, os.path.join(output_dir, base_name + "_NCI_LX.tif"), extent, nci_resolution, crs)
        if save_intermediate:
            save_geotiff(cv_gf, os.path.join(output_dir, base_name + "_CV.tif"), extent, nci_resolution, crs)

    if "CC" in methods:
        large_gap_nci = calculate_large_gap_fraction(las_data, extent, resolution=nci_resolution, canopy_height=canopy_height)
        nci_cc, fmr = calculate_cc_nci(gf_nci, large_gap_nci, min_gap_fraction)
        save_geotiff(nci_cc, os.path.join(output_dir, base_name + "_NCI_CC.tif"), extent, nci_resolution, crs)

    if "CLX" in methods:
        large_gap_segment = calculate_large_gap_fraction(las_data, extent, resolution=segment_size, canopy_height=canopy_height)
        cc_segment, _ = calculate_cc_nci(gf_segment, large_gap_segment, min_gap_fraction)
        nci_clx = calculate_clx_nci(gf_segment, cc_segment, nci_resolution, segment_size, min_gap_fraction)
        save_geotiff(nci_clx, os.path.join(output_dir, base_name + "_NCI_CLX.tif"), extent, nci_resolution, crs)

    print("Finished:", las_path)
    print("Outputs saved to:", output_dir)


# =============================================================================
# User parameters
# =============================================================================

if __name__ == "__main__":
    # Input files
    LAS_FILE = r"sample_PCD.laz"
    DEM_FILE = r"H:sample_DEM.tif"        # Set to None if LAS is already height normalized
    OUTPUT_DIR = r"outputs"

    # Processing parameters
    ALREADY_NORMALIZED = False           # True if LAS z is already normalized height
    EPSG_CODE = None                     # Example: 32618, only needed if CRS is missing
    LPI = "BL"                           # D2, ACI, FCI, LCI, SCI, RI, FCI_RI, BL
    METHODS = ("LX", "CC", "CLX")

    NCI_RESOLUTION = 30                  # Output NCI resolution, in meters
    SEGMENT_SIZE = 2                     # Segment size for LX and CLX, in meters
    CANOPY_HEIGHT = 2                    # Height threshold separating canopy and ground returns
    SCAN_ANGLE_THRESHOLD = 14            # Set to None if no scan angle filtering is needed
    MIN_GAP_FRACTION = 0.001
    MIN_RETURNS_PER_CELL = 4
    SAVE_INTERMEDIATE = True

    calculate_nci_from_las(
        las_path=LAS_FILE,
        output_dir=OUTPUT_DIR,
        dem_path=DEM_FILE,
        already_normalized=ALREADY_NORMALIZED,
        epsg_code=EPSG_CODE,
        lpi=LPI,
        methods=METHODS,
        nci_resolution=NCI_RESOLUTION,
        segment_size=SEGMENT_SIZE,
        canopy_height=CANOPY_HEIGHT,
        scan_angle_threshold=SCAN_ANGLE_THRESHOLD,
        min_gap_fraction=MIN_GAP_FRACTION,
        min_returns_per_cell=MIN_RETURNS_PER_CELL,
        save_intermediate=SAVE_INTERMEDIATE,
    )
