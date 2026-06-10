# ALS NCI Estimation

This repository provides a simple Python script for estimating nadir clumping index (NCI) from a single discrete return airborne laser scanning (ALS) LAS/LAZ file.

The script implements three NCI estimation methods:

* LX: logarithmic gap fraction averaging method
* CC: gap size distribution method
* CLX: combined gap size and logarithmic averaging method

The output products are saved as GeoTIFF files.

## 1. Overview

The workflow includes the following steps:

1. Read a LAS/LAZ point cloud.
2. Normalize point heights using a DEM, or use an already height normalized point cloud.
3. Optionally filter points by scan angle threshold.
4. Calculate gap fraction using a selected laser penetration index.
5. Estimate NCI using LX, CC, and/or CLX methods.
6. Save NCI maps and optional intermediate products as GeoTIFF files.

This script is designed for processing a single LAS/LAZ file. Batch processing can be implemented by calling the main function for multiple files.

## 2. Requirements

The code requires Python 3 and the following packages:

```bash
os
numpy
laspy
rasterio
```

## 3. Input data

The script requires:

* A LAS or LAZ point cloud file
* A DEM raster file if the point cloud is not height normalized

If the LAS/LAZ file already stores normalized heights, set:

```python
ALREADY_NORMALIZED = True
DEM_FILE = None
```

If the LAS/LAZ file stores original elevations, provide a DEM:

```python
ALREADY_NORMALIZED = False
DEM_FILE = r"path_to_dem.tif"
```

The LAS/LAZ file and DEM should use the same projected coordinate system.

## 4. Main parameters

The user parameters are set at the bottom of the script.

```python
LAS_FILE = r"sample_PCD.laz"
DEM_FILE = r"sample_DEM.tif"
OUTPUT_DIR = r"outputs"

ALREADY_NORMALIZED = False
EPSG_CODE = None
LPI = "BL"
METHODS = ("LX", "CC", "CLX")

NCI_RESOLUTION = 30
SEGMENT_SIZE = 2
CANOPY_HEIGHT = 2
SCAN_ANGLE_THRESHOLD = 14
MIN_GAP_FRACTION = 0.001
MIN_RETURNS_PER_CELL = 4
SAVE_INTERMEDIATE = True
```

### Parameter descriptions

| Parameter              | Description                                                     |
| ---------------------- | --------------------------------------------------------------- |
| `LAS_FILE`             | Path to the input LAS/LAZ file                                  |
| `DEM_FILE`             | Path to the DEM used for height normalization                   |
| `OUTPUT_DIR`           | Folder for saving output GeoTIFF files                          |
| `ALREADY_NORMALIZED`   | Set to `True` if LAS z values are already normalized heights    |
| `EPSG_CODE`            | EPSG code used when CRS information is missing                  |
| `LPI`                  | Laser penetration index used to calculate gap fraction          |
| `METHODS`              | NCI methods to run: `LX`, `CC`, and/or `CLX`                    |
| `NCI_RESOLUTION`       | Output NCI resolution, in meters                                |
| `SEGMENT_SIZE`         | Segment size for LX and CLX, in meters                          |
| `CANOPY_HEIGHT`        | Height threshold separating canopy and ground returns           |
| `SCAN_ANGLE_THRESHOLD` | Maximum absolute scan angle used for filtering                  |
| `MIN_GAP_FRACTION`     | Minimum gap fraction used to avoid log zero                     |
| `MIN_RETURNS_PER_CELL` | Minimum number of returns required for gap fraction calculation |
| `SAVE_INTERMEDIATE`    | Whether to save intermediate products                           |

## 5. Supported laser penetration indices

The script supports the following laser penetration indices for gap fraction estimation:

| LPI      |
| -------- |
| `D2`     |
| `ACI`    |
| `FCI`    |
| `LCI`    |
| `SCI`    |
| `RI`     |
| `FCI_RI` |
| `BL`     |

The default option is:

```python
LPI = "BL"
```

For clarity and comparison among indices, the function `calculate_gap_fraction()` calculates intermediate components for all supported LPIs. If only one LPI is needed, users can comment out unrelated component calculations to improve efficiency.

## 6. Usage

Edit the user parameters at the bottom of `als_nci.py`, then run:

```bash
python als_nci.py
```

Example:

```python
LAS_FILE = r"data/sample_PCD.laz"
DEM_FILE = r"data/sample_DEM.tif"
OUTPUT_DIR = r"outputs"

ALREADY_NORMALIZED = False
LPI = "BL"
METHODS = ("LX", "CC", "CLX")

NCI_RESOLUTION = 30
SEGMENT_SIZE = 2
CANOPY_HEIGHT = 2
SCAN_ANGLE_THRESHOLD = 14
SAVE_INTERMEDIATE = True
```

## 7. Outputs

The script saves GeoTIFF files to the output folder.

Main outputs:

| Output          | Description                        |
| --------------- | ---------------------------------- |
| `*_NCI_LX.tif`  | NCI estimated using the LX method  |
| `*_NCI_CC.tif`  | NCI estimated using the CC method  |
| `*_NCI_CLX.tif` | NCI estimated using the CLX method |

Optional intermediate outputs:

| Output               | Description                                            |
| -------------------- | ------------------------------------------------------ |
| `*_GF.tif`           | Gap fraction at NCI resolution                         |
| `*_PulseDensity.tif` | Pulse density at NCI resolution                        |
| `*_CV.tif`           | Coefficient of variation of segment level gap fraction |

Intermediate outputs are saved only when:

```python
SAVE_INTERMEDIATE = True
```

## 8. Important notes

* The LAS/LAZ point cloud and DEM should be in the same projected coordinate system.
* The NCI resolution should be an integer multiple of the segment size.
* A DEM is required unless the point cloud is already height normalized.
* The scan angle threshold can be set to `None` if scan angle filtering is not needed.
* Cells with fewer returns than `MIN_RETURNS_PER_CELL` are set to `NaN` during gap fraction calculation.
* The output GeoTIFF files use `-9999` as the NoData value.
