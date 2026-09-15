"""Final v2.6.8 defaults; runtime paths are set explicitly by configuration."""

import re
from pathlib import Path

import numpy as np

OCSSWROOT = Path("~/ocssw").expanduser().resolve()

OCSSW_INSTALLER_DIR = Path("~/ocssw-installer").expanduser().resolve()

RUN_OCSSW_INSTALL = True

PREFERRED_OCSSW_TAG = "V2026.3"

FALLBACK_TAG_COUNT = 3

PROJECT_ROOT = Path("data").resolve()

L2_MODE = "full_swath_rhos"

NDCI_SOURCE = "rhos"

PROCESS_FULL_SCENE = False

GET_ANCILLARY = True

ALLOW_CLIMATOLOGY_FALLBACK = False

SAVE_FULL_REFLECTANCE_SPECTRA = False

COMPUTE_SHA256 = False

CYAN_CI_DETECTION_LIMIT = 0.0001

CI_DETECTION_LIMIT = 0.0001

MASK_FLAG_CHOICES = (
    "ATMFAIL",
    "NAVFAIL",
    "NAVWARN",
    "HILT",
    "CLDICE",
    "STRAYLIGHT",
    "HIGLINT",
    "HISATZEN",
    "HISOLZEN",
)

MASK_PRESETS = {
    "bloom_retaining": {
        "label": "Bloom-retaining inland (recommended)",
        "qa_flags": ("NAVFAIL",),
        "land_adjacency_pixels": 0,
        "note": "Excludes only navigation failure. HILT and CLDICE are retained "
        "because bright blooms can trigger them; inspect flags and reason "
        "maps before publication.",
    },
    "balanced_inland": {
        "label": "Balanced inland (recommended)",
        "qa_flags": ("ATMFAIL", "NAVFAIL", "HILT", "CLDICE", "STRAYLIGHT"),
        "land_adjacency_pixels": 0,
        "note": "Excludes hard failures, saturation, cloud/ice, and straylight; no "
        "automatic shoreline dilation.",
    },
    "relaxed_inland": {
        "label": "Relaxed inland (exploratory)",
        "qa_flags": ("NAVFAIL", "HILT", "CLDICE"),
        "land_adjacency_pixels": 0,
        "note": "Retains ATMFAIL, warning, glint, geometry, and straylight pixels for "
        "sensitivity analysis.",
    },
    "diagnostic": {
        "label": "Diagnostic finite-water only",
        "qa_flags": (),
        "land_adjacency_pixels": 0,
        "note": "Troubleshooting only: applies LAND/geolocation and finite-input tests but "
        "no additional QA exclusions.",
    },
    "conservative": {
        "label": "Conservative",
        "qa_flags": (
            "ATMFAIL",
            "NAVFAIL",
            "NAVWARN",
            "HILT",
            "CLDICE",
            "STRAYLIGHT",
            "HIGLINT",
            "HISATZEN",
            "HISOLZEN",
        ),
        "land_adjacency_pixels": 1,
        "note": "Original strict mask: all listed QA flags plus a one-pixel shoreline buffer.",
    },
    "cyan_comparable": {
        "label": "CyAN-comparable quantitative (LAND + CLDICE + HISATZEN)",
        "qa_flags": ("CLDICE", "HISATZEN"),
        "land_adjacency_pixels": 0,
        "note": "Closest documented post-L2 flag configuration to NASA/OBPG CyAN "
        "Level-3 binning. LAND is enforced separately through "
        "REQUIRE_L2_NONLAND=True.",
    },
    "bloom_aware_cyan": {
        "label": "Bloom-aware quantitative (cloud masked + CLDICE bloom recovery)",
        "qa_flags": ("HISATZEN", "NAVFAIL"),
        "land_adjacency_pixels": 0,
        "note": "Primary quantitative mask. LAND, HISATZEN and NAVFAIL are excluded. "
        "CLDICE is excluded unless the pixel independently passes the "
        "standard CIcyano spectral gates, has positive OLCI AFAI, and shows "
        "754-nm reflectance above blue 490-nm reflectance. The recovery "
        "concept is informed by Wang & Jiang (2025), but is not claimed to "
        "reproduce NOAA MSL12 thresholds exactly.",
    },
}

DEFAULT_MASK_PROFILE = "bloom_aware_cyan"

QA_EXCLUDE_FLAGS = ("HISATZEN", "NAVFAIL")

LAND_ADJACENCY_PIXELS = 0

NETCDF_DEFLATE = 4

ROW_BLOCK_SIZE = 256

DEFAULT_WATER_DEFINITION = "hydrolakes"

DEFAULT_SHORE_BUFFER_M = 0.0

DEFAULT_BRIGHT_SCREEN = False

DEFAULT_BRIGHT_RHOS865 = 0.08

DEFAULT_MAP_RESOLUTION_M = 300.0

MAP_MAX_CELLS = 60000000

INTERACTIVE_MAX_CELLS = 60000000

RAW_L1_DIR = Path("scratch/raw_l1_fallback").resolve()

STAGING_DIR = Path("scratch/staging").resolve()

L2_DIR = Path("scratch/l2").resolve()

DERIVED_DIR = Path("scratch/legacy_dense_derived").resolve()

EXPORT_DIR = Path("data/exports").resolve()

PAR_DIR = Path("data/par").resolve()

LOG_DIR = Path("data/logs").resolve()

HYDRO_ROOT = Path("data/hydrolakes").resolve()

CDSE_CATALOGUE_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"

CDSE_TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
)

CDSE_DOWNLOAD_ROOT = "https://download.dataspace.copernicus.eu/odata/v1"

CORE_RHOS_WAVELENGTHS = (620, 665, 681, 709)

SPECTRAL_RHOS_WAVELENGTHS = (
    400,
    412,
    443,
    490,
    510,
    560,
    620,
    665,
    674,
    681,
    709,
    754,
    779,
    865,
    885,
    900,
)

AQUATIC_RRS_WAVELENGTHS = (400, 412, 443, 490, 510, 560, 620, 665, 674, 681, 709)

TRUE_COLOR_PRODUCTS = ("rhot_443", "rhot_560", "rhot_665")

PRODUCT_PATTERN = re.compile("^(?P<prefix>[A-Za-z]+)_(?P<wavelength>\\d+(?:\\.\\d+)?)$", 32)

DERIVED_FILL = np.float32(-32767.0)

HYDROLAKES_URL = "https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip"

US_STATES_URL = "https://www2.census.gov/geo/tiger/GENZ2025/shp/cb_2025_us_state_5m.zip"

USE_GOOGLE_DRIVE_FOR_PERSISTENT = False

PERSISTENT_ROOT = Path("data").resolve()

SCRATCH_ROOT = Path("scratch").resolve()

COMPACT_DIR = Path("data/lake_only").resolve()

COMPACT_NC_DIR = Path("data/lake_only/netcdf").resolve()

STATS_SCENE_DIR = Path("data/lake_only/stats_by_scene").resolve()

CATALOG_DIR = Path("data/catalog").resolve()

HTML_DIR = Path("data/html_maps").resolve()

ADMIN_ROOT = Path("data/admin_boundaries").resolve()

REGISTRY_DB = Path("data/scene_registry.sqlite").resolve()

TARGET_LAKES_GPKG = Path("data/catalog/target_hydrolakes.gpkg").resolve()

DOWNLOAD_METHOD = "auto"

KEEP_STAGED_L1 = False

KEEP_RAW_ZIP_FALLBACK = False

KEEP_L2_INTERMEDIATE = False

KEEP_FAILED_SCRATCH = True

SAVE_TRUE_COLOR_L2 = False

SAVE_RHOS865_FOR_BRIGHT_SCREEN = False

SAVE_COMPACT_LAKE_NETCDF = True

SAVE_PIXEL_PARQUET = False

COMPACT_DEFLATE = 6

COMPACT_CHUNK_OBS = 65536

PROCESSING_STRATEGY = "adaptive_lake_windows"

L2_WINDOW_MARGIN_KM = 3.0

MAX_L2_WINDOWS_PER_SCENE = 6

FULL_SCENE_FRACTION_THRESHOLD = 0.35

REQUIRE_L2_NONLAND = True

CDSE_S3_ENDPOINT = "https://eodata.dataspace.copernicus.eu"

CDSE_S3_KEYS_URL = "https://s3-keys-manager.cloudferro.com/api/user/credentials"

GLOBAL_ANALYSIS_RHOS_WAVELENGTHS = (490, 560, 620, 665, 681, 709, 754, 865, 885)

TRUE_COLOR_RHOS_WAVELENGTHS = (665, 560, 490)

ADDITIONAL_INDEX_NAMES = ("MPH", "FAI")

OLCI_L2GEN_RHOS_NAME_OVERRIDES = {885: 884}

S3_FRAME_RE = re.compile(
    "^(S3[AB])_OL_1_EFR____(\\d{8}T\\d{6})_(\\d{8}T\\d{6})_(\\d{8}T\\d{6})_(\\d{4})_(\\d{3})_(\\d{3})_(\\d{4})_",
    32,
)

GEOB_API = "https://www.geoboundaries.org/api/current/gbOpen/{iso3}/ADM1/"

BASEMAPS = ["Esri Light Gray", "Esri Dark Gray", "Esri World Imagery", "OpenStreetMap", "None"]

SCALAR_METRICS = [
    "CI",
    "CI_cyano",
    "NDCI",
    "MPH",
    "FAI",
    "rhos_490",
    "rhos_560",
    "rhos_620",
    "rhos_665",
    "rhos_681",
    "rhos_709",
    "rhos_754",
    "rhos_865",
    "rhos_885",
]

MAP_METRICS = [
    "true_color",
    "CI",
    "CI_cyano",
    "NDCI",
    "MPH",
    "FAI",
    "rhos_490",
    "rhos_560",
    "rhos_620",
    "rhos_665",
    "rhos_681",
    "rhos_709",
    "rhos_754",
    "rhos_865",
    "rhos_885",
]

PLOT_METRICS = [
    "ci_cyano_detection_fraction_of_ci_valid",
    "ci_cyano_detection_fraction_of_retained",
    "CI_mean_candidates",
    "CI_median_candidates",
    "CIcyano_mean_detections",
    "CIcyano_median_detections",
    "NDCI_mean_retained",
    "NDCI_median_retained",
    "MPH_mean_retained",
    "FAI_mean_retained",
    "retained_percent",
]

CMAPS = [
    "viridis",
    "plasma",
    "inferno",
    "magma",
    "cividis",
    "turbo",
    "BrBG",
    "Spectral",
    "RdYlBu_r",
    "YlGnBu",
]

DISPLAY_MAX_CELLS = 8000000

SCIENCE_EXPORT_MAX_CELLS = 25000000

SCIENCE_EXPORT_DEFAULT_RESOLUTION_M = 300

SCIENCE_EXPORT_RHOS = [
    "rhos_490",
    "rhos_560",
    "rhos_620",
    "rhos_665",
    "rhos_681",
    "rhos_709",
    "rhos_754",
    "rhos_865",
    "rhos_885",
]

SCIENCE_EXPORT_INDICES = ["CI", "CI_cyano", "NDCI", "MPH", "FAI"]

SCIENCE_EXPORT_VARIABLES = [
    "rhos_490",
    "rhos_560",
    "rhos_620",
    "rhos_665",
    "rhos_681",
    "rhos_709",
    "rhos_754",
    "rhos_865",
    "rhos_885",
    "CI",
    "CI_cyano",
    "NDCI",
    "MPH",
    "FAI",
]

GRID_ID_RE = re.compile(
    "^G(?P<size>\\d{2})_(?P<ns>[NS])(?P<lat>\\d{2})_(?P<ew>[EW])(?P<lon>\\d{3})$", 34
)

ARCHIVE_INDEX_CSV = Path("data/catalog/global_archive_index.csv").resolve()

APP_CSS_GLOBAL = (
    ".gradio-container {max-width: 1550px !important}.small-note{font-size:.9rem;color:#4b5563}"
)

CYAN_STYLE_GRID_METHOD = "spherical_sinusoidal_sparse_native_observation_bin_v2"

CYAN_STYLE_SPHERE_RADIUS_M = 6371007.181

V264_SCHEMA = "global-efficient-v2.6.4"

CYAN_STRICT_QA_FLAGS = ("CLDICE", "HISATZEN", "NAVFAIL")

BLOOM_AWARE_HARD_FLAGS = ("HISATZEN", "NAVFAIL")

BLOOM_RESCUE_REQUIRE_CI_CYANO = True

BLOOM_RESCUE_REQUIRE_POSITIVE_AFAI = True

BLOOM_RESCUE_REQUIRE_NIR_OVER_BLUE = True

BLOOM_RESCUE_AFAI_MIN = 0.0

SCIENCE_EXPORT_CHUNK = 256

V265_SCHEMA = "global-efficient-v2.6.5"

V265_DISPLAY_MAX_CELLS = 30000000

V265_EXPORT_CHUNK = 128

V265_EXPORT_NEAREST_RADIUS_M = 600.0

V267_SCHEMA = "global-efficient-v2.6.7-efficient-production"

V267_GRID_METHOD = "epsg4326_same_lake_native_nearest_v2"

V267_DISPLAY_METHOD = "webmercator_same_lake_native_nearest_v2"

V267_EXPORT_CHUNK = 128

V267_NEAREST_RADIUS_M = 600.0

V268_SCHEMA = "global-efficient-v2.6.8-native-first-cyan-comparison"

L2GEN_BIN = None

GETANC_BIN = None

OCSSW_SUBPROCESS_ENV = {}

PYTHON_PROJ_DATA_DIR = None

PYTHON_PROJ_DATABASE = None

IS_COLAB = False
