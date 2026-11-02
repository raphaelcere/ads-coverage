"""Theoretical reception coverage of a fixed ADS-B / ADS-L ground receiver (line of sight, terrain, k-factor)."""
from .core import (  # noqa: F401
    GeoTiffMosaic, LinkBudget, SwissAlti3D, SwissSurface3D, Terrain, compute_coverage, compute_reception, continuous_range,
    copernicus_urls, decimate_xyz_zip, destination, main, max_visible_range, reception_probability, slant_range,
    visible_runs, wgs84_to_lv95, write_kml,
)

__version__ = "0.1.0"
