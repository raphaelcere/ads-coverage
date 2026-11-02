"""
Theoretical reception coverage of a fixed ADS-B / ADS-L ground receiver.

Model
-----
Line-of-sight (LOS) visibility with an effective Earth radius k*R (k = 4/3 by default, the standard
atmosphere [2]), obstructed by terrain along radial profiles. For every azimuth the tool
decides, at each distance, whether an aircraft at a given level is visible from the antenna. Levels are
either constant altitudes (AMSL, the default) or heights above the local ground (AGL, ``--agl``).

Terrain
-------
- swissALTI3D (swisstopo, 2 m), downloaded automatically, for the near zone (``--near-km``).
- Copernicus GLO-30 (30 m, AWS open data) for the rest, read at low resolution from the cloud-optimised
  GeoTIFF overviews. Requires rasterio.
- or local GeoTIFF terrain models in EPSG:4326 (``--dem``).

Fading (optional, ``--fading``)
-------------------------------
On top of the geometry, the tool can estimate the probability that a message is received. The mean received
power comes from the free-space link budget (Friis transmission formula [5]) at the slant range. Around that mean
the envelope fluctuates as in Clarke's scattering model: the field is the sum of many randomly phased waves, so its
components are complex Gaussian and the envelope is Rayleigh distributed when there is no direct wave, and Rice
distributed when a direct wave of significant magnitude is added to the scattered ones [1, section 2.1,
pp. 961-962, which refers to Davenport and Root [4], pp. 160-161 (Rayleigh) and 165-167 (Rice), and to Rice [3]].
The power is therefore exponentially distributed (Rayleigh) or non-central chi-square distributed with two degrees
of freedom (Rice, K-factor = direct power / scattered power); these two closed forms are standard consequences of
the envelope distributions and are not written out in [1]. A message is received when the instantaneous power
exceeds the receiver sensitivity. Only these first-order statistics of the paper are used; its Doppler spectrum,
frequency correlation and time-delay results are not.

References
----------
[1] R. H. Clarke, "A Statistical Theory of Mobile-Radio Reception", Bell System Technical Journal, vol. 47,
    no. 6, pp. 957-1000, July-August 1968.
[2] Effective Earth radius and k-factor (k = 4/3): Recommendation ITU-R P.834, "Effects of tropospheric
    refraction on radiowave propagation", Annex 1, section 2 (Re = k a; equation (3) for k; k = 4/3 as the
    standard value). Read in version P.834-3.
[3] S. O. Rice, "Statistical Properties of a Sine Wave Plus Random Noise", Bell System Technical Journal,
    vol. 27, no. 1, pp. 109-157, January 1948 (cited by [1]).
[4] W. B. Davenport and W. L. Root, An Introduction to the Theory of Random Signals and Noise, McGraw-Hill,
    New York, 1958, pp. 160-161 and 165-167 (as cited by [1], p. 962).
[5] H. T. Friis, "A Note on a Simple Transmission Formula", Proceedings of the IRE, vol. 34, no. 5,
    pp. 254-256, May 1946.

Output: a KML file (one layer per level) and a summary table of the ranges on the console.

Examples
--------
    adsb-coverage 46.9480 7.4474 --ant-h 15 --alt-ft 1000,3000,10000,30000 -o bern.kml
    adsb-coverage 46.9480 7.4474 --agl -o low.kml        # low altitudes, heights above ground

Dependencies: numpy, requests; rasterio for the far field (Copernicus) and local GeoTIFFs; scipy for --fading rice.
"""
import argparse
import io
import math
import sys
import warnings
import zipfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from xml.sax.saxutils import escape

import numpy as np

R_EARTH = 6371008.8  # mean Earth radius (m)
FT = 0.3048  # metres per foot
C_LIGHT = 299_792_458.0  # speed of light (m/s)
STAC_ITEMS = "https://data.geo.admin.ch/api/stac/v0.9/collections/ch.swisstopo.swissalti3d/items"
COP_URL = "https://copernicus-dem-30m.s3.amazonaws.com/{n}/{n}.tif"


# ----------------------------------------------------------------------------- network
# Settings shared by all downloads; filled in from the command line. With proxies=None, requests falls back
# on the HTTPS_PROXY / HTTP_PROXY environment variables.
HTTP = {"proxies": None, "verify": True, "timeout": 30}


def http_get(url, **kwargs):
    """requests.get with the proxy, CA bundle and timeout from HTTP."""
    import requests

    kwargs.setdefault("timeout", HTTP["timeout"])
    return requests.get(url, proxies=HTTP["proxies"], verify=HTTP["verify"], **kwargs)


def network_help(exc):
    return (f"Connection failed: {exc}\n"
            "Options: (1) proxy: --proxy http://user:password@host:port or the HTTPS_PROXY variable; "
            "(2) SSL inspection: --ca-bundle path/to/ca.pem; "
            "(3) offline: download the swissALTI3D tiles (2 m, .xyz.zip) on a machine with access and pass "
            "their folder with --tile-dir; (4) --near-km 0 --dem terrain.tif.")


# ----------------------------------------------------------------------------- geodesy
def wgs84_to_lv95(lat, lon):
    """WGS84 -> Swiss LV95 (E, N) with the approximate swisstopo formulae (about 1 m in Switzerland)."""
    p = (np.asarray(lat, float) * 3600 - 169028.66) / 10000
    lam = (np.asarray(lon, float) * 3600 - 26782.5) / 10000
    e = 2600072.37 + 211455.93 * lam - 10938.51 * lam * p - 0.36 * lam * p**2 - 44.54 * lam**3
    n = 1200147.07 + 308807.95 * p + 3745.25 * lam**2 + 76.63 * p**2 - 194.56 * lam**2 * p + 119.79 * p**3
    return e, n


def destination(lat, lon, az_deg, dist_m):
    """Point reached from (lat, lon) along the given azimuth and distance on a sphere (numpy broadcasting)."""
    p1, l1 = math.radians(lat), math.radians(lon)
    th = np.radians(az_deg)
    d = dist_m / R_EARTH
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(th))
    l2 = l1 + np.arctan2(np.sin(th) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2))
    return np.degrees(p2), (np.degrees(l2) + 540.0) % 360.0 - 180.0


# ----------------------------------------------------------------------------- terrain: swissALTI3D
_CELLS = 500  # 2 m cells along one side of a 1 km tile
_CELL_M = 2.0
_FACTORS = (1, 2, 4, 5, 10, 20, 25, 50)  # decimation factors that divide a tile exactly
_KEY = 100_000  # tile key = E_km * _KEY + N_km (both below 100 000)


def decimate_xyz_zip(data, factor):
    """Read a swissALTI3D 2 m .xyz.zip tile and average it over factor x factor cells (north at row 0)."""
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        name = next(n for n in z.namelist() if n.lower().endswith(".xyz"))
        lines = z.read(name).decode("ascii", "ignore").splitlines()
    if lines and not lines[0].lstrip()[:1].isdigit():
        lines = lines[1:]  # "X Y Z" header
    xyz = np.array(" ".join(lines).split(), dtype=float).reshape(-1, 3)
    col = np.rint((xyz[:, 0] - xyz[:, 0].min()) / _CELL_M).astype(int)
    row = np.rint((xyz[:, 1].max() - xyz[:, 1]) / _CELL_M).astype(int)
    grid = np.full((_CELLS, _CELLS), np.nan)
    grid[np.clip(row, 0, _CELLS - 1), np.clip(col, 0, _CELLS - 1)] = xyz[:, 2]
    grid[grid < -500] = np.nan  # no-data
    n = _CELLS // factor
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN blocks
        return np.nanmean(grid.reshape(n, factor, n, factor), axis=(1, 3)).astype(np.float32)


class SwissAlti3D:
    """swissALTI3D 1 km tiles (2 m), decimated to about `res_m` and cached on disk as .npy files.

    Missing tiles are downloaded through the swisstopo STAC API, or, with `tile_dir`, read from a folder of
    .xyz.zip files downloaded by hand (offline mode, nothing is downloaded).
    """

    def __init__(self, lat, lon, radius_m, res_m=50, cache_dir=None, workers=8, tile_dir=None):
        self.factor = min(_FACTORS, key=lambda f: abs(f * _CELL_M - res_m))
        self.res = self.factor * _CELL_M
        self.n = _CELLS // self.factor
        self.cache = Path(cache_dir or Path.home() / ".cache" / "adsb_coverage" / "swissalti3d")
        self.cache.mkdir(parents=True, exist_ok=True)
        self.workers = workers
        self.tile_dir = Path(tile_dir) if tile_dir else None
        self.tiles = {}
        self._load(lat, lon, radius_m)

    @staticmethod
    def _tiles_around(lat, lon, radius_m):
        """(E_km, N_km) of the tiles that may lie within radius_m; this is the swisstopo tile id, e.g. 2485-1109."""
        e0, n0 = wgs84_to_lv95(lat, lon)
        half_diagonal = 710  # m, so that corner tiles are kept
        return [(tx, ty)
                for tx in range(int((e0 - radius_m) // 1000), int((e0 + radius_m) // 1000) + 1)
                for ty in range(int((n0 - radius_m) // 1000), int((n0 + radius_m) // 1000) + 1)
                if math.hypot(tx * 1000 + 500 - e0, ty * 1000 + 500 - n0) <= radius_m + half_diagonal]

    def _cache_path(self, tx, ty):
        return self.cache / f"{tx}-{ty}_{int(self.res)}m.npy"

    def _store(self, tx, ty, array):
        np.save(self._cache_path(tx, ty), array)
        self.tiles[tx * _KEY + ty] = array

    def _load(self, lat, lon, radius_m):
        wanted = self._tiles_around(lat, lon, radius_m)
        missing = []
        for tx, ty in wanted:
            path = self._cache_path(tx, ty)
            if path.exists():
                self.tiles[tx * _KEY + ty] = np.load(path)
            else:
                missing.append((tx, ty))
        if missing and self.tile_dir:
            missing = self._load_local(missing)
            if missing:
                print(f"swissALTI3D: {len(missing)} tile(s) not found in --tile-dir", file=sys.stderr)
        elif missing:
            self._download(lat, lon, radius_m, missing)
        print(f"swissALTI3D: {len(self.tiles)}/{len(wanted)} tiles available", file=sys.stderr)

    def _load_local(self, missing):
        """Read tiles from tile_dir; return those that were not found."""
        not_found = []
        for tx, ty in missing:
            path = next(iter(sorted(self.tile_dir.glob(f"*{tx}-{ty}_2_2056*.xyz.zip"))), None)
            if path is None:
                not_found.append((tx, ty))
            else:
                self._store(tx, ty, decimate_xyz_zip(path.read_bytes(), self.factor))
        return not_found

    def _download(self, lat, lon, radius_m, missing):
        print(f"swissALTI3D: downloading {len(missing)} tiles...", file=sys.stderr)
        dlat = radius_m / 111_200
        dlon = radius_m / (111_200 * math.cos(math.radians(lat)))
        try:
            items = self._list_items((lon - dlon, lat - dlat, lon + dlon, lat + dlat))
        except Exception as exc:  # noqa: BLE001
            raise SystemExit(network_help(exc)) from None
        todo = [((tx, ty), items[f"{tx}-{ty}"]) for tx, ty in missing if f"{tx}-{ty}" in items]
        with ThreadPoolExecutor(self.workers) as pool:
            for (tx, ty), array in pool.map(self._fetch, todo):
                if array is not None:
                    self._store(tx, ty, array)

    @staticmethod
    def _list_items(bbox):
        """Most recent 2 m tile of each id in the bounding box: {"2485-1109": (year, url)}."""
        url, params, tiles = STAC_ITEMS, {"bbox": ",".join(f"{v:.6f}" for v in bbox), "limit": 100}, {}
        while url:
            response = http_get(url, params=params, timeout=60)
            response.raise_for_status()
            page = response.json()
            for item in page.get("features", []):
                _, year, tile_id = item["id"].split("_")
                href = next((a["href"] for name, a in item["assets"].items()
                             if "_2_2056" in name and name.endswith(".xyz.zip")), None)
                if href and (tile_id not in tiles or year > tiles[tile_id][0]):
                    tiles[tile_id] = (year, href)
            url = next((link["href"] for link in page.get("links", []) if link.get("rel") == "next"), None)
            params = None  # the "next" link already carries the query
        return tiles

    def _fetch(self, tile):
        (tx, ty), (year, href) = tile
        for _ in range(3):
            try:
                response = http_get(href, timeout=120)
                response.raise_for_status()
                return (tx, ty), decimate_xyz_zip(response.content, self.factor)
            except Exception as exc:  # noqa: BLE001
                error = exc
        print(f"  ! tile {tx}-{ty} ({year}) failed: {error}", file=sys.stderr)
        return (tx, ty), None

    def sample(self, lat, lon):
        """Elevation (m) at the given points; NaN where no tile is loaded."""
        e, n = wgs84_to_lv95(lat, lon)
        out = np.full(e.shape, np.nan)
        keys = np.floor(e / 1000).astype(np.int64) * _KEY + np.floor(n / 1000).astype(np.int64)
        for key in np.unique(keys):
            tile = self.tiles.get(int(key))
            if tile is None:
                continue
            tx, ty = divmod(int(key), _KEY)
            m = keys == key
            col = np.clip(((e[m] - tx * 1000) / self.res).astype(int), 0, self.n - 1)
            row = np.clip(((ty * 1000 + 1000 - n[m]) / self.res).astype(int), 0, self.n - 1)
            out[m] = tile[row, col]
        return out


# ----------------------------------------------------------------------------- terrain: GeoTIFF (WGS84)
class GeoTiffMosaic:
    """GeoTIFF terrain models in EPSG:4326 (local files or URLs), read at low resolution.

    Cloud-optimised GeoTIFFs are read through their overviews, so a 400 km radius stays fast. Each block is
    reduced with Resampling.max (conservative for obstruction), or with the average if rasterio refuses.
    """

    def __init__(self, paths, bbox, res_m=300):
        import rasterio
        from rasterio.enums import Resampling
        from rasterio.windows import Window

        lon0, lat0, lon1, lat1 = bbox
        self.chunks = []  # (west, north, cell width, cell height, array)
        unavailable = 0
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif"):
            for path in paths:
                try:
                    ds = rasterio.open(path)
                except Exception:  # noqa: BLE001  (missing tile: sea, outside coverage, network error)
                    unavailable += 1
                    continue
                with ds:
                    if not (ds.crs and ds.crs.is_geographic):
                        raise SystemExit(f"{path}: the terrain model must use geographic coordinates (EPSG:4326).")
                    t = ds.transform
                    west, east = max(lon0, ds.bounds.left), min(lon1, ds.bounds.right)
                    south, north = max(lat0, ds.bounds.bottom), min(lat1, ds.bounds.top)
                    if west >= east or south >= north:
                        continue
                    c0, r0 = (~t) * (west, north)
                    c1, r1 = (~t) * (east, south)
                    c0, r0 = max(int(math.floor(c0)), 0), max(int(math.floor(r0)), 0)
                    c1, r1 = min(int(math.ceil(c1)), ds.width), min(int(math.ceil(r1)), ds.height)
                    step_deg = res_m / 111_320
                    out_w = max(1, int(round((c1 - c0) * abs(t.a) / step_deg)))
                    out_h = max(1, int(round((r1 - r0) * abs(t.e) / step_deg)))
                    window = Window(c0, r0, c1 - c0, r1 - r0)
                    try:
                        array = ds.read(1, window=window, out_shape=(out_h, out_w), resampling=Resampling.max)
                    except Exception:  # noqa: BLE001
                        array = ds.read(1, window=window, out_shape=(out_h, out_w), resampling=Resampling.average)
                    array = array.astype(np.float32)
                    if ds.nodata is not None:
                        array[array == ds.nodata] = np.nan
                    self.chunks.append((t.c + c0 * t.a, t.f + r0 * t.e, (c1 - c0) * t.a / out_w,
                                        (r1 - r0) * abs(t.e) / out_h, array))
        print(f"GeoTIFF terrain: {len(self.chunks)} tile(s) read, {unavailable} unavailable", file=sys.stderr)

    def sample(self, lat, lon):
        """Elevation (m) at the given points; NaN outside the tiles read."""
        out = np.full(lat.shape, np.nan)
        for west, north, dx, dy, array in self.chunks:
            col = np.floor((lon - west) / dx).astype(int)
            row = np.floor((north - lat) / dy).astype(int)
            m = (col >= 0) & (col < array.shape[1]) & (row >= 0) & (row < array.shape[0]) & np.isnan(out)
            out[m] = array[row[m], col[m]]
        return out


def copernicus_urls(bbox):
    """URLs of the Copernicus GLO-30 1-degree tiles that cover the bounding box (lon0, lat0, lon1, lat1)."""
    lon0, lat0, lon1, lat1 = bbox
    return [COP_URL.format(n=f"Copernicus_DSM_COG_10_{'N' if la >= 0 else 'S'}{abs(la):02d}_00_"
                             f"{'E' if lo >= 0 else 'W'}{abs(lo):03d}_00_DEM")
            for la in range(math.floor(lat0), math.floor(lat1) + 1)
            for lo in range(math.floor(lon0), math.floor(lon1) + 1)]


class Terrain:
    """Chain of terrain sources: the first one that has a value wins; unknown points count as 0 m (sea)."""

    def __init__(self, sources):
        self.sources = list(sources)
        self.unknown = 0.0  # share of the last sampled points that no source knew

    def sample(self, lat, lon):
        lat, lon = np.asarray(lat), np.asarray(lon)
        z = np.full(lat.shape, np.nan)
        for source in self.sources:
            m = np.isnan(z)
            if not m.any():
                break
            z[m] = source.sample(lat[m], lon[m])
        self.unknown = float(np.isnan(z).mean())
        return np.nan_to_num(z, nan=0.0)


# ----------------------------------------------------------------------------- propagation model
def _terrain_profiles(terrain, lat, lon, ant_h, rmax, step, az_step, k, ground):
    """Terrain along radial profiles and the steepest slope to it, seen from the antenna."""
    az = np.arange(0, 360, az_step, dtype=float)
    d = np.arange(step, rmax + step, step, dtype=float)
    points_lat, points_lon = destination(lat, lon, az[:, None], d[None, :])
    ground_m = float(terrain.sample(np.array([lat]), np.array([lon]))[0]) if ground is None else ground
    z = terrain.sample(points_lat, points_lon)
    antenna = ground_m + ant_h
    bulge = d**2 / (2 * k * R_EARTH)
    steepest_terrain = np.maximum.accumulate((z - bulge - antenna) / d, axis=1)
    return az, d, z, ground_m, antenna, bulge, steepest_terrain


def _level_visibility(z, height, agl, antenna, bulge, d, steepest_terrain):
    """Altitude AMSL of the aircraft (scalar, or [azimuth, distance] in AGL) and its visibility mask."""
    aircraft = z + height if agl else height
    return aircraft, (aircraft - bulge - antenna) / d >= steepest_terrain


def compute_coverage(terrain, lat, lon, ant_h, heights_m, rmax, step, az_step, k, ground=None, agl=False):
    """
    Visibility of aircraft along radial profiles around the antenna.

    Returns (azimuths, distances, visible, ground_m) where visible[level, azimuth, distance] is True when an
    aircraft at that level and distance is in line of sight of the antenna.

    The effective Earth radius is k R [2]. In the "flat Earth" frame (the usual parabolic approximation of the Earth
    bulge) a point at distance d and height h has the apparent height h - d^2 / (2 k R).
    An aircraft is visible when the slope from the antenna to it is at least the steepest slope to any terrain
    met so far (including the terrain under the aircraft).
    - AMSL levels (agl=False): constant altitudes. Visibility ends at the first obstruction.
    - AGL levels (agl=True): heights above the local ground, the aircraft follows the relief. Visibility is
      no longer monotonic with distance: an aircraft can reappear beyond a ridge, over a higher valley.
    """
    az, d, z, ground_m, antenna, bulge, steepest = _terrain_profiles(terrain, lat, lon, ant_h, rmax, step,
                                                                     az_step, k, ground)
    visible = np.stack([_level_visibility(z, h, agl, antenna, bulge, d, steepest)[1] for h in heights_m])
    return az, d, visible, ground_m


@dataclass(frozen=True)
class LinkBudget:
    """Free-space link budget and fading statistics of the transmission.

    k_db is the Rice K-factor in dB (direct power / scattered power); None means Rayleigh, i.e. Clarke's case of a
    purely scattered field.
    """

    freq_mhz: float
    eirp_dbm: float  # equivalent isotropic radiated power of the aircraft transmitter
    rx_gain_dbi: float  # receiving antenna gain, assumed constant in elevation and azimuth
    sensitivity_dbm: float  # minimum received power for a message to be decoded
    losses_db: float = 0.0  # cable, connector and other losses at the receiver
    k_db: Optional[float] = None

    def mean_power_dbm(self, slant_m):
        """Mean received power (dBm) at the slant range: EIRP + gain - losses - free-space path loss [5]."""
        path_loss = 20 * np.log10(4 * np.pi * slant_m * self.freq_mhz * 1e6 / C_LIGHT)
        return self.eirp_dbm + self.rx_gain_dbi - self.losses_db - path_loss

    def reception_probability(self, slant_m):
        """Probability that the instantaneous power exceeds the sensitivity at the slant range."""
        return reception_probability(self.mean_power_dbm(slant_m) - self.sensitivity_dbm, self.k_db)


def reception_probability(margin_db, k_db=None):
    """P(instantaneous power >= sensitivity) for a mean power `margin_db` above the sensitivity.

    With s = sensitivity / mean power (linear), the power of a Rayleigh envelope (scattered field only, Clarke [1],
    section 2.1, p. 962) is exponentially distributed: P = exp(-s).

    With a direct wave (Rice [3, 4]) and K-factor K (linear), 2(K+1) times the power over its mean is non-central
    chi-square with 2 degrees of freedom and non-centrality 2K, so P = ncx2.sf(2(K+1) s, 2, 2K), the Marcum
    Q-function Q1(sqrt(2K), sqrt(2(K+1) s)). K = 0 (linear, i.e. k_db = -inf) gives the Rayleigh result.
    Rice requires scipy.
    """
    s = 10.0 ** (-np.asarray(margin_db, dtype=float) / 10)
    if k_db is None:
        return np.exp(-s)
    from scipy.stats import ncx2

    k = 10.0 ** (k_db / 10)
    return ncx2.sf(2 * (k + 1) * s, 2, 2 * k)


def slant_range(antenna, aircraft, d):
    """Straight-line distance (m) between the antenna and an aircraft at ground distance d (arc length)."""
    r1, r2 = R_EARTH + antenna, R_EARTH + aircraft
    return np.sqrt((r2 - r1) ** 2 + 4 * r1 * r2 * np.sin(d / (2 * R_EARTH)) ** 2)


def compute_reception(terrain, lat, lon, ant_h, heights_m, rmax, step, az_step, k, link, ground=None, agl=False):
    """
    Probability of receiving a message, along the same profiles as compute_coverage.

    Returns (azimuths, distances, probability, ground_m); probability[level, azimuth, distance] is 0 where the
    aircraft is not in line of sight, otherwise the reception probability given by the link budget `link`.
    """
    az, d, z, ground_m, antenna, bulge, steepest = _terrain_profiles(terrain, lat, lon, ant_h, rmax, step,
                                                                     az_step, k, ground)
    levels = []
    for h in heights_m:
        aircraft, visible = _level_visibility(z, h, agl, antenna, bulge, d, steepest)
        probability = link.reception_probability(slant_range(antenna, aircraft, d))
        levels.append(np.where(visible, probability, 0.0))
    return az, d, np.stack(levels), ground_m


def continuous_range(visible, d):
    """Distance to the first obstruction along each azimuth: the continuous coverage from the receiver."""
    blocked = ~visible
    first = np.where(blocked.any(axis=-1), blocked.argmax(axis=-1), visible.shape[-1])
    return np.where(first > 0, d[np.clip(first - 1, 0, len(d) - 1)], 0.0)


def max_visible_range(visible, d):
    """Distance of the farthest visible point along each azimuth (0 if none), even behind a ridge."""
    last = visible.shape[-1] - 1 - visible[..., ::-1].argmax(axis=-1)
    return np.where(visible.any(axis=-1), d[last], 0.0)


def visible_runs(visible, d, step, min_len=2):
    """Contiguous visible stretches for one level (array [azimuth, distance]).

    Returns (azimuth indices, start distances, end distances) in metres, for stretches of at least
    min_len samples.
    """
    pad = np.zeros((visible.shape[0], 1), dtype=int)
    edges = np.diff(np.concatenate([pad, visible.astype(int), pad], axis=1), axis=1)
    (rows, starts), (_, ends) = np.nonzero(edges == 1), np.nonzero(edges == -1)  # row-major: starts and ends pair up
    keep = (ends - starts) >= min_len
    rows, starts, ends = rows[keep], starts[keep], ends[keep]
    return rows, d[starts] - step / 2, d[ends - 1] + step / 2


# ----------------------------------------------------------------------------- KML
def _polygon(ring, tessellate=False, height_m=None):
    """KML polygon; ring is a list of (lon, lat) and is closed automatically.

    Clamped to the ground, or, with height_m, drawn that many metres above the ground (relativeToGround), which
    keeps it visible above the 3D buildings of Google Earth.
    """
    coords = " ".join(f"{x:.6f},{y:.6f},{0 if height_m is None else height_m:.1f}" for x, y in list(ring) + [ring[0]])
    mode = "clampToGround" if height_m is None else "relativeToGround"
    return (f"<Polygon>{'<tessellate>1</tessellate>' if tessellate else ''}<altitudeMode>{mode}</altitudeMode>"
            f"<outerBoundaryIs><LinearRing><coordinates>{coords}</coordinates></LinearRing></outerBoundaryIs></Polygon>")


def _level_colours(n):
    """n colours (r, g, b) from red (first level) through green to blue (last level)."""
    for i in range(n):
        t = i / max(n - 1, 1)
        yield int(255 * t), int(200 * (1 - abs(2 * t - 1)) + 55), int(255 * (1 - t))


def write_kml(path, lat, lon, ground_m, ant_h, az, d, visible, levels_ft, k, agl=False, az_step=1.0, step=250.0,
              quantity="visible", note="", clamp=False):
    """Write one KML layer per level; return (continuous range, farthest range), each [level, azimuth] in m.

    `visible` is the boolean mask [level, azimuth, distance] to draw; `quantity` names it in the labels
    ("visible", or "received" when it comes from a reception probability) and `note` is added to the description.

    AMSL levels are drawn as the outline of the continuous coverage; AGL levels as the visible sectors,
    which may be fragmented in mountainous terrain. AGL sectors are drawn at their real height above the ground
    (the aircraft height), unless `clamp` is True, which drapes them on the ground.
    """
    reference = "AGL" if agl else "AMSL"
    continuous = np.stack([continuous_range(v, d) for v in visible])
    farthest = np.stack([max_visible_range(v, d) for v in visible])
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
           f"<name>ADS-B/ADS-L coverage {lat:.4f}, {lon:.4f}</name>",
           f"<description>Theoretical line-of-sight coverage, k={k:.4f}, antenna {ant_h:g} m AGL "
           f"(ground {ground_m:.0f} m). Levels in {reference}"
           f"{' (height above ground)' if agl else ''}.{' ' + note if note else ''}</description>"]
    for i, (r, g, b) in enumerate(_level_colours(len(levels_ft))):
        out.append(f'<Style id="s{i}"><LineStyle><color>ff{b:02x}{g:02x}{r:02x}</color><width>2</width></LineStyle>'
                   f'<PolyStyle><color>{"55" if agl else "33"}{b:02x}{g:02x}{r:02x}</color>'
                   f'<outline>{0 if agl else 1}</outline></PolyStyle></Style>')
    out.append('<Style id="rx"><IconStyle><scale>1.2</scale></IconStyle></Style>')
    out.append(f'<Placemark><name>Receiver</name><styleUrl>#rx</styleUrl><description>Ground {ground_m:.0f} m, '
               f'antenna +{ant_h:g} m</description><Point><coordinates>{lon:.6f},{lat:.6f},0</coordinates></Point></Placemark>')
    out.append("<Folder><name>Coverage by level</name>")
    for i in reversed(range(len(levels_ft))):  # high levels underneath, low levels on top
        label = (f"{levels_ft[i]:g} ft {reference} - continuous max {continuous[i].max() / 1000:.0f} km, "
                 f"{quantity} max {farthest[i].max() / 1000:.0f} km")
        if agl:
            rows, starts, ends = visible_runs(visible[i], d, step)
            if len(rows) == 0:
                continue
            half = az_step / 2
            sectors = []
            for row, r0, r1 in zip(rows, starts, ends):
                a0, a1 = az[row] - half, az[row] + half
                corners = [destination(lat, lon, a, r) for a, r in ((a0, r0), (a1, r0), (a1, r1), (a0, r1))]
                sectors.append(_polygon([(float(c[1]), float(c[0])) for c in corners],
                                        height_m=None if clamp else levels_ft[i] * FT))
            geometry = f"<MultiGeometry>{''.join(sectors)}</MultiGeometry>"
            description = f"<description>{visible[i].mean():.1%} of sampled points {quantity}</description>"
        else:
            reach = continuous[i]
            plat, plon = destination(lat, lon, az, reach)
            plat, plon = np.where(reach > 0, plat, lat), np.where(reach > 0, plon, lon)
            geometry = _polygon(list(zip(plon, plat)), tessellate=True)
            description = ""
        out.append(f"<Placemark><name>{escape(label)}</name>{description}<styleUrl>#s{i}</styleUrl>{geometry}</Placemark>")
    out.append("</Folder></Document></kml>")
    Path(path).write_text("\n".join(out), encoding="utf-8")
    return continuous, farthest


# ----------------------------------------------------------------------------- command line
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lat", type=float, help="latitude, WGS84 decimal degrees")
    ap.add_argument("lon", type=float, help="longitude, WGS84 decimal degrees")
    ap.add_argument("-o", "--output", default="coverage.kml", help="output KML file")
    ap.add_argument("--ant-h", type=float, default=10.0, help="antenna height above ground (m)")
    ap.add_argument("--ground-m", type=float, default=None, help="ground elevation at the receiver (default: read from the terrain model)")
    ap.add_argument("--agl", action="store_true",
                    help="levels are heights above the local ground (low altitudes) instead of altitudes AMSL")
    ap.add_argument("--alt-ft", default=None,
                    help="comma-separated levels in ft (default AMSL: 1000..30000; default AGL: 100..3000)")
    ap.add_argument("--k", type=float, default=4 / 3, help="effective Earth radius factor (default 4/3)")
    ap.add_argument("--rmax-km", type=float, default=None, help="maximum range computed (default 400 AMSL, 150 AGL)")
    ap.add_argument("--step-m", type=float, default=None, help="sampling step along each profile (default 250 m AMSL, 100 m AGL)")
    ap.add_argument("--az-step", type=float, default=1.0, help="azimuth step (degrees)")
    ap.add_argument("--near-km", type=float, default=10.0, help="radius of the swissALTI3D zone (0 = disabled)")
    ap.add_argument("--near-res-m", type=float, default=50.0, help="swissALTI3D resolution after decimation (m); 20 = finer")
    ap.add_argument("--far", choices=["copernicus", "none"], default="copernicus", help="terrain beyond the near zone")
    ap.add_argument("--far-res-m", type=float, default=300.0, help="resolution of the far terrain (m)")
    ap.add_argument("--dem", nargs="*", default=None, help="local EPSG:4326 GeoTIFF(s), used instead of Copernicus")
    ap.add_argument("--cache-dir", default=None, help="cache folder for swissALTI3D tiles")
    ap.add_argument("--tile-dir", default=None, help="folder of swissALTI3D .xyz.zip tiles downloaded by hand (offline)")
    ap.add_argument("--proxy", default=None, help="HTTP(S) proxy, e.g. http://host:8080")
    ap.add_argument("--ca-bundle", default=None, help="CA certificate(s) in PEM format (proxy with SSL inspection)")
    ap.add_argument("--timeout", type=float, default=30.0, help="connection timeout (s)")
    ap.add_argument("--clamp-to-ground", action="store_true",
                    help="KML: with --agl, drape the zones on the ground instead of drawing them at their real height")
    ap.add_argument("--flat", action="store_true", help="flat terrain at 0 m (theoretical radio horizon, for testing)")

    fading = ap.add_argument_group("fading and link budget (Clarke 1968)",
                                   "With --fading, the drawn zone is where a message is received with at least the "
                                   "probability --prob, instead of the pure line-of-sight zone. All the link budget "
                                   "values below must then be given.")
    fading.add_argument("--fading", choices=["none", "rayleigh", "rice"], default="none",
                        help="rayleigh: scattered field only (Clarke); rice: direct wave plus scattered field, needs --k-db and scipy")
    fading.add_argument("--freq-mhz", type=float, default=None, help="carrier frequency (MHz)")
    fading.add_argument("--eirp-dbm", type=float, default=None, help="aircraft transmitter EIRP (dBm)")
    fading.add_argument("--rx-gain-dbi", type=float, default=None, help="receiving antenna gain (dBi), constant in all directions")
    fading.add_argument("--sensitivity-dbm", type=float, default=None, help="receiver sensitivity, minimum power to decode a message (dBm)")
    fading.add_argument("--losses-db", type=float, default=0.0, help="cable and connector losses at the receiver (dB)")
    fading.add_argument("--k-db", type=float, default=None, help="Rice K-factor, direct over scattered power (dB); for the purely scattered case use --fading rayleigh")
    fading.add_argument("--prob", type=float, default=0.9, help="minimum reception probability of the drawn zone, between 0 and 1 (default 0.9)")
    return ap


def link_from_args(ap, a):
    """LinkBudget from the command line, or None without --fading; reports missing or invalid values with ap.error."""
    if a.fading == "none":
        return None
    required = [("--freq-mhz", a.freq_mhz), ("--eirp-dbm", a.eirp_dbm), ("--rx-gain-dbi", a.rx_gain_dbi),
                ("--sensitivity-dbm", a.sensitivity_dbm)] + ([("--k-db", a.k_db)] if a.fading == "rice" else [])
    missing = [name for name, value in required if value is None]
    if missing:
        ap.error(f"--fading {a.fading} requires: {', '.join(missing)}")
    if not 0 < a.prob < 1:
        ap.error("--prob must be between 0 and 1 (exclusive)")
    if a.fading == "rice":
        try:
            import scipy.stats  # noqa: F401
        except ImportError:
            ap.error("--fading rice requires scipy: pip install scipy")
    return LinkBudget(a.freq_mhz, a.eirp_dbm, a.rx_gain_dbi, a.sensitivity_dbm, a.losses_db,
                      a.k_db if a.fading == "rice" else None)


def main(argv=None):
    ap = build_parser()
    a = ap.parse_args(argv)
    link = link_from_args(ap, a)

    if a.proxy:
        HTTP["proxies"] = {"http": a.proxy, "https": a.proxy}
    if a.ca_bundle:
        HTTP["verify"] = a.ca_bundle
    HTTP["timeout"] = a.timeout

    levels_ft = [float(x) for x in (a.alt_ft or ("100,300,500,1000,2000,3000" if a.agl else "1000,3000,6000,10000,20000,30000")).split(",")]
    rmax = (a.rmax_km or (150.0 if a.agl else 400.0)) * 1000
    step = a.step_m or (100.0 if a.agl else 250.0)
    dlat = rmax / 111_200
    dlon = rmax / (111_200 * math.cos(math.radians(a.lat)))
    bbox = (a.lon - dlon, a.lat - dlat, a.lon + dlon, a.lat + dlat)

    sources = []
    if not a.flat:
        if a.near_km > 0:
            sources.append(SwissAlti3D(a.lat, a.lon, a.near_km * 1000, a.near_res_m, a.cache_dir, tile_dir=a.tile_dir))
        if a.dem or a.far == "copernicus":
            try:
                sources.append(GeoTiffMosaic(a.dem or copernicus_urls(bbox), bbox, a.far_res_m))
            except ImportError:
                raise SystemExit("rasterio is required for the far terrain: pip install rasterio "
                                 "(or use --far none to rely on swissALTI3D only).") from None
    terrain = Terrain(sources)

    heights_m = [f * FT for f in levels_ft]
    if link:
        az, d, probability, ground_m = compute_reception(terrain, a.lat, a.lon, a.ant_h, heights_m, rmax, step,
                                                         a.az_step, a.k, link, a.ground_m, a.agl)
        visible = probability >= a.prob
        quantity = "received"
        statistic = f"Rice K={link.k_db:g} dB" if link.k_db is not None else "Rayleigh"
        note = (f"Reception probability at least {a.prob:.0%} ({statistic}; {link.freq_mhz:g} MHz, EIRP {link.eirp_dbm:g} dBm, "
                f"receiver gain {link.rx_gain_dbi:g} dBi, losses {link.losses_db:g} dB, sensitivity {link.sensitivity_dbm:g} dBm).")
    else:
        az, d, visible, ground_m = compute_coverage(terrain, a.lat, a.lon, a.ant_h, heights_m, rmax, step, a.az_step,
                                                    a.k, a.ground_m, a.agl)
        quantity, note = "visible", ""
    if terrain.unknown > 0.02 and not a.flat:
        print(f"Warning: {terrain.unknown:.0%} of the points have no terrain data (treated as 0 m).", file=sys.stderr)

    continuous, farthest = write_kml(a.output, a.lat, a.lon, ground_m, a.ant_h, az, d, visible, levels_ft, a.k,
                                     a.agl, a.az_step, step, quantity, note, a.clamp_to_ground)
    reference = "AGL" if a.agl else "AMSL"
    print(f"Receiver ground: {ground_m:.0f} m; antenna: {ground_m + a.ant_h:.0f} m AMSL; k = {a.k:.4f}; levels in {reference}")
    if note:
        print(note)
    print(f"{'Level':>9} | continuous range (km): {'min':>5} {'mean':>5} {'max':>5} | {'max ' + quantity:>12} | {quantity:>8}")
    for ft, cont, far, vis in zip(levels_ft, continuous, farthest, visible):
        print(f"{ft:>6g} ft | {'':>22} {cont.min() / 1000:5.1f} {cont.mean() / 1000:5.1f} {cont.max() / 1000:5.1f} | "
              f"{far.max() / 1000:9.1f} km | {vis.mean():8.1%}")
    print(f"KML: {a.output}")


if __name__ == "__main__":
    main()
