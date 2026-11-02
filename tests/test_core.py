import io
import math
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pytest

import adsb_coverage as ac
from adsb_coverage.core import FT, R_EARTH

K_EARTH = 4 / 3 * R_EARTH  # effective Earth radius (m)
LAT, LON = 46.9, 7.4
KML = "{http://www.opengis.net/kml/2.2}"


class Flat:
    def sample(self, lat, lon):
        return np.zeros(lat.shape)


class EastProfile:
    """Terrain that depends only on the distance east of (LAT, LON): height_at(distance in m)."""

    def __init__(self, height_at):
        self.height_at = height_at

    def sample(self, lat, lon):
        dx = (lon - LON) * 111320 * math.cos(math.radians(LAT))
        return np.where(np.abs(lat - LAT) < 0.05, self.height_at(dx), 0.0)


def continuous(terrain, heights_m, rmax, step, az_step, agl=False, ant_h=10.0):
    az, d, visible, _ = ac.compute_coverage(ac.Terrain([terrain]), LAT, LON, ant_h, heights_m, rmax, step, az_step,
                                            4 / 3, agl=agl)
    return az, np.stack([ac.continuous_range(v, d) for v in visible])


def radio_horizon(antenna_m, aircraft_m):
    return math.sqrt(2 * K_EARTH * antenna_m) + math.sqrt(2 * K_EARTH * aircraft_m)


def test_lv95_bern_reference_point():
    e, n = ac.wgs84_to_lv95(46.9510827861, 7.4386324167)  # old Bern observatory
    assert abs(e - 2600000) < 2 and abs(n - 1200000) < 2


@pytest.mark.parametrize("altitude_m", [305.0, 914.0, 9144.0])
def test_flat_terrain_amsl_matches_radio_horizon(altitude_m):
    _, ranges = continuous(Flat(), [altitude_m], 450e3, 100.0, 10.0)
    assert abs(ranges[0].mean() - radio_horizon(10, altitude_m)) < 250


@pytest.mark.parametrize("height_ft", [100, 300, 1000])
def test_flat_terrain_agl_matches_radio_horizon(height_ft):
    _, ranges = continuous(Flat(), [height_ft * FT], 150e3, 100.0, 10.0, agl=True)
    assert abs(ranges[0].mean() - radio_horizon(10, height_ft * FT)) < 250


def test_ridge_shadows_the_east():
    ridge = EastProfile(lambda dx: 2000 * np.exp(-((dx - 20000) ** 2) / (2 * 1500**2)))  # 2000 m ridge, 20 km east
    az, ranges = continuous(ridge, [3000.0], 400e3, 100.0, 10.0)
    assert ranges[0][list(az).index(90)] < 40e3 < 150e3 < ranges[0][list(az).index(270)]


def test_agl_aircraft_reappears_beyond_a_ridge():
    """Low aircraft: visible, hidden behind a 300 m wall at 5 km, visible again over a 600 m plateau from 8 km."""
    terrain = EastProfile(lambda dx: np.where((dx > 5000) & (dx < 5500), 300.0, np.where(dx > 8000, 600.0, 0.0)))
    az, d, visible, _ = ac.compute_coverage(ac.Terrain([terrain]), LAT, LON, 10, [30.0], 15e3, 100.0, 90.0, 4 / 3, agl=True)
    east = list(az).index(90)
    rows, starts, ends = ac.visible_runs(visible[0], d, 100.0)
    assert (rows == east).sum() >= 2  # at least two separate visible stretches
    # the plateau edge hides the aircraft again a few hundred metres further, so the gain is modest
    assert ac.max_visible_range(visible, d)[0, east] > ac.continuous_range(visible, d)[0, east] + 2000


def test_visible_runs():
    visible = np.array([[0, 1, 1, 1, 0, 1, 0, 1, 1]], dtype=bool)
    d = np.arange(1, 10) * 100.0
    rows, starts, ends = ac.visible_runs(visible, d, 100.0, min_len=2)
    assert list(rows) == [0, 0]
    assert list(starts) == [150.0, 750.0] and list(ends) == [450.0, 950.0]


def test_terrain_chain_falls_back_to_the_next_source():
    class Half:  # knows only the northern points
        def sample(self, lat, lon):
            return np.where(lat > 0, 100.0, np.nan)

    class Fixed:
        def sample(self, lat, lon):
            return np.full(lat.shape, 7.0)

    lat, lon = np.array([1.0, -1.0]), np.array([0.0, 0.0])
    assert list(ac.Terrain([Half(), Fixed()]).sample(lat, lon)) == [100.0, 7.0]
    chain = ac.Terrain([Half()])
    assert list(chain.sample(lat, lon)) == [100.0, 0.0] and chain.unknown == 0.5


def make_swissalti3d_tile(tx_km, ty_km):
    """A fake 1 km swissALTI3D 2 m .xyz.zip tile whose elevation is 500 m + 0.1 m per metre eastwards."""
    rows = ["X Y Z"] + [f"{tx_km * 1000 + i * 2 + 1}.0 {ty_km * 1000 + 999 - j * 2}.0 {500 + i * 0.2:.2f}"
                        for j in range(500) for i in range(500)]
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as z:
        z.writestr("tile.xyz", "\n".join(rows))
    return buffer.getvalue()


def test_swissalti3d_offline_tile_and_cache(tmp_path):
    tiles = tmp_path / "tiles"
    tiles.mkdir()
    (tiles / "swissalti3d_2019_2600-1200_2_2056_5728.xyz.zip").write_bytes(make_swissalti3d_tile(2600, 1200))
    lat, lon = np.array([46.9555]), np.array([7.4452])  # about E 2600500, N 1200500: elevation about 550 m

    swiss = ac.SwissAlti3D(lat[0], lon[0], 300, res_m=50, cache_dir=tmp_path / "cache", tile_dir=tiles)
    assert 520 < swiss.sample(lat, lon)[0] < 580
    assert np.isnan(swiss.sample(np.array([47.5]), np.array([8.5]))[0])  # outside the loaded tiles

    cached = ac.SwissAlti3D(lat[0], lon[0], 300, res_m=50, cache_dir=tmp_path / "cache")  # no tile_dir, no network
    assert cached.sample(lat, lon)[0] == swiss.sample(lat, lon)[0]


def test_geotiff_mosaic_reads_a_local_model(tmp_path):
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    path = tmp_path / "dem.tif"
    with rasterio.open(path, "w", driver="GTiff", height=60, width=60, count=1, dtype="float32", crs="EPSG:4326",
                       transform=from_origin(7.0, 46.06, 0.001, 0.001), nodata=-9999) as ds:
        ds.write(np.full((60, 60), 800.0, dtype="float32"), 1)

    mosaic = ac.GeoTiffMosaic([str(path), str(tmp_path / "missing.tif")], (7.0, 46.0, 7.06, 46.06), res_m=300)
    z = mosaic.sample(np.array([46.03, 47.0]), np.array([7.03, 8.0]))
    assert z[0] == 800.0 and np.isnan(z[1])


def test_copernicus_tile_urls():
    urls = ac.copernicus_urls((-0.5, -0.5, 0.5, 0.5))
    assert len(urls) == 4
    assert any("S01_00_W001_00" in u for u in urls) and any("N00_00_E000_00" in u for u in urls)


@pytest.mark.parametrize("flags", [[], ["--agl"]])
def test_cli_writes_a_valid_kml(tmp_path, flags, capsys):
    out = tmp_path / "coverage.kml"
    ac.main([str(LAT), str(LON), "--flat", "--az-step", "10", "-o", str(out), *flags])
    root = ET.parse(out).getroot()
    assert len(root.findall(f".//{KML}Polygon")) > 0
    assert "KML:" in capsys.readouterr().out


# ----------------------------------------------------------------------------- fading (Clarke 1968)
def test_free_space_loss_at_1090_mhz_and_100_km():
    link = ac.LinkBudget(freq_mhz=1090, eirp_dbm=0, rx_gain_dbi=0, sensitivity_dbm=-84)
    assert abs(link.mean_power_dbm(100e3) - (-133.2)) < 0.1


def test_link_budget_adds_gain_and_subtracts_losses():
    base = ac.LinkBudget(1090, 51, 0, -84).mean_power_dbm(50e3)
    assert ac.LinkBudget(1090, 51, 3, -84, losses_db=2).mean_power_dbm(50e3) == pytest.approx(base + 1)


def test_rayleigh_probability_closed_form():
    for margin_db in (-10.0, 0.0, 10.0):
        assert ac.reception_probability(margin_db) == pytest.approx(math.exp(-(10 ** (-margin_db / 10))))


def test_rice_probability_against_rayleigh():
    pytest.importorskip("scipy")
    margins = np.array([-5.0, 0.0, 5.0])
    assert ac.reception_probability(margins, k_db=-60) == pytest.approx(ac.reception_probability(margins), rel=1e-3)
    assert ac.reception_probability(3.0, k_db=10) > ac.reception_probability(3.0)  # a direct wave helps above the threshold
    assert ac.reception_probability(-10.0, k_db=10) < ac.reception_probability(-10.0)  # and hurts below it
    assert np.all(np.diff(ac.reception_probability(np.linspace(-20, 20, 50), k_db=6)) > 0)


def test_slant_range():
    assert ac.slant_range(10.0, 1010.0, 0.0) == pytest.approx(1000.0)
    assert ac.slant_range(0.0, 0.0, 30e3) == pytest.approx(30e3, rel=1e-4)  # chord slightly shorter than the arc
    assert ac.slant_range(0.0, 1000.0, 30e3) == pytest.approx(math.hypot(30e3, 1000.0), rel=1e-3)


def test_reception_zone_matches_the_link_budget_on_flat_terrain():
    """Margin 0 dB at a slant range of 100 km; Rayleigh needs +9.77 dB for 90 % and +1.59 dB for 50 %."""
    link = ac.LinkBudget(freq_mhz=1090, eirp_dbm=51, rx_gain_dbi=0, sensitivity_dbm=51 - 133.2)
    az, d, probability, _ = ac.compute_reception(ac.Terrain([Flat()]), LAT, LON, 10.0, [1000.0], 450e3, 100.0, 90.0,
                                                 4 / 3, link)
    _, _, visible, _ = ac.compute_coverage(ac.Terrain([Flat()]), LAT, LON, 10.0, [1000.0], 450e3, 100.0, 90.0, 4 / 3)
    r50 = ac.continuous_range(probability >= 0.5, d)[0, 0]
    r90 = ac.continuous_range(probability >= 0.9, d)[0, 0]
    los = ac.continuous_range(visible, d)[0, 0]
    assert r90 == pytest.approx(math.sqrt((100e3 * 10 ** (-9.77 / 20)) ** 2 - 990.0**2), abs=300)
    assert r50 == pytest.approx(math.sqrt((100e3 * 10 ** (-1.59 / 20)) ** 2 - 990.0**2), abs=300)
    assert r90 < r50 < los
    assert np.all(probability[~visible] == 0)


def test_reception_is_zero_behind_the_terrain():
    link = ac.LinkBudget(1090, 51, 0, -100)  # generous budget: only the terrain can stop the signal
    ridge = EastProfile(lambda dx: 2000 * np.exp(-((dx - 20000) ** 2) / (2 * 1500**2)))
    az, d, probability, _ = ac.compute_reception(ac.Terrain([ridge]), LAT, LON, 10.0, [3000.0], 100e3, 100.0, 90.0,
                                                 4 / 3, link)
    east, west = list(az).index(90), list(az).index(270)
    assert probability[0, east, 400] == 0 and probability[0, west, 400] > 0.9  # 40 km: shadow vs open


def test_cli_fading_requires_the_link_budget(capsys):
    with pytest.raises(SystemExit) as error:
        ac.main([str(LAT), str(LON), "--flat", "--fading", "rice", "--freq-mhz", "1090"])
    assert error.value.code == 2
    message = capsys.readouterr().err
    assert "--eirp-dbm" in message and "--sensitivity-dbm" in message and "--k-db" in message


@pytest.mark.parametrize("fading", ["rayleigh", "rice"])
def test_cli_fading_writes_a_kml_of_the_received_zone(tmp_path, capsys, fading):
    pytest.importorskip("scipy")
    out = tmp_path / "reception.kml"
    ac.main([str(LAT), str(LON), "--flat", "--agl", "--az-step", "10", "-o", str(out), "--fading", fading,
             "--freq-mhz", "1090", "--eirp-dbm", "51", "--rx-gain-dbi", "3", "--sensitivity-dbm", "-84", "--k-db", "10"])
    text = out.read_text(encoding="utf-8")
    assert "received max" in text and "Reception probability at least 90%" in text
    assert len(ET.fromstring(text.encode()).findall(f".//{KML}Polygon")) > 0
    assert "max received" in capsys.readouterr().out


@pytest.mark.parametrize("clamp", [False, True])
def test_agl_kml_is_drawn_at_the_real_height_unless_clamped(tmp_path, capsys, clamp):
    out = tmp_path / "agl.kml"
    ac.main([str(LAT), str(LON), "--flat", "--agl", "--alt-ft", "100", "--az-step", "10", "-o", str(out),
             *(["--clamp-to-ground"] if clamp else [])])
    text = out.read_text(encoding="utf-8")
    if clamp:
        assert "<altitudeMode>clampToGround</altitudeMode>" in text and "relativeToGround" not in text
    else:
        assert "<altitudeMode>relativeToGround</altitudeMode>" in text and "clampToGround" not in text
        assert f",{100 * FT:.1f} " in text  # 30.5 m: the altitude of the vertices
