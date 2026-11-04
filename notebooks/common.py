"""Shared settings and functions for the figure notebooks.

Set two environment variables before starting Jupyter:
    AMOC_DATA_DIR    folder holding the intermediate data files listed in data/MANIFEST.csv (required)
    AMOC_FIGURE_DIR  folder for the saved figures (optional; default: ../output next to this folder)
"""
import os
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.signal import butter, filtfilt

BRANCH_YEAR = 420        # model year of the abrupt CO2 increase, branched from the control run
RECOVERY_YEAR = 2428     # year since 8xCO2 when the AMOC recovery begins (the vertical line in time series)
SMOOTH_CUTOFF = 15       # low-pass filter cutoff period [years]
REF_PRESSURE = 1000      # reference pressure of sigma_1 [dbar]
PCT_HI, PCT_LO = 95, 5   # upper and lower area-weighted percentiles of surface density

HLNA_LAT = (50, 65)      # high-latitude North Atlantic band [degrees N]
STNA_LAT = (30, 40)      # subtropical North Atlantic band [degrees N]
SO_LAT = (-75, -55)      # Southern Ocean channel band [degrees N]

X_START_PAD = 50         # years shown before year 0, so that year-0 values stand clear of the axis


def data_path(filename):
    """Full path of an intermediate data file inside AMOC_DATA_DIR."""
    data_dir = os.environ.get('AMOC_DATA_DIR')
    if not data_dir:
        raise RuntimeError('Set AMOC_DATA_DIR to the folder holding the files listed in data/MANIFEST.csv')
    path = Path(data_dir) / filename
    if not path.exists():
        raise FileNotFoundError(f'{filename} is not in {data_dir}; see data/MANIFEST.csv')
    return str(path)


def figure_path(filename):
    """Full path for a saved figure, creating the output folder if needed."""
    default = Path(__file__).resolve().parent.parent / 'output'
    folder = Path(os.environ.get('AMOC_FIGURE_DIR', default))
    folder.mkdir(parents=True, exist_ok=True)
    return str(folder / filename)


def butterworth_lowpass(series, cutoff_period=SMOOTH_CUTOFF, order=4, sampling_rate=1.0):
    """Zero-phase low-pass filter of a yearly series; cutoff_period in years. Same units as the input."""
    nyquist = 0.5 * sampling_rate
    b, a = butter(order, (1.0 / cutoff_period) / nyquist, btype='low', analog=False)
    return filtfilt(b, a, series)


def weighted_percentile(values, weights, percentile):
    """Area-weighted percentile (0-100) of values (year, cell) for each year; same units as values.

    Each cell's value is placed at the middle of its share of the total area, and the percentile is interpolated
    between neighboring cells. Larger cells therefore count for more.
    """
    order = np.argsort(values, axis=1)
    sorted_values = np.take_along_axis(values, order, axis=1)
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights, axis=1)
    cumulative = (cumulative - 0.5 * sorted_weights) / cumulative[:, -1:]
    return np.array([np.interp(percentile / 100., cumulative[i], sorted_values[i])
                     for i in range(sorted_values.shape[0])])


def eos_sigma1(salt, temp):
    """sigma_1 [kg/m^3, minus 1000] from salinity [g/kg] and temperature [degC], with POP's equation of state."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)  # pop_tools warns about pkg_resources when imported
        import pop_tools
    return pop_tools.eos(salt=np.asarray(salt, dtype='float64'), temp=np.asarray(temp, dtype='float64'),
                         pressure=REF_PRESSURE) - 1000.


def amoc_strength(filename):
    """Annual AMOC strength [Sv]: maximum of the Eulerian-mean streamfunction north of 45N, at any depth.

    Plot against the returned .time coordinate. The 2x and control files count from 0, the 4x and 8x files from 1.
    """
    moc = (xr.open_dataset(data_path(filename), decode_times=False)
           .isel(moc_comp=0, drop=True).mean(dim='transport_reg').moc)
    return moc.sel(lat_aux_grid=slice(45, None), drop=True).max(dim=('lat_aux_grid', 'moc_z'))


def sigma_cache(run):
    """Per-cell late-winter surface sigma_1 of a run: FMA in the N. Atlantic, ASO in the Southern Ocean; (year, cell)."""
    return xr.open_dataset(data_path(f'FMAASOmean_sigma1_cells_{run}_natl20-80N_so75-50S.nc')).load()


def _smoothed_percentiles(cache, selected):
    """Smoothed upper and lower percentiles of sigma_1 over the selected cells, skipping cells with any NaN."""
    values, area = cache.season_mean.values[:, selected], cache.area.values[selected]
    has_data = ~np.isnan(values).any(axis=0)
    return {q: butterworth_lowpass(weighted_percentile(values[:, has_data], area[has_data], q))
            for q in (PCT_HI, PCT_LO)}


def density_ranges(run, natl_lat):
    """Yearly smoothed p95 and p5 of late-winter surface sigma_1 in a N. Atlantic band and the S.O. channel region.

    Returns a dict with 'years' and the series 'natl_hi', 'natl_lo', 'so_hi', 'so_lo' [kg/m^3, minus 1000].
    """
    cache = sigma_cache(run)
    in_natl = ((cache.basin == 'natl') & (cache.lat >= natl_lat[0]) & (cache.lat <= natl_lat[1])).values
    in_so = ((cache.basin == 'so') & (cache.lat >= SO_LAT[0]) & (cache.lat <= SO_LAT[1])).values
    northern, southern = _smoothed_percentiles(cache, in_natl), _smoothed_percentiles(cache, in_so)
    return {'years': cache.year.values,
            'natl_hi': northern[PCT_HI],
            'natl_lo': northern[PCT_LO],
            'so_hi': southern[PCT_HI],
            'so_lo': southern[PCT_LO]}


def split_density_change(temp, salt, temp_ref, salt_ref):
    """Change in sigma_1 from a control reference, split into temperature and salinity parts [kg/m^3].

    temp [degC] and salt [g/kg] are yearly series. Each part changes one tracer and holds the other at its control
    value; 'total' differs slightly from their sum because the equation of state is nonlinear.
    """
    sigma_ref = eos_sigma1(salt_ref, temp_ref)
    return {'total': eos_sigma1(salt, temp) - sigma_ref,
            'temperature': eos_sigma1(np.full_like(temp, salt_ref), temp) - sigma_ref,
            'salinity': eos_sigma1(salt, np.full_like(salt, temp_ref)) - sigma_ref}


def _densest_mean_ts(sigma, temp, salt, area):
    """Area-weighted T [degC] and S [g/kg] over each year's densest (100 - PCT_HI)% of cells by area."""
    threshold = weighted_percentile(sigma, area, PCT_HI)
    top_area = (sigma >= threshold[:, None]) * area   # zero outside each year's densest set
    return (top_area * temp).sum(1) / top_area.sum(1), (top_area * salt).sum(1) / top_area.sum(1)


def densest_water_breakdown(run):
    """Density change of the densest 5% of HLNA water relative to the control; returns years and the split.

    Each year's densest cells are chosen from FMA-mean sigma_1, and their FMA-mean T and S are averaged by area. The
    control reference applies the same selection to every control year and then averages over those years.
    """
    ts = xr.open_dataset(data_path(f'FMAmean_TS_cells_{run}_natl20-80N.nc')).load()
    ctrl = xr.open_dataset(data_path('surface_TS_sigma1_cells_ctrl_natl20-80N.nc')).load()
    # The sigma_1 cache stores the same N. Atlantic cells first, in the same order, then the Southern Ocean cells
    sigma = sigma_cache(run).isel(cell=slice(0, ts.sizes['cell']))
    assert (sigma.basin == 'natl').all() and np.allclose(sigma.lat, ts.lat) and np.allclose(sigma.lon, ts.lon)
    assert np.allclose(ctrl.lat, ts.lat) and np.allclose(ctrl.lon, ts.lon)
    assert np.array_equal(ts.year.values, sigma.year.values)

    in_hlna = ((ts.lat >= HLNA_LAT[0]) & (ts.lat <= HLNA_LAT[1])).values
    run_fields = [sigma.season_mean.values[:, in_hlna], ts.TEMP.values[:, in_hlna], ts.SALT.values[:, in_hlna]]
    ctrl_fields = [ctrl.SIGMA1_FMA.values[:, in_hlna], ctrl.TEMP_FMA.values[:, in_hlna], ctrl.SALT_FMA.values[:, in_hlna]]

    # Keep cells with data in every year of both runs, so that both use the same set of cells
    has_data = ~np.any([np.isnan(field).any(0) for field in run_fields + ctrl_fields], axis=0)
    area = ts.area.values[in_hlna][has_data]
    temp_top, salt_top = _densest_mean_ts(*[field[:, has_data] for field in run_fields], area)
    temp_top_c, salt_top_c = _densest_mean_ts(*[field[:, has_data] for field in ctrl_fields], area)
    return ts.year.values, split_density_change(temp_top, salt_top, temp_top_c.mean(), salt_top_c.mean())


def band_mean_breakdown(lat_band=(60, 65)):
    """As densest_water_breakdown, for annual-mean T and S averaged over all 8xCO2 cells in a latitude band."""
    annual = xr.open_dataset(data_path('ANNmean_TS_cells_8x_natl20-80N.nc')).load()
    ctrl = xr.open_dataset(data_path('surface_TS_sigma1_cells_ctrl_natl20-80N.nc')).load()
    assert np.allclose(ctrl.lat, annual.lat) and np.allclose(ctrl.lon, annual.lon)
    in_band = ((annual.lat >= lat_band[0]) & (annual.lat <= lat_band[1])).values
    temp, salt = annual.TEMP.values[:, in_band], annual.SALT.values[:, in_band]
    temp_c, salt_c = ctrl.TEMP_ANN.values[:, in_band], ctrl.SALT_ANN.values[:, in_band]

    has_data = ~(np.isnan(temp).any(0) | np.isnan(salt).any(0) | np.isnan(temp_c).any(0) | np.isnan(salt_c).any(0))
    area = annual.area.values[in_band][has_data]

    def area_mean(field):
        """Area-weighted mean over the band for each year."""
        return (field[:, has_data] * area).sum(1) / area.sum()

    return annual.year.values, split_density_change(area_mean(temp), area_mean(salt),
                                                    area_mean(temp_c).mean(), area_mean(salt_c).mean())


def freshwater_terms(run):
    """Annual freshwater flux terms over the Atlantic + GIN 60-65N box [mm/day], positive into the ocean.

    Returns years (model years) and a dict. 'Sea ice' is melt plus the salt flux expressed as freshwater, and 'Runoff'
    is liquid plus ice runoff. A flux in kg/m^2/s times 86400 gives mm/day, since 1 kg/m^2 of water is a 1 mm layer.
    """
    terms = xr.open_dataset(data_path(f'SFWF_terms_ann_atlGIN60-65N_{run}.nc')).sel(region='box')
    in_mm_per_day = {name: terms[name].values * 86400 for name in terms.data_vars}
    return terms.year.values, {'Precipitation': in_mm_per_day['PREC_F'],
                               'Evaporation': in_mm_per_day['EVAP_F'],
                               'Runoff': in_mm_per_day['ROFF_F'] + in_mm_per_day['IOFF_F'],
                               'Sea ice': in_mm_per_day['MELT_F'] + in_mm_per_day['SALT_F_fw'],
                               'Total': in_mm_per_day['SFWF']}
