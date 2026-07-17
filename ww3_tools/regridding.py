"""
Regridding utilities for WW3 output.

Mirrors the API and conventions of crocotools_py.regridding.regrid_tier3
so that the front-end consumes WW3 fields on the same regular lat/lon
grid as CROCO fields, with the same chunking and metadata shape.
"""

import os
import re
import shutil
import subprocess
from datetime import datetime
from glob import glob

import dask
import dask.array as da
import matplotlib.path as mplPath
import numpy as np
import xarray as xr
from dask import delayed
from scipy.interpolate import griddata
from scipy.spatial import Delaunay


DEFAULT_VARS = ['hs', 't0m1', 'dir', 'dp']

# Wave directions, in degrees. These are circular and must not be interpolated
# as plain scalars: a linear blend of 350 and 10 gives 180 — the opposite
# direction. They are interpolated via their unit-vector components instead.
DIRECTIONAL_VARS = ('dir', 'dp')


def _zarr_compressor_encoding(chunks):
    """Return the (key, value) pair that goes into a per-variable encoding
    dict to request Blosc/Zstd compression. zarr v2 uses
    'compressor': numcodecs.Blosc; zarr v3 uses 'compressors': BloscCodec."""
    import zarr
    major = int(zarr.__version__.split('.')[0])
    if major >= 3:
        from zarr.codecs import BloscCodec
        return ('compressors', BloscCodec(cname='zstd', clevel=3, shuffle='shuffle'))
    else:
        from numcodecs import Blosc
        return ('compressor', Blosc(cname='zstd', clevel=3, shuffle=Blosc.SHUFFLE))


def regrid_tier3(fname_in, dir_out, spacing=0.01, varList=DEFAULT_VARS,
                 output_format='netcdf', doi_link=None, derive_tp=True):
    """
    Regrid WW3 ww3_ounf.nc output (curvilinear lon/lat) onto a regular
    lat/lon grid suitable for web visualisation.

    Parameters
    ----------
    fname_in : str or list of str
        Path or glob to WW3 ww3_ounf.nc file(s).
    dir_out : str
        Output directory.
    spacing : float
        Output grid spacing in degrees.
    varList : list of str
        Native WW3 variables to interpolate. All assumed (time, lat, lon).
    output_format : {'netcdf', 'zarr'}
        Zarr uses Blosc/Zstd compression with time=1 chunking, suited
        for per-map fetch from a web tier.
    doi_link : str, optional
        DOI added to dataset attrs.
    derive_tp : bool
        If True and `fp` is present, also publish `tp = 1/fp` (peak
        period, seconds).

    Notes
    -----
    Output basename swaps `ww3_ounf` for `ww3_t3` and ends in `.zarr`
    when output_format='zarr'.
    """

    if output_format not in ('netcdf', 'zarr'):
        raise ValueError(f"output_format must be 'netcdf' or 'zarr', got '{output_format}'")

    if isinstance(fname_in, str):
        fname_in = sorted(glob(fname_in)) if '*' in fname_in else [fname_in]

    os.makedirs(dir_out, exist_ok=True)

    for file in fname_in:
        print(f'Opening: {file}')
        ds = xr.open_dataset(file)

        # Optionally derive Tp from fp; do it before deciding which vars survive.
        if derive_tp and 'fp' in ds and 'tp' not in ds:
            tp = 1.0 / ds['fp'].where(ds['fp'] > 0)
            tp.attrs = {'units': 's', 'long_name': 'peak wave period'}
            ds = ds.assign(tp=tp)

        # Filter varList to those actually present, plus tp if it was just added.
        active_vars = [v for v in varList if v in ds]
        if derive_tp and 'tp' in ds and 'tp' not in active_vars:
            active_vars.append('tp')
        if not active_vars:
            raise RuntimeError(f'None of {varList} found in {file}')

        Nt = ds.time.size
        lon_2d = ds['longitude'].values
        lat_2d = ds['latitude'].values

        # Drop land cells so nearest-neighbour interp can't pick land values
        # near the coast. MAPSTA==1 means active sea cell.
        if 'MAPSTA' in ds:
            wet = ds['MAPSTA'].values == 1
        else:
            wet = np.ones_like(lat_2d, dtype=bool)

        # Build the boundary polygon by walking the four edges of the
        # curvilinear grid in order (same idiom as CROCO regrid_tier3).
        lon_b = np.hstack(
            (lon_2d[:, 0], lon_2d[-1, 1:], lon_2d[-1::-1, -1], lon_2d[0, -2::-1])
        )
        lat_b = np.hstack(
            (lat_2d[:, 0], lat_2d[-1, 1:], lat_2d[-1::-1, -1], lat_2d[0, -2::-1])
        )

        lon_min = np.floor(np.min(lon_b) / spacing) * spacing
        lon_max = np.ceil(np.max(lon_b) / spacing) * spacing
        lat_min = np.floor(np.min(lat_b) / spacing) * spacing
        lat_max = np.ceil(np.max(lat_b) / spacing) * spacing

        Nlon = int(np.rint((lon_max - lon_min) / spacing)) + 1
        Nlat = int(np.rint((lat_max - lat_min) / spacing)) + 1
        lon_out = np.linspace(lon_min, lon_max, Nlon, endpoint=True)
        lat_out = np.linspace(lat_min, lat_max, Nlat, endpoint=True)
        lon_grd, lat_grd = np.meshgrid(lon_out, lat_out)

        # Polygon mask (NaN outside the curvilinear boundary)
        poly = mplPath.Path(np.column_stack([lon_b, lat_b]))
        mask_out = np.zeros((Nlat, Nlon))
        for y in range(Nlat):
            for x in range(Nlon):
                if poly.contains_point((lon_grd[y, x], lat_grd[y, x])):
                    mask_out[y, x] = 1.0

        # Exclude land from the output grid. The boundary polygon above only
        # clips to the domain's outer edge, so land inside it would otherwise
        # be filled by interpolation across the coastline. Nearest-interpolate
        # the native land/sea mask (all cells, wet=1 / dry=0) onto the output
        # grid and fold dry cells into mask_out so they become NaN in the
        # `* mask_out / mask_out` step.
        wet_out = griddata(
            np.column_stack([lon_2d.ravel(), lat_2d.ravel()]),
            wet.ravel().astype(float),
            (lon_grd, lat_grd), method='nearest',
        )
        mask_out = mask_out * (wet_out > 0.5)

        # Source points restricted to wet cells
        wet_flat = wet.ravel()
        lonlat_input = np.column_stack(
            [lon_2d.ravel()[wet_flat], lat_2d.ravel()[wet_flat]]
        )

        # Precompute the Delaunay triangulation and barycentric weights ONCE.
        # The source mesh and the output points are identical for every
        # (time, variable) tile, so the expensive triangulation and per-point
        # simplex lookup are hoisted out of the per-chunk loop — each field
        # then collapses to a gather + weighted sum. Same idiom as
        # crocotools_py.regridding.regrid_tier3.
        #
        # Linear (rather than nearest) matters here: the source grid is much
        # coarser than the output, so nearest-neighbour paints Voronoi wedges
        # and leaves radial streaks fanning out to the domain edge. Linear
        # returns NaN outside the convex hull of the wet points instead.
        xi_flat = np.column_stack([lon_grd.ravel(), lat_grd.ravel()])
        tri = Delaunay(lonlat_input)
        simplex = tri.find_simplex(xi_flat)          # -1 outside the hull
        valid = simplex >= 0
        safe_simplex = np.where(valid, simplex, 0)
        vtx = tri.simplices[safe_simplex]            # (Nout, 3)
        T = tri.transform[safe_simplex]              # (Nout, 3, 2)
        delta = xi_flat - T[:, 2, :]
        bary = np.einsum('nij,nj->ni', T[:, :2, :], delta)
        weights = np.column_stack([bary, 1.0 - bary.sum(axis=1)])
        weights[~valid] = 0.0
        valid_2d = valid.reshape(Nlat, Nlon)

        def _interp(values_flat):
            gathered = values_flat[vtx]              # (Nout, 3)
            out = np.einsum('ni,ni->n', gathered, weights).reshape(Nlat, Nlon)
            return np.where(valid_2d, out, np.nan)

        @delayed
        def compute_2d_chunk(t, variable, circular=False):
            vals = np.asarray(variable[t]).ravel()[wet_flat]
            if circular:
                rad = np.deg2rad(vals)
                out = np.rad2deg(
                    np.arctan2(_interp(np.sin(rad)), _interp(np.cos(rad)))
                ) % 360.0
            else:
                out = _interp(vals)
            return out * mask_out / mask_out

        out_time = {v: [] for v in active_vars}
        for t in range(Nt):
            for v in active_vars:
                out_time[v].append(
                    da.from_delayed(
                        compute_2d_chunk(t, ds[v].values,
                                         circular=v in DIRECTIONAL_VARS),
                        shape=(Nlat, Nlon),
                        dtype=float,
                    )
                )
        for v in active_vars:
            out_time[v] = da.stack(out_time[v], axis=0)

        # Build output dataset
        data_vars = {
            v: xr.Variable(['time', 'latitude', 'longitude'], out_time[v], ds[v].attrs)
            for v in active_vars
        }
        coords = {
            'longitude': xr.Variable(['longitude'], lon_out,
                                     {'units': 'degrees_east',
                                      'standard_name': 'longitude'}),
            'latitude':  xr.Variable(['latitude'], lat_out,
                                     {'units': 'degrees_north',
                                      'standard_name': 'latitude'}),
            'time':      xr.Variable(['time'], ds.time.values, ds.time.attrs),
        }
        data_out = xr.Dataset(data_vars=data_vars, coords=coords)
        data_out.attrs['title'] = 'Regridded WW3 output created by ww3_tools.regridding.regrid_tier3'
        data_out.attrs['source'] = file
        data_out.attrs['history'] = 'Created on ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        data_out.attrs['conventions'] = 'CF-1.8'
        if doi_link is not None:
            data_out.attrs['doi'] = f'https://doi.org/{doi_link}'

        chunksizes = {'time': 1 if output_format == 'zarr' else ds.time.size}
        default_chunksizes = {d: data_out.sizes[d] for d in data_out.sizes}

        if output_format == 'zarr':
            comp_key, comp_val = _zarr_compressor_encoding(None)
            encoding = {
                var: {
                    'dtype': 'float32',
                    'chunks': [chunksizes.get(d, default_chunksizes[d])
                               for d in data_out[var].dims],
                    comp_key: comp_val,
                }
                for var in data_out.data_vars
            }
        else:
            encoding = {
                var: {
                    'dtype': 'float32',
                    'chunksizes': [chunksizes.get(d, default_chunksizes[d])
                                   for d in data_out[var].dims],
                }
                for var in data_out.data_vars
            }
        encoding['latitude'] = {'dtype': 'float32'}
        encoding['longitude'] = {'dtype': 'float32'}

        # Output filename: replace 'ww3_ounf' with 'ww3_t3'; switch ext if zarr
        basename = os.path.basename(file)
        match = re.search(r'(\.nc\.\d+)$|(\.nc)$', basename)
        if match is None:
            raise RuntimeError(f'Cannot parse extension of {basename}')
        no_ext = basename[:match.start()].replace('ww3_ounf', 'ww3_t3')
        ext = match.group(0)
        if output_format == 'zarr':
            ext = ext.replace('.nc', '.zarr')
        fname_out = os.path.abspath(os.path.join(dir_out, no_ext + ext))

        if os.path.exists(fname_out):
            if output_format == 'zarr':
                shutil.rmtree(fname_out)
            else:
                os.remove(fname_out)

        if output_format == 'zarr':
            print('Generating Zarr store')
            write_op = data_out.to_zarr(fname_out, encoding=encoding,
                                        mode='w', compute=False)
            print('Writing Zarr store')
            dask.compute(write_op)
        else:
            print('Generating NetCDF data')
            write_op = data_out.to_netcdf(fname_out, encoding=encoding,
                                          mode='w', compute=False)
            print('Writing NetCDF file')
            dask.compute(write_op)

        subprocess.call(['chmod', '-R', '775', fname_out])
        print(f'Created: {fname_out}')
