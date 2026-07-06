import numpy as np
import xarray as xr
import dask
from datetime import timedelta, datetime
from glob import glob
from scipy.spatial import cKDTree

def get_ds_map(fname,
               time=slice(None) # Can be a single integer, or a slice. Default it to select all time
               ):
    '''
    flexible method to get the xarray dataset for either a
    single or multiple WW3 map files 
    The dimensions latitude,longitude are renamed to y,x as we have a curvilinear grid and so
    latitude,longitude are 2D variables, not dimensions
    '''
    if ('*' in fname) or ('?' in fname) or ('[' in fname):
        
        ds = xr.open_mfdataset(
            fname,
            combine='nested',
            concat_dim='time',
            preprocess=lambda ds: ds.rename_dims({'latitude': 'y', 'longitude': 'x'}),
            coords='minimal',
            compat='override'
        )
    else:
        ds = xr.open_dataset(fname)
        ds = ds.rename_dims({'latitude': 'y', 'longitude': 'x'})
    # subset on time
    if isinstance(time, int): # select only a single timestep if specified as an integer
        ds=ds.isel(time=time)
    else: # otherwise use the input time slice 
        ds=ds.sel(time=time)
    return ds

def get_boundary(fname):
    '''
    Return lon,lat of perimeter around the curvilinear grid (i.e. coordinates of bounding cells)
    '''
    ds=get_ds_map(fname)
    lon=ds.longitude.values
    lat = ds.latitude.values
    lon_out = np.hstack((lon[0:, 0], lon[-1, 1:-1],
                     lon[-1::-1, -1], lon[0, -2::-1]))
    lat_out = np.hstack((lat[0:, 0], lat[-1, 1:-1],
                     lat[-1::-1, -1], lat[0, -2::-1]))
    return lon_out, lat_out

def find_nearest_point(ds, lon_ts, lat_ts):
    """
    Find the nearest indices of the model curvilinear grid to a specified lon, lat coordinate.

    If the dataset carries a MAPSTA mask, only active sea cells
    (MAPSTA == 1) are considered, so a request very close to a coast
    returns the nearest wet point rather than a dry cell with junk
    values.

    Returns:
    - j :the nearest y index
    - i :the nearest x index

    j,i can be used in xarrays built-in xr.isel() function to extract data at this grid point

    """

    lon = ds.longitude.values
    lat = ds.latitude.values

    # Calculate the distance between (lon_ts, lat_ts) and all grid points
    distance = ((lon - lon_ts) ** 2 +
                (lat - lat_ts) ** 2) ** 0.5

    # Restrict the search to active sea cells if a MAPSTA mask is present.
    if 'MAPSTA' in ds:
        distance = np.where(ds['MAPSTA'].values == 1, distance, np.inf)

    # Find the indices of the minimum distance
    # unravel_index method Converts a flat index or array of flat indices into a tuple of coordinate
    # arrays: https://numpy.org/doc/stable/reference/generated/numpy.unravel_index.html
    min_index = np.unravel_index(distance.argmin(), distance.shape)

    j, i = min_index

    return j, i

def get_ts(fname, lon_ts, lat_ts,
                time=slice(None),
                nc_out=None,
                ):
    """
           Extract a ts from the map output
    """
    
    ds = get_ds_map(fname, time=time)
    
    j, i = find_nearest_point(ds, lon_ts, lat_ts) 
    
    ds = ds.isel(y=j,x=i)
    
    if nc_out is not None:
        print('writing the netcdf file')
        ds.to_netcdf(nc_out)

    return ds

def build_kdtree(ds):
    '''
    Build a scipy cKDTree from the 2D lon/lat arrays of a WW3 dataset.

    Parameters:
    - ds: xarray dataset from get_ds_map()

    Returns:
    - tree: scipy.spatial.cKDTree built from the flattened grid coordinates
    - grid_shape: tuple (ny, nx) for unraveling flat indices back to 2D
    '''
    lon = ds.longitude.values
    lat = ds.latitude.values
    grid_shape = lon.shape
    coords = np.column_stack([lon.ravel(), lat.ravel()])
    tree = cKDTree(coords)
    return tree, grid_shape

def extract_at_obs_points(fname, obs_lon, obs_lat, obs_time, var='hs',
                          max_dist_deg=0.1, tree=None, grid_shape=None):
    '''
    Extract model values at multiple observation points (lon, lat, time).
    Uses KDTree for efficient spatial nearest-neighbor lookup on the curvilinear grid,
    and nearest-in-time matching to model output.

    Parameters:
    - fname: path to a single monthly WW3 netCDF file
    - obs_lon: 1D array of observation longitudes
    - obs_lat: 1D array of observation latitudes
    - obs_time: 1D array of observation times (numpy datetime64)
    - var: variable name to extract (default 'hs')
    - max_dist_deg: maximum distance in degrees for spatial matching (default 0.1 ~ 10km)
    - tree: pre-built cKDTree (optional, built from file if None)
    - grid_shape: tuple (ny, nx) required if tree is provided

    Returns:
    - result: 1D numpy array of model values (NaN where tolerance exceeded)
    '''
    ds = get_ds_map(fname)

    if tree is None:
        tree, grid_shape = build_kdtree(ds)

    obs_lon = np.asarray(obs_lon)
    obs_lat = np.asarray(obs_lat)
    obs_time = np.asarray(obs_time, dtype='datetime64[ns]')

    # Spatial: query KDTree for nearest grid points
    dist, idx = tree.query(np.column_stack([obs_lon, obs_lat]))
    j_arr, i_arr = np.unravel_index(idx, grid_shape)
    valid_space = dist < max_dist_deg

    # Temporal: find nearest model timestep
    model_times = ds.time.values.astype('datetime64[ns]')
    time_idx = np.searchsorted(model_times, obs_time)
    # clamp to valid range
    time_idx = np.clip(time_idx, 0, len(model_times) - 1)
    # check if the previous timestep is closer
    prev_idx = np.clip(time_idx - 1, 0, len(model_times) - 1)
    diff_curr = np.abs(model_times[time_idx] - obs_time)
    diff_prev = np.abs(model_times[prev_idx] - obs_time)
    use_prev = diff_prev < diff_curr
    time_idx[use_prev] = prev_idx[use_prev]

    # Apply temporal tolerance (30 minutes)
    time_diff = np.abs(model_times[time_idx] - obs_time)
    valid_time = time_diff < np.timedelta64(30, 'm')

    valid = valid_space & valid_time

    # Load only needed timesteps
    unique_tidx = np.unique(time_idx[valid])
    if len(unique_tidx) == 0:
        ds.close()
        return np.full(len(obs_lon), np.nan)

    data_subset = ds[var].isel(time=unique_tidx).values  # (n_unique_times, ny, nx)

    # Map global time indices to local indices in the subset
    tidx_map = {int(t): k for k, t in enumerate(unique_tidx)}
    local_tidx = np.array([tidx_map.get(int(time_idx[n]), -1) for n in range(len(obs_lon))])

    result = np.full(len(obs_lon), np.nan)
    valid_idx = np.where(valid)[0]
    result[valid_idx] = data_subset[local_tidx[valid_idx], j_arr[valid_idx], i_arr[valid_idx]]

    ds.close()
    return result

