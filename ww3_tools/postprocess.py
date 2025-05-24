import numpy as np
import xarray as xr
import dask
from datetime import timedelta, datetime
from glob import glob

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

def find_nearest_point(ds, lon_ts, lat_ts):
    """
    Find the nearest indices of the model curvilinear grid to a specified lon, lat coordinate:

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

