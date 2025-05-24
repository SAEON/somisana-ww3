import os
import numpy as np
import xarray as xr
import wavespectra
from wavespectra.construct import construct_partition

def croco_grd_2_ww3(croco_grdfile,
                    lat_out='lat.dat',
                    lon_out='lon.dat',
                    depth_out='depth.dat',
                    mask_out='mask.dat'):

    '''
    concert a croco grid file to separate files which can be ingested by WW3
    '''
        
    ds = xr.open_dataset(croco_grdfile)

    lon=ds.lon_rho.values
    lat=ds.lat_rho.values
    mask=ds.mask_rho.values.astype(np.int64)
    h=ds.h.values*-1

    # convert the sea points on the model boundaries to active wave boundary conditions - i.e. make equal 2
    # Top edge
    mask[0, :] = np.where(mask[0, :] == 1, 2, mask[0, :])
    # Bottom edge
    mask[-1, :] = np.where(mask[-1, :] == 1, 2, mask[-1, :])
    # Left edge (excluding the top and bottom corners already processed)
    mask[1:-1, 0] = np.where(mask[1:-1, 0] == 1, 2, mask[1:-1, 0])
    # Right edge (excluding the top and bottom corners already processed)
    mask[1:-1, -1] = np.where(mask[1:-1, -1] == 1, 2, mask[1:-1, -1])

    # write output files
    np.savetxt(lon_out, np.flip(lon, axis=0), fmt='%.6f')
    np.savetxt(lat_out, np.flip(lat, axis=0), fmt='%.6f')
    np.savetxt(depth_out, np.flip(h, axis=0), fmt='%.3f')
    np.savetxt(mask_out, np.flip(mask, axis=0), fmt="%d")


def gamma_from_tp(Tp, gamma_min=1.5, gamma_max=6.0, Tp_mid=9.0, k=0.8):
    """
    Estimate JONSWAP gamma (peakedness) as a smooth function of peak wave period (Tp).
    
    This was Chatgpt's suggested relationship, and so I'm sure we need to adjust it
    But the values seem reasonable and either way we'll be guessing the gamma value.
    I'm opreferring this dependence on Tp rather than using a constant value since gamma
    should physically depend on Tp
    
    Parameters:
    - Tp: Peak wave period (in seconds)
    - gamma_min: Minimum gamma value (broad spectrum, typical wind sea)
    - gamma_max: Maximum gamma value (narrow spectrum, mature swell)
    - Tp_mid: Tp at which gamma is halfway between min and max
    - k: Controls steepness of transition (higher = steeper)
    
    Returns:
    - gamma: Estimated JONSWAP peakedness parameter
    """
    return gamma_min + (gamma_max - gamma_min) / (1 + np.exp(-k * (Tp - Tp_mid)))


def get_bry_lon_lat(lon_file,
                        lat_file,
                        mask_file,
                        cmems_file):
    '''
    extract lon and lat coordinates from a .nc file with a regular grid file 
    (e.g. CMEMS or other I guess) which are within a tolerence (in degrees)
    of the open boundaries of your WW3 grid. This is useful to identify which
    coordinates to extract for generating spectral boundary files for the WW3 model
    '''
    
    # Read the lon,lat and mask for the ww3 grid
    lon_ww = np.loadtxt(lon_file)
    lat_ww = np.loadtxt(lat_file)
    mask_ww = np.loadtxt(mask_file)

    # Find the indices where mask is equal to 2 i.e. boundary points where spectral boundary conditions will be applied
    masked_indices = np.where(mask_ww == 2)

    # Extract the corresponding longitude and latitude values
    lon_ww = lon_ww[masked_indices]
    lat_ww = lat_ww[masked_indices]
    
    # extract the grid from the cmems file
    if isinstance(cmems_file, xr.Dataset): # handles the case of using an already extracted dataset as input
        ds = cmems_file.copy()
    else:
        ds = xr.open_mfdataset(cmems_file)
    lon = ds.longitude.values
    lat = ds.latitude.values
    ds.close()
    
    # Create a meshgrid of all (lon, lat) combinations from the file
    lon_grid, lat_grid = np.meshgrid(lon, lat)
    
    # Flatten the grid for easier comparison
    lon_flat = lon_grid.ravel()
    lat_flat = lat_grid.ravel()
    
    # Store matching coordinates
    matched_coords = []
    
    for mlon, mlat in zip(lon_ww, lat_ww):
        # Compute distance (in degrees) from current masked point
        dlon = np.abs(lon_flat - mlon)
        dlat = np.abs(lat_flat - mlat)
        close_idx = np.where((dlon <= 0.2) & (dlat <= 0.2))[0]
    
        for idx in close_idx:
            matched_coords.append((lon_flat[idx], lat_flat[idx]))
    
    # Convert to numpy array or process further
    matched_coords = np.array(matched_coords)
    
    # remove duplicates
    matched_coords = np.unique(matched_coords, axis=0)
    
    return matched_coords
    
def make_bry_spec_cmems(lon_file,
                        lat_file,
                        mask_file,
                        cmems_file,
                        output_dir,
                        num_freqs=32,
                        start_freq=0.0373,
                        freq_factor=1.1,
                        num_dirs=24,
                        fp_2_dspr = 120, # factor to get directional spreading from peak frequency
                        ):
    '''
    extract CMEMS data which are within a tolerence (in degrees) of the open boundaries 
    of your WW3 grid, and write a spectral boundary nc file for each point to be used
    as boundary conditions for a WW3 simulation
    '''
    matched_coords = get_bry_lon_lat(lon_file,
                            lat_file,
                            mask_file,
                            cmems_file)
    
    with open(os.path.join(output_dir,"spec.list"), "w") as f:
        for lon, lat in matched_coords:
            spec_file = f"cmems.lon{lon:.2f}.lat{lat:.2f}.spec.nc"
            cmems_2_spec(cmems_file,
                                  lon,
                                  lat,
                                  spec_file,
                                  num_freqs=num_freqs,
                                  start_freq=start_freq,
                                  freq_factor=freq_factor,
                                  num_dirs=num_dirs,
                                  fp_2_dspr = fp_2_dspr, # factor to get directional spreading from peak frequency
                                  )
            f.write(spec_file + "\n")
    
def cmems_2_spec(file_in,
                      lon_out,
                      lat_out,
                      file_out,
                      num_freqs=32,
                      start_freq=0.0373,
                      freq_factor=1.1,
                      num_dirs=24,
                      fp_2_dspr = 120, # factor to get directional spreading from peak frequency
                      ):
    '''
    this function is designed to work on the CMEMS Global Ocean Waves Reanalysis product:
    https://doi.org/10.48670/moi-00022
    As well as the CMEMS Global Ocean Waves Analysis and Forecast product:
    https://doi.org/10.48670/moi-00017
    
    It extracts a specified single point and recreates the full freq-dir spectrum from the wave partition data
    And writes a file which can be read by WW3 as a boundary condition file
    '''
    
    # setup the frequency axis of the wave spectrum
    freq = np.zeros(num_freqs)
    freq[0] = start_freq
    for i in range(1, num_freqs):
      freq[i] = freq[i-1] * freq_factor
    
    # set up the direction axis of the wave spectrum
    dir = np.linspace(0, 360, num_dirs, endpoint=False)

    # extract data for the specified lon,lat coordinate
    if isinstance(file_in, xr.Dataset): # handles the case of using an already extracted dataset as input
        ds = file_in.copy()
    else:
        ds = xr.open_mfdataset(file_in)
    ds = ds.interp(longitude=("site", [lon_out]), latitude=("site", [lat_out]),method='nearest')
    
    # read the partition data from the file
    # for each variable, concatenate the 2 swells and wind sea along a new 'partition' dimension
    #
    # significant wave height
    hs = xr.concat([ds.VHM0_SW1, ds.VHM0_SW2, ds.VHM0_WW], dim='partition')
    # mean wave direction
    dm = xr.concat([ds.VMDR_SW1, ds.VMDR_SW2, ds.VMDR_WW], dim='partition')
    #
    # To recreate the 2D spectrum we need Tp, not Tm01, for each partition
    # Unfortunately we only have Tm01 for the partitions, but we do have the bulk Tp (i.e. Tp for the whole spectrum)
    # So let's use the relationship between Tp and the primary swell Tm01 to get the primary swell Tp
    # (it is not correct to apply this to all partitions, only the one representative of Tp
    # which we are assuming will be the primary swell)
    Tm01_2_tp_SW1 = ds.VTPK/ds.VTM01_SW1 # conversion for the primary swell (also tested using a constant 1.2 but this dynamic approach gave better results)
    Tm01_2_tp_SW2 = 1.2 # conversion for the secondary swell
    Tm01_2_tp_WW = 1.1 # conversion for the wind wave
    #
    Tp = xr.concat([Tm01_2_tp_SW1 * ds.VTM01_SW1, 
                    Tm01_2_tp_SW2 * ds.VTM01_SW2, 
                    Tm01_2_tp_WW * ds.VTM01_WW], dim='partition')
    fp = 1 / Tp
    
    # spreading factor (standard deviation in the wave direction)
    # no spreading factor is provided with the partition data, so we need to use some estimate
    # Physically, lower period waves (higher frequencies) have higher spreading factors
    # Allowing for a rough linear estimate:
    dspr = fp_2_dspr * fp
    
    # ideally gamma, the peakedness parameter should also vary with Tp
    gamma = gamma_from_tp(Tp, gamma_min=1.5, gamma_max=6.0, Tp_mid=9.0, k=0.8)
    
    # Here's where the wavespectra magic happens to
    # recreate the full spectrum for each of the partitions
    efth = construct_partition(
        freq_name="jonswap",
        freq_kwargs={"freq": freq, "fp": fp, "gamma": gamma, "hs": hs},
        dir_name="cartwright",
        dir_kwargs={"dir": dir, "dm": dm, "dspr": dspr}
    )
    efth = efth.fillna(0)
    # and then sum over the partitions to create the full 2D spectrum
    efth = efth.sum(dim='partition')
    
    # plot just to check it looks OK
    # efth.isel(time=0,site=0).spec.plot()
    
    # Now format the data for writing to the output file
    # make a dataset
    ds_efth = xr.Dataset({'efth': efth})
    
    # Before writing we need to change the names of the coords to match the 
    # conventions in wavespectra.core.attributes
    ds_efth = ds_efth.rename({
        'longitude': 'lon',
        'latitude': 'lat'
    })
    # ds_efth['time'].encoding['unlimited_dims'] = ('time',)
    ds_efth['time'].attrs['units'] = 'hours since 1990-01-01'
    ds_efth['time'].encoding['dtype'] = 'float64'
    
    # make a specDataset object to get access to the built in functions to 
    # write spectral output files for specific models like WW3
    ds_spec = wavespectra.SpecDataset(ds_efth)
    
    # Finally write the WW3 spectral file
    ds_spec.to_ww3(file_out)
    
    ds.close()
    ds_efth.close()
    ds_spec.close()

def ncep_hind_2_spec():
    '''
    test function designed to work on NCEP's WW3 phase 2 hindcast wave partition output:
    https://polar.ncep.noaa.gov/waves/hindcasts/nopp-phase2.php
    It was never completed, but should be based largely on cmems_rean_2_spec
    I'm just sticking some lines here in case its useful in the future
    '''
    
    # spec_file='/mnt/c/Users/GilesF/Downloads/multi_reanal.partition.glo_30m.200911.nc'
    
    # # I'm specifically not decoding times below since it tries to decode period since it is units seconds
    # ds = xr.open_mfdataset(spec_file,decode_times=False)
    
    # # Manually decode the 'date' variable
    # ds['date'] = xr.decode_cf(ds[['date']].assign_coords(date=ds['date']))['date']
    
    # # read the partition data from the file
    # hs = ds.significant_wave_height
    # fp = 1/ds.peak_period
    # dm = ds.wave_direction
    # dspr = ds.direction_spreading
    
    # # Here's where the wavespectra magic happens - recreate the full spectrum from the partitions
    # # the first partition is actually the bulk paramter! so start at 1
    # efth = construct_partition(
    #     freq_name="jonswap",
    #     freq_kwargs={"freq": freq, "fp": fp[:,1:], "gamma": 2.0, "hs": hs[:,1:]},
    #     dir_name="cartwright",
    #     dir_kwargs={"dir": dir, "dm": dm[:,1:], "dspr": dspr[:,1:]}
    # )
    # efth = efth.sum(dim='partition')
    # # efth.isel(date=1,site=0).spec.plot()
    
    # # make a dataset
    # ds_efth = xr.Dataset({'efth': efth})
    
    # # Before writing I need to change the names of the coords to match the conventions in wavespectra.core.attributes
    # ds_efth = ds_efth.rename({
    #     'date': 'time',
    #     'longitude': 'lon',
    #     'latitude': 'lat'
    # })
    # ds_efth['time'].encoding['unlimited_dims'] = ('time',)
    
    # # make a specDataset object to get access to the built in functions to write output
    # ds_spec = wavespectra.SpecDataset(ds_efth)
    
    # ds_spec.to_ww3('test.nc')

if __name__ == '__main__':
    
    file_in='/home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/CMEMS_WAV/eez/2010_01.nc'
    grid_dir='/home/gfearon/code/somisana-ww3/configs/sa_west_02/GRID/'
    output_dir='/home/gfearon/code/somisana-ww3/configs/sa_west_02/SPEC_CMEMS/2010_01/'
    make_bry_spec_cmems(grid_dir+'lon.dat',grid_dir+'lat.dat',grid_dir+'mask.dat',file_in,output_dir)
    
    
    file_out='/mnt/c/Users/GilesF/Downloads/CMEMS_WAV/eez/test.nc'
    # cmems_2_spec(file_in,14,-34,file_out)
    # ncep_hind_2_spec()
