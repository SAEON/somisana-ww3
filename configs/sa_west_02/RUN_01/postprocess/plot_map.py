from datetime import datetime,timedelta
import ww3_tools.plotting as ww3plt
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean.cm as cmo
import xarray as xr
from matplotlib.path import Path

time_plt=datetime(2010,5,10)
figsize=(12,8) # (hz,vt)
extents = [14,20.5,-37,-29]
cmap = 'Spectral_r'
ticks = np.linspace(0,7,num=8)
cbar_label = 'Hm0 (m)'
scale = 80

fig = plt.figure(figsize=figsize) 
ax = plt.axes(projection=ccrs.Mercator())

fname_cmems='/home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/CMEMS_WAV/eez/2010_05.nc'
ds_cmems=xr.open_dataset(fname_cmems)
ds_cmems=ds_cmems.sel(time=time_plt,method='nearest')
lon_cmems=ds_cmems.longitude.values
lat_cmems=ds_cmems.latitude.values
lon_cmems, lat_cmems = np.meshgrid(lon_cmems, lat_cmems)
hs_cmems = ds_cmems['VHM0'].values
dir_cmems = ds_cmems['VMDR'].values
# Convert direction-from to direction-to (in degrees)
dir_cmems = (dir_cmems + 180) % 360
# Convert degrees to radians
dir_cmems = np.deg2rad(dir_cmems)
# Compute vector components (direction *towards*)
u_cmems = hs_cmems * np.sin(dir_cmems)  # x-component (eastward)
v_cmems = hs_cmems * np.cos(dir_cmems)  # y-component (northward)

# plot the data
levs = np.array(ticks)
cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)
cmems_plt = ax.pcolormesh(lon_cmems,
                          lat_cmems,
                          hs_cmems,
                          cmap=cmap,
                          norm=cmap_norm,
                          transform=ccrs.PlateCarree())

# add vectors
width=0.002  # Explicitly set the vector width for consistency between data and reference vector
skip_uv=3
uv_plt = ax.quiver(lon_cmems[::skip_uv,::skip_uv],
                  lat_cmems[::skip_uv,::skip_uv],
                  u_cmems[::skip_uv, ::skip_uv],
                  v_cmems[::skip_uv, ::skip_uv],
                  scale=scale,
                  color='k',
                  width=width,
                  transform=ccrs.PlateCarree(), zorder=1)

ww3plt.setup_plot(ax,extents = extents,lscale='i')
# add the time
ax.text(0.5,0.95,  datetime.strftime(time_plt, '%Y-%m-%d %H:%M'),
    ha='center', fontsize=12,
    transform=ax.transAxes)

fname='../output/ww3.201005.nc'

# ww3plt.plot(fname,time=slice(datetime(2010,5,1),datetime(2010,5,15)),skip_time=12,
#             gif_out='ww3.201005.gif')


ww3plt.plot(fname,time=slice(time_plt,time_plt+timedelta(hours=0.5)),
            ax=ax,
            ticks=ticks,
            scale_uv=scale,
            )