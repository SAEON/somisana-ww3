import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from datetime import datetime
from matplotlib import gridspec
 
fname_mod='ww3_CP01.nc'
ds_mod = xr.open_dataset(fname_mod)

fname_obs='/home/gfearon/wave_obs/processed/CP01.nc'
ds_obs=xr.open_dataset(fname_obs)

fname_cmems='/home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/CMEMS_WAV/eez/timeseries/cmems_CP01.nc'
ds_cmems=xr.open_dataset(fname_cmems)

start_time = datetime(2010,5,1)#ds_mod.time[0]
end_time = datetime(2010,6,1)#ds_mod.time[-1]

ds_mod = ds_mod.sel(time=slice(start_time,end_time))
ds_obs = ds_obs.sel(time=slice(start_time,end_time))
ds_cmems = ds_cmems.sel(time=slice(start_time,end_time))

hs_obs=ds_obs['Hmo']
tp_obs=ds_obs['Tp']
dir_obs=ds_obs['Directn']
hs_mod = ds_mod['hs']
tp_mod = 1/ds_mod['fp']
dir_mod = ds_mod['dir']
hs_cmems = ds_cmems['VHM0']
tp_cmems = ds_cmems['VTPK']
dir_cmems = ds_cmems['VMDR']

fig = plt.figure(figsize=(12,8),dpi=300, facecolor='white')

gs=gridspec.GridSpec(3,1)
gs.update(hspace=0.15, wspace=0.)

# all data
ax1 = fig.add_subplot(gs[0,:])
ax1.plot(ds_obs.time, hs_obs, label='observations', color='black', linestyle='-', linewidth=2)
ax1.plot(ds_mod.time, hs_mod, label='WW3', color='blue', linestyle='-', linewidth=2)
ax1.plot(ds_cmems.time, hs_cmems, label='CMEMS', color='red', linestyle='-', linewidth=2)
# Customize the plot
# ax1.set_xlabel('Time')
ax1.set_ylabel('Hm0 (m)')
ax1.legend(loc='best') # Let matplotlib decide the best location for the legend
ax1.grid(True, alpha=0.3)  # Add a subtle grid for readability
ax1.set_ylim(0,9.)
ax1.set_xlim(start_time,end_time)

ax2 = fig.add_subplot(gs[1,:])
ax2.plot(ds_obs.time, tp_obs, label='observations', color='black', linestyle='-', linewidth=2)
ax2.plot(ds_mod.time, tp_mod, label='WW3', color='blue', linestyle='-', linewidth=2)
ax2.plot(ds_cmems.time, tp_cmems, label='CMEMS', color='red', linestyle='-', linewidth=2)
# Customize the plot
# ax2.set_xlabel('Time')
ax2.set_ylabel('Tp (s)')
ax2.grid(True, alpha=0.3)  # Add a subtle grid for readability
ax2.set_ylim(4,20)
ax2.set_xlim(start_time,end_time)

ax3 = fig.add_subplot(gs[2,:])
ax3.plot(ds_obs.time, dir_obs, label='observations', color='black', linestyle='-', linewidth=2)
ax3.plot(ds_mod.time, dir_mod, label='WW3', color='blue', linestyle='-', linewidth=2)
ax3.plot(ds_cmems.time, dir_cmems, label='CMEMS', color='red', linestyle='-', linewidth=2)
# Customize the plot
ax3.set_xlabel('Time')
ax3.set_ylabel('Mean Dir ($\deg$)')
ax3.grid(True, alpha=0.3)  # Add a subtle grid for readability
ax3.set_ylim(120,360)
ax3.set_yticks(np.linspace(120,360,5))
ax3.set_xlim(start_time,end_time)

# jpg_out='plot_ww3_CP01.jpg'
# plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')



