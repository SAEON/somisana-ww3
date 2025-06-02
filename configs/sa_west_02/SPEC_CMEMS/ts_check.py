import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from datetime import datetime
from matplotlib import gridspec
from wavespectra import read_ww3

fname_bry=['2010_04/cmems.lon14.40.lat-29.80.spec.nc','2010_05/cmems.lon14.40.lat-29.80.spec.nc']
ds_bry=read_ww3(fname_bry)

fname_cmems='/home/gfearon/code/somisana-croco/DATASETS_CROCOTOOLS/CMEMS_WAV/eez/timeseries/cmems.lon14.40.lat-29.80.nc'
ds_cmems=xr.open_dataset(fname_cmems)

start_time = datetime(2010,4,1)#ds_mod.time[0]
end_time = datetime(2010,6,1)#ds_mod.time[-1]

ds_bry = ds_bry.sel(time=slice(start_time,end_time))
ds_cmems = ds_cmems.sel(time=slice(start_time,end_time))

hs_bry=ds_bry.spec.hs()
tp_bry=ds_bry.spec.tp()
dir_bry=ds_bry.spec.dm()
hs_cmems = ds_cmems['VHM0']
tp_cmems = ds_cmems['VTPK']
dir_cmems = ds_cmems['VMDR']

fig = plt.figure(figsize=(12,8),dpi=300, facecolor='white')

gs=gridspec.GridSpec(3,1)
gs.update(hspace=0.15, wspace=0.)

# all data
ax1 = fig.add_subplot(gs[0,:])
ax1.plot(ds_cmems.time, hs_cmems, label='CMEMS', color='blue', linestyle='-', linewidth=2)
ax1.plot(ds_bry.time, hs_bry, label='bry', color='red', linestyle='-', linewidth=1)
# Customize the plot
# ax1.set_xlabel('Time')
ax1.set_ylabel('Hm0 (m)')
ax1.legend(loc='best') # Let matplotlib decide the best location for the legend
ax1.grid(True, alpha=0.3)  # Add a subtle grid for readability
ax1.set_ylim(0,9.)
ax1.set_xlim(start_time,end_time)

ax2 = fig.add_subplot(gs[1,:])
ax2.plot(ds_cmems.time, tp_cmems, label='CMEMS', color='blue', linestyle='-', linewidth=2)
ax2.plot(ds_bry.time, tp_bry, label='bry', color='red', linestyle='-', linewidth=1)
# Customize the plot
# ax2.set_xlabel('Time')
ax2.set_ylabel('Tp (s)')
ax2.grid(True, alpha=0.3)  # Add a subtle grid for readability
ax2.set_ylim(4,20)
ax2.set_xlim(start_time,end_time)

ax3 = fig.add_subplot(gs[2,:])
ax3.plot(ds_cmems.time, dir_cmems, label='CMEMS', color='blue', linestyle='-', linewidth=2)
ax3.plot(ds_bry.time, dir_bry, label='bry', color='red', linestyle='-', linewidth=1)
# Customize the plot
ax3.set_xlabel('Time')
ax3.set_ylabel('Mean Dir ($\deg$)')
ax3.grid(True, alpha=0.3)  # Add a subtle grid for readability
ax3.set_ylim(120,360)
ax3.set_yticks(np.linspace(120,360,5))
ax3.set_xlim(start_time,end_time)

# jpg_out='plot_ww3_CP01.jpg'
# plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')



