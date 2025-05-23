import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.animation import FuncAnimation
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import ww3_tools.postprocess as post

class LandmaskFeature(cfeature.GSHHSFeature):
    """from the OpenDrift code"""
    def __init__(self, scale='auto', globe=None, **kwargs):
        super().__init__(scale, **kwargs)

        if globe is not None:
            self._crs = ccrs.PlateCarree(globe=globe)

    def geometries(self):
        self.intersecting_geometries(extent=None)

    def intersecting_geometries(self, extent):
        global __polys__

        if self._scale == 'auto':
            scale = self._scale_from_extent(extent)
        else:
            scale = self._scale[0]
        return super().intersecting_geometries(extent)

def plot_land(ax, ocean_color = 'white', land_color = cfeature.COLORS['land'], lscale = 'auto', globe=None):
    """
    Plot the landmask or the shapes from GSHHG.
    (from the OpenDrift code)
    lscale = resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
    """
    land = LandmaskFeature(scale=lscale, facecolor=land_color, globe=globe)

    ax.add_feature(land, zorder=2,
                   facecolor=land_color,
                   edgecolor='black')

def setup_plot(ax, lon=None, lat=None, extents=None, land_color=('k', 0), lscale='h'):
    '''
    generic stuff applicable to all 2D plots
    
    lscale = resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto')
    '''
    # extents = [lon_min, lon_max, lat_min, lat_max]
    #
    # first need to get the domain extents if it's not set autmatically
    if extents is None:
        lon_min = min(np.ravel(lon))
        lon_max = max(np.ravel(lon))
        lat_min = min(np.ravel(lat))
        lat_max = max(np.ravel(lat))
        factor=0.05 # factor of domain size used to get dl
        dl = 0.5 * (lon_max - lon_min + lat_max - lat_min) * factor
        extents=[lon_min-dl,lon_max+dl,lat_min-dl,lat_max+dl]
    
    ax.set_extent(extents)
    plot_land(ax,land_color=land_color,lscale=lscale)
    # ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='dimgrey', alpha=0.5, linestyle=':')
    gl.right_labels = False
    gl.top_labels = False
    
    return extents

def plot_var(ax,var_data,lon,lat,
             ticks = [], # the ticks to plot
             cmap = 'Spectral_r'
             ):
    '''
    Add a variable to a 2D plot
    '''
    
    # set up the cmap to handle non-uniform input ticks
    if len(ticks)==0:
        ticks = np.linspace(min(np.ravel(var_data)),max(np.ravel(var_data)),num=20)
    levs = np.array(ticks)
    cmap_norm = mplc.BoundaryNorm(boundaries=levs, ncolors=256)
    
    # plot the data
    var_plt = ax.pcolormesh(lon,
                              lat,
                              var_data,
                              cmap=cmap,
                              norm=cmap_norm,
                              transform=ccrs.PlateCarree())
    
    return var_plt

def plot_cbar(ax,var_plt,
             ticks=[],
             tick_font = 13,
             label='values',
             label_font=15,
             loc=None, # [left, bottom, width, height]
             aspect_ratio=1,
             orientation='vertical'):
    
    '''
    Add a colorbar to a plot
    '''
    if loc is None:
        # this can be hard coded because we took care to set up the figsize and
        # axis location to accomodate the colorbar
        x_position = 0.8
        x_thickness = 0.015
        loc = [x_position, 0.2, x_thickness, 0.6]
    
    cbarax = plt.gcf().add_axes(loc) 
    
    cbar_plt = plt.colorbar(var_plt, cbarax,
                        ticks=ticks,
                        orientation=orientation)
    cbar_plt.set_label(label, fontsize=label_font)
    cbar_plt.ax.tick_params(labelsize=tick_font)
    
    return cbar_plt

def plot_time(ax,time,
             loc=[0.5,1.01],
             tstep=0,
             ref_date = datetime(2000, 1, 1, 0, 0, 0),
             time_fmt = '%Y-%m-%d %H:%M',
             time_font=15):
    '''
    Add time text to a 2D plot
    '''
    
    time_plt = ax.text(loc[0], loc[1],  pd.Timestamp(time).strftime(time_fmt),#datetime.strftime(time, time_fmt),
        ha='center', fontsize=time_font,
        transform=ax.transAxes)
    
    return time_plt

def plot_uv(ax,u,v,lon,lat,
              extents = None,
              skip_uv = 5,
              scale = 10,
              ref_vector = None,
              col = 'k'
              ):
    '''
    Add vectors to a 2D plot
    '''
        
    # plot the data
    width=0.0035  # Explicitly set the vector width for consistency between data and reference vector
    uv_plt = ax.quiver(lon[::skip_uv, ::skip_uv],
                      lat[::skip_uv, ::skip_uv],
                      u[::skip_uv, ::skip_uv],
                      v[::skip_uv, ::skip_uv],
                      scale=scale,
                      color=col,
                      width=width,
                      transform=ccrs.PlateCarree(), zorder=1)
    
    if ref_vector is not None:
        # add a reference vector
        d_y=extents[3]-extents[2]
        d_x=extents[1]-extents[0]
        loc_x=extents[0]+d_x*0.04
        loc_y=extents[3]-d_y*0.04
        ax.quiver(np.array([loc_x]), np.array([loc_y]),
                  np.array([ref_vector]), np.array([0]), 
                  scale=scale,
                  color=col,
                  width=width,
                  transform=ccrs.PlateCarree(), zorder=1,
                  )
        # add the label for it
        # try to define the y location of the text by shifting it down by a defined fraction of the y-range in the data
        loc_y_txt=loc_y-d_y*0.04
        ax.text(loc_x, loc_y_txt,  str(ref_vector)+' m',
            ha='left', 
            #fontsize=time_font,
            transform=ccrs.PlateCarree())
    return uv_plt

def get_uv_params(max_hs,scale_uv,ref_vector,skip_uv,num_vectors,aspect_ratio):
    # utility function to dynamically define params for scaling vectors    
    
    # Define thresholds for scaling and reference vector size
    # based on the maximum hs being plotted
    if max_hs < 1.5:
        scale = 10
        ref = 0.5
    elif max_hs < 2.5:
        scale = 20
        ref = 1
    elif max_hs < 4:
        scale = 30
        ref = 2
    elif max_hs < 5:
        scale = 40
        ref = 3
    else:
        scale = 60
        ref = 4
    
    # only update if they haven't already been specified by the user
    if scale_uv is None:
        scale_uv =scale
    if ref_vector is None:
        ref_vector =ref
        
    return scale_uv, skip_uv, ref_vector

def plot(fname,
        ax=None, # allowing for adding to an existing axis
        var='hs', # variable to plot
        time=slice(None), # If a single value, then a plot is made, if a slice, then an animation between those times is made
        ticks = None, #np.linspace(12,22,num=11), (gets set automatically if None)
        cmap = 'Spectral_r',
        extents = None, # [lon0,lon1,lat0,lat1] whole domain plotted if None
        add_cbar = True, # add a colorbar?
        cbar_loc = None, # [left, bottom, width, height] (gets set automatically if None)
        cbar_label = None, # 'temperature ($\degree$C)', we just use 'var' is None
        add_vectors = True, # add horizontal vectors?
        dir_var = 'dir', # allow for specifying other direction parameters like dp
        scale_uv = None, # define the horizontal vector scaling (gets set automatically if None)
        ref_vector = None, # value of reference vector (gets set automatically if None)
        skip_uv = None, # only every nth vector will get plotted (automatically defined in None)
        num_vectors=25, # baseline number of vectors in each direction, given an aspect ratio of 1
        skip_time = 1, # every nth time-step will be animated (if provided)
        add_time_label = True,
        isobaths = None, # optional list of isobaths to overlay over plot
        jpg_out=None, # full path to jpg output
        gif_out=None, # full path to gif output
        mp4_out=None, # option to rather write an mp4
        ):
    '''
    this is a convenience function for doing a quick 2D plot with minimal coding.
    this might also be used as example code for doing your own plots 
    there's also an option to turn the plot into an animation
    '''
    
    # get the data we want to 
    print('extracting the data to plot')    
    ds = post.get_ds_map(fname,time=time)
    lon = ds.longitude.values
    lat = ds.latitude.values
    
    # get the variable dataarray
    time_var=np.atleast_1d(ds.time.values)
    da_var=ds[var].squeeze()
    
    if len(time_var)==1:
        data_plt=da_var.values
    else:
        # this will be an animation, starting with the first time-step
        data_plt=da_var.isel(time=0).values
    
    if ticks is None:
        # get the range of the data to plot (using 5th and 95th percentiles)
        vmin=np.nanpercentile(da_var, 1)
        vmax=np.nanpercentile(da_var, 99)
        # round these to two significant figures
        vmin=round(vmin, 2 - int(np.floor(np.log10(abs(vmin)))) - 1)
        vmax=round(vmax, 2 - int(np.floor(np.log10(abs(vmax)))) - 1)
        num_ticks = 10
        step = (vmax - vmin) / num_ticks
        step = round(step, 2 - int(np.floor(np.log10(abs(step)))) - 1)
        # update vmax based on the rounded step
        vmax = vmin + num_ticks * step
        # Generate the ticks using the rounded step size
        ticks = np.arange(vmin, vmax + step/10, step) # Add a small value to ensure new_vmax is included
    
    # compute the extents from the grid if not explicitly defined
    if extents is None:
        lon_min = min(np.ravel(lon))
        lon_max = max(np.ravel(lon))
        lat_min = min(np.ravel(lat))
        lat_max = max(np.ravel(lat))
        factor=0.05 # factor of domain size used to get dl
        dl = 0.5 * (lon_max - lon_min + lat_max - lat_min) * factor
        extents=[lon_min-dl,lon_max+dl,lat_min-dl,lat_max+dl]
    
    # Create a Mercator projection
    proj = ccrs.Mercator()
    # Convert the corner points of the extents to Mercator projected coordinates
    # this is needed to compute the plot aspect ratio properly
    x0, y0 = proj.transform_point(extents[0], extents[2], ccrs.PlateCarree())  # lon0, lat0
    x1, y1 = proj.transform_point(extents[1], extents[3], ccrs.PlateCarree())  # lon1, lat1
    # Calculate the width and height of the domain in projected coordinates
    width = abs(x1 - x0)
    height = abs(y1 - y0)
    aspect_ratio = width/height
    
    if ax is None:        
        # set figsize according to the plot aspect ratio
        if aspect_ratio>1:
            fig_width = 6*aspect_ratio
            fig_height = 6
        else:
            fig_width = 6
            fig_height = 6/aspect_ratio
        # cbar_ax_width = 0.2 * fig_width if add_cbar else 0
        buffer_left=0.1 * fig_width
        buffer_right=0.2 * fig_width if add_cbar else 0.1 * fig_width
        figsize = (buffer_left + fig_width + buffer_right, fig_height)
        
        fig = plt.figure(figsize=figsize) 
        # [left, bottom, width, height] in fractions of figure dimensions
        width=0.7 if add_cbar else 0.8
        ax = fig.add_axes([0.1, 0.1, width, 0.8], projection=proj)
    
    # set up the plot
    extents = setup_plot(ax, lon, lat, extents = extents)
    
    # plot the data
    var_plt = plot_var(ax,data_plt,lon,lat, 
             ticks=ticks,
             cmap=cmap
             )
    
    # add a time label
    if add_time_label:
        time_plt = plot_time(ax,time_var[0])
    
    if add_cbar:
        if cbar_label is None:
            cbar_label = var
        plot_cbar(ax,var_plt,label=cbar_label,ticks=ticks,loc=cbar_loc,aspect_ratio=aspect_ratio)
    
    if add_vectors:
        
        print('getting the u/v vectors')
        
        # get vector components from Hs and wave direction
        da_hs=ds['hs'].squeeze()
        da_dir=ds[dir_var].squeeze()
        # Convert direction-from to direction-to (in degrees)
        da_dir = (da_dir + 180) % 360
        # Convert degrees to radians
        da_dir = np.deg2rad(da_dir)
        # Compute vector components (direction *towards*)
        da_u = da_hs * np.sin(da_dir)  # x-component (eastward)
        da_v = da_hs * np.cos(da_dir)  # y-component (northward)
        
        if len(time_var)==1:
            u_plt=da_u.values
            v_plt=da_v.values
        else:
            # this will be an animation, starting with the first time-step
            u_plt=da_u.isel(time=0).values
            v_plt=da_v.isel(time=0).values
        
        if skip_uv is None:
            # dynamically define how many vectors we want to plot
            [Ny,Nx]=np.shape(u_plt)
            skip_uv=int(Ny/num_vectors*aspect_ratio)
            
        hs_max=np.nanmax(da_hs.values.ravel())
        
        scale_uv, skip_uv, ref_vector = get_uv_params(hs_max,scale_uv,ref_vector,skip_uv,num_vectors,aspect_ratio)
        
        uv_plt = plot_uv(ax,u_plt,v_plt,lon,lat,
                      scale = scale_uv,
                      skip_uv = skip_uv,
                      ref_vector = ref_vector,
                      extents=extents
                      )
    
    # write a jpg if specified
    if len(time_var)==1: # single time-step specified, so a plot, not an animation
        if jpg_out is not None:
            print('writing '+jpg_out)
            plt.savefig(jpg_out,dpi=500,bbox_inches = 'tight')
            
        return ax
    
    else: # do the animation
        def plot_tstep(i):
            # get the data for this time-step
            var_i=da_var.isel(time=i).values
            
            # update the figure for this time-step
            if add_time_label:
                time_plt.set_text(pd.Timestamp(time_var[i]).strftime('%Y-%m-%d %H:%M'))
                var_plt.set_array(var_i.ravel())
            
            if add_vectors:
                u_i=da_u.isel(time=i).values
                v_i=da_v.isel(time=i).values
                uv_plt.set_UVC(u_i[::skip_uv, ::skip_uv],
                                        v_i[::skip_uv, ::skip_uv])
        
        # animate
        print('making animation')
        tstep_end=len(time_var)
        
        anim = FuncAnimation(
            fig, plot_tstep, frames=range(0,tstep_end,skip_time)) 
        
        if gif_out is not None:
            print('writing '+gif_out)
            anim.save(gif_out, writer='imagemagick')
        if mp4_out is not None:
            print('writing '+mp4_out)
            anim.save(mp4_out, writer="ffmpeg")

if __name__ == "__main__":
    
    fname='/home/gfearon/code/somisana-ww3/configs/sa_west_02/RUN_01/output/ww3.20101*.nc'
    plot(fname,time=1200)
    
        
