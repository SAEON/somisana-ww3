from datetime import datetime
import ww3_tools.postprocess as post

fname='../output/ww3.201005.nc'

# Cape Point
lon_ts=18.28667
lat_ts=-34.204
post.get_ts(fname,lon_ts=lon_ts,lat_ts=lat_ts,nc_out='ww3_CP01.nc')

# Gordons Bay
lon_ts=18.84
lat_ts=-34.159
post.get_ts(fname,lon_ts=lon_ts,lat_ts=lat_ts,nc_out='ww3_GB01.nc')