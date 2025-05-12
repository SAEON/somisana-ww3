'''
this serves as a command line interface (CLI) to execute functions from 
within this python repo directly from the command line.
The intended use is to allow python functions to be run from the cli docker image for this repo
The only functions I'm adding here are ones which produce an output e.g. a netcdf file
Feel free to add more functions from the repo as we need them in the cli
'''
import argparse
import os
import warnings
import xarray as xr
from datetime import datetime, timedelta
from ww3_tools.preprocess import make_bry_spec_cmems

# functions to help parsing string input to object types needed by python functions
def parse_datetime(value):
    try:
        return datetime.strptime(value, '%Y-%m-%d %H:%M:%S')
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid datetime format. Please use 'YYYY-MM-DD HH:MM:SS'.")

def parse_int(value):
    if value is None:
        return None
    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid integer value: {value}")

def parse_list(value):
    return [float(x) for x in value.split(',')]

def parse_list_str(value):
    if value is None or value == 'None':
        return None
    else:
        return [x.strip() for x in value.split(',')]
    
def parse_bool(s: str) -> bool:
    try:
        return {'true':True, 'false':False}[s.lower()]
    except KeyError:
        raise argparse.ArgumentTypeError(f'expect true/false, got: {s}')

def main():
    
    parser = argparse.ArgumentParser(description='Command-line interface for selected functions in the somisana-croco repo')
    subparsers = parser.add_subparsers(dest='function', help='Select the function to run')

    # just keep adding new subparsers for each new function as we go...
    
    # ----------------------
    # make_bry_cmems_monthy
    # ----------------------
    parser_make_bry_cmems_monthy = subparsers.add_parser('make_bry_cmems_monthy',
            help='Make monthly spectral boundary condition files for WW3 interannual runs')
    parser_make_bry_cmems_monthy.add_argument('--input_dir', required=True, type=str, 
            help='Path to directory containing the monthly files downloaded from CMEMS')
    parser_make_bry_cmems_monthy.add_argument('--lon_file', required=True, type=str,
            help='Path to longitude file in ww3 format')
    parser_make_bry_cmems_monthy.add_argument('--lat_file', required=True, type=str,
            help='Path to latitude file in ww3 format')
    parser_make_bry_cmems_monthy.add_argument('--mask_file', required=True, type=str,
            help='Path to mask file in ww3 format')
    parser_make_bry_cmems_monthy.add_argument('--output_dir', required=True, type=str,
            help='Path to where the forcing files will be saved')
    parser_make_bry_cmems_monthy.add_argument('--month_start', required=True, type=str, 
            help='first month in the interannual run in format "YYYY-MM"')
    parser_make_bry_cmems_monthy.add_argument('--month_end', required=True, type=str,
            help='last month in the interannual run in format "YYYY-MM"')
    def make_bry_cmems_monthy_handler(args):
        
        month_now = datetime.strptime(args.month_start+'-01','%Y-%m-%d')
        month_end = datetime.strptime(args.month_end+'-01','%Y-%m-%d')
        
        while month_now <= month_end:
            
            print('working on '+month_now.strftime('%Y-%m'))
            
            # create a directory for the spectral files for this month
            output_dir_month = os.path.join(args.output_dir,month_now.strftime('%Y_%m'))
            if not os.path.exists(output_dir_month):
                os.makedirs(output_dir_month)
            
            # define a list of 2 input files - last month and next month
            # this is to ensure we have forcing data over the last few hours of the month
            month_next = month_now + timedelta(days=32) # 32 days ensures we at least get to the next month
            fname_month_now = os.path.join(args.input_dir, month_now.strftime('%Y_%m.nc'))
            fname_month_next = os.path.join(args.input_dir, month_next.strftime('%Y_%m.nc'))
            # We could handle the case where the user doesn't have files for the previous and next month
            if not os.path.exists(fname_month_now):
                raise ValueError("Processing of "+month_now.strftime('%Y-%m')+" requires "+fname_month_now)
            if not os.path.exists(fname_month_next):
                raise ValueError("Processing of "+month_now.strftime('%Y-%m')+" requires "+fname_month_next)
            input_file=[fname_month_now, fname_month_next]
            
            cmems_ds = xr.open_mfdataset(input_file)
            cmems_ds = cmems_ds.sel(time=slice(month_now,datetime(month_next.year,month_next.month,1))) # we only need the first time-step of next month
            
            # make the boundary file for this month
            # Suppress some pesky warnings which were cluttering up the console 
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            warnings.filterwarnings("ignore", message=".*", category=UserWarning, module="wavespectra.output.ww3")
            make_bry_spec_cmems(args.lon_file,
                                    args.lat_file,
                                    args.mask_file,
                                    cmems_ds,
                                    output_dir_month,
                                    num_freqs=32,
                                    start_freq=0.0373,
                                    freq_factor=1.1,
                                    num_dirs=24,
                                    fp_2_dspr = 120, # factor to get directional spreading from peak frequency
                                    )
            cmems_ds.close()
            
            month_now=datetime(month_next.year, month_next.month, 1) # set month_now to the first day of the next month
        
    parser_make_bry_cmems_monthy.set_defaults(func=make_bry_cmems_monthy_handler)
    
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        print("Please specify a function.")

if __name__ == "__main__":
    main()
