#!/bin/env python

from config import DATA_DIR, LATS, LONS, GRIDS, MONTHS_YEAR
from netCDF4 import Dataset
from local_utils import get_month_slice, get_grid_subset, days_year, get_day_of_year_index, is_leap_year 
import numpy as np

from mpi4py import MPI
import argparse
import sys

comm = MPI.COMM_WORLD
mpi_rank = comm.Get_rank()
mpi_size = comm.Get_size()
grids_size = len(GRIDS)

diff_mode = False

def main(year, ensemble, diff_filename):
    #If the wrong number of cluster nodes have been allocated, crash the program 
    if (grids_size % mpi_size) > 0:
        print 'Error: number of grid points not evenly divisible by allocated nodes'
        sys.exit(1)

    nc = Dataset(DATA_DIR + 'z500_{y}.nc'.format(y=year))
    diff = Dataset(diff_filename) if diff_mode else None

    #Initialize reduction products arrays
    daily_results_reduce = np.zeros((days_year(year), LATS, LONS))
    anom_results_reduce = np.zeros((days_year(year), LATS, LONS))
    
    #Do the mean calculations
    daily_results, anom_results = calc_mean(year, ensemble, nc, diff) 

    #Stitch the distributed results back together
    comm.Reduce(daily_results, daily_results_reduce, op=MPI.SUM, root=0)
    comm.Reduce(anom_results, anom_results_reduce, op=MPI.SUM, root=0)
    
    #Write out the results
    if mpi_rank == 0:
        write_nc_file(daily_results_reduce, 'z500_{y}_daily_mean_en{e}.nc'.format(y=year, e=ensemble), nc)

        if diff_mode:
            write_nc_file(anom_results_reduce, 'z500_{y}_daily_anom_en{e}.nc'.format(y=year, e=ensemble), nc, anom_mode=True)


def calc_mean(year, ensemble, nc, diff):
    time_var = nc.variables['time']
    z500_var = nc.variables['z500']    

    daily_results = np.zeros((days_year(year), LATS, LONS))
    anom_results = np.zeros((days_year(year), LATS, LONS))

    #Do some quick calculations to determine which subset of the array is ours
    start_index, end_index = get_grid_subset(grids_size, mpi_rank, mpi_size)

    for lat, lon in GRIDS[start_index:end_index]:
        for day in range(days_year(year)):
            day_idx = get_day_of_year_index(day + 1)

            daily_value = z500_var[day_idx, ensemble, lat, lon]
            daily_results[day, lat, lon] = daily_value 
            
            if diff:
                baseperiod_array = diff.variables['1981_2010_daily_mean_leap'] if is_leap_year(year) else diff.variables['1981_2010_daily_mean'] 
                anom_results[day, lat, lon] = daily_value - baseperiod_array[day, lat, lon]
    
    return daily_results, anom_results


#Write out results as netCDF
def write_nc_file(daily_results, filename, nc, anom_mode=False):
    #Grab every 4th time value to represent daily
    daily_time_var = nc.variables['time'][::4]

    nc_out = Dataset(filename, mode='w', format='NETCDF4') 
    nc_out.createDimension('lon', LONS)
    nc_out.createDimension('lat', LATS)
    nc_out.createDimension('time', None) #UNLIMITED
    nc_out.createDimension('month', MONTHS_YEAR)
    nc_out.title = ''
    nc_out.institution = ''
    nc_out.project = ''
    nc_out.contact = 'ken.kunkel@noaa.gov'
    nc_out.Conventions = "CF-1.6"
    
    longitude = nc_out.createVariable('lon', 'f8', ('lon',))
    longitude.standard_name = 'longitude'
    longitude.long_name = 'longitude'
    longitude.units = 'degrees_east'
    longitude.modulo = 360.0
    longitude.axis = 'X'
    longitude[:] = np.arange(0, 360.0, 2.0)
    
    latitude = nc_out.createVariable('lat', 'f8', ('lat',))
    latitude.standard_name = 'latitude'
    latitude.long_name = 'latitude'
    latitude.units = 'degrees_north'
    latitude.axis = 'Y'
    latitude[:] = np.arange(-90.0, 92.0, 2.0)
    
    time = nc_out.createVariable('time', 'f8', ('time',))
    time.units = 'hours since 1-1-1 0:0:0' 
    time.calendar = 'standard' #Gregorian
    time[:] = daily_time_var 
    
    if anom_mode:
        daily_mean = nc_out.createVariable('daily_anom', 'f8', ('time', 'lat', 'lon'))
        daily_mean.long_name = 'z500 daily anomaly vs 1981-2010'
    else:
        daily_mean = nc_out.createVariable('daily_mean', 'f8', ('time', 'lat', 'lon'))
        daily_mean.long_name = 'z500 daily mean'

    daily_mean[:] = daily_results
    nc_out.close()
    

#Setup command-line arguments parser options
parser = argparse.ArgumentParser()
parser.add_argument('year', help='Year to process', type=int)
parser.add_argument('ensemble', help='Ensemble index to process', type=int)
parser.add_argument('-d', '--diff', help='File that contains base period grid')
parser.add_argument('-o', '--outdir', help='Directory to which to write output files')
args = parser.parse_args()

#output_dir = args.outdir if args.outdir else '../dist/'

if args.diff:
    diff_mode = True    
main(args.year, args.ensemble, args.diff)
