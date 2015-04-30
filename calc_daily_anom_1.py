#!/bin/env python

from netCDF4 import Dataset, num2date
from local_utils import get_month_slice, get_day_slice
import numpy as np

DATA_DIR = '/snfs4/assessment/brooke/20Century/z500_Members/'
LATS = 91
LONS = 180
GRIDS = [(lat, lon) for lat in range(LATS) for lon in range(LONS)]
ENSEMBLE = 0

for year in range(2012, 2013):
    nc = Dataset(DATA_DIR + 'z500_{y}.nc'.format(y=year))
    time_var = nc.variables['time']
    z500_var = nc.variables['z500']    
    
    results = np.zeros((time_var.size, LATS, LONS), 
            dtype=[('daily_mean','f8'), ('daily_anom', 'f8')])

    for lat, lon in GRIDS:
        for day_idx, date in enumerate(dates_in_month):
            day_slicer = get_day_slice(time_var, year, month + 1, date.day)
            daily_mean = z500_var[day_slicer, ENSEMBLE, lat, lon].mean()
            daily_anom = daily_mean - monthly_mean

            results['daily_mean'][day_idx, lat, lon] = daily_mean
            results['daily_anom'][day_idx, lat, lon] = daily_anom

    #Write out results as netCDF
    with Dataset('z500_daily_anom_{y}'.format(y=year), mode='w', format='NETCDF4') as nc_out:
        nc_out.createDimension('lon', nc.variables['lon'].size)
        nc_out.createDimension('lat', nc.variables['lat'].size)
        nc_out.createDimension('time', time_var.size)
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
        longitude[:] = nc.variables['lon']
    
        latitude = nc_out.createVariable('lat', 'f8', ('lat',))
        latitude.standard_name = 'latitude'
        latitude.long_name = 'latitude'
        latitude.units = 'degrees_north'
        latitude.axis = 'Y'
        latitude[:] = nc.variables['lat']
    
        time = nc_out.createVariable('time', 'f8', ('time',))
        time.units = 'hours since 1-1-1 0:0:0' 
        time.calendar = 'standard' #Gregorian
    
        daily_mean = nc_out.createVariable('daily_mean', 'f8', ('time', 'lat', 'lon'))
        daily_mean[:] = results['daily_mean']
    
        daily_anom = nc_out.createVariable('daily_anom', 'f8', ('time', 'lat', 'lon'))
        daily_anom[:] = results['daily_anom']

