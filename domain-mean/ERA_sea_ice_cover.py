# #This python script download sea_ice cover from ERA 5 archive
# Behrooz Keshtgar, KIT

import cdsapi
c = cdsapi.Client()

years = ['2012','2013','2014','2015','2016']

months = ['09','10']

days = ['01','02','03','04','05','06','07','08','09','10',
        '11','12','13','14','15','16','17','18','19','20',
        '21','22','23','24','25','26','27','28','29','30']

for year in years:

    for month in months:

        for day in days:

            sel_day = year+'-'+month+'-'+day

            c.retrieve('reanalysis-era5-single-levels', {
            'date': [sel_day],
            'param': '31',
            'product_type': 'reanalysis',
            'area': '80/-78/23/40',  #N/W/S/E
            'grid': '0.25/0.25',
            'format': 'netcdf'},
            'ERA5_sea_ice_'+sel_day)
#-----------------------------------------------------------
# Download last day of October

for year in years:

           da = year+'-10-31'

           c.retrieve('reanalysis-era5-single-levels', {
           'date': [da],
           'param': '31',
           'product_type': 'reanalysis',
           'area': '80/-78/23/40',  #N/W/S/E
           'grid': '0.25/0.25',
           'format': 'netcdf'},
           'ERA5_sea_ice_'+da)

