import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sys, os
sys.path.append(os.path.abspath("/xdisk/sylvia/nawdex-hackathon/shared/"))
import dict_nawdexsims
# simulations dictionary
simdict = dict_nawdexsims.simdictionary()
# dictionary for colors
colordict = dict_nawdexsims.colordictionary() 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
            spine.set_smart_bounds(True)


def select_analysis_days(ds, expid):
    
    import sys
    sys.path.append('/xdisk/sylvia/nawdex-hackathon/shared')

    import dict_nawdexsims
    simdict     = dict_nawdexsims.simdictionary()
    anadaysdict = dict_nawdexsims.anadaysdictionary()
    
    startday  = simdict[expid]['start']
    anadays   = anadaysdict[startday]
    
    return ds.sel(time=slice(anadays[0], anadays[-1]))

def seldays():
    # Selecting the same days as ICON simulations
    seldays = ['2016-09-21', '2016-09-22',
               '2016-09-23', '2016-09-24', '2016-09-25',
               '2016-09-30', '2016-10-01', '2016-10-02',
               '2016-10-03', '2016-10-04', '2016-10-05',
               '2016-10-14', '2016-10-15', '2016-10-16']

    # Period of simulations where all ICON setups at different resolution exist
    seldays_2 = ['2016-09-23', '2016-09-24', '2016-09-25']

    # 5 year climatology for september and october
    seldays_sep = []
    for year in ['2012','2013','2014','2015','2016']:
        for i in range(1,31,1): 
            for day in [str(i)]:
                if i in np.arange(1,10,1):
                   seldays_sep.append(year+'-09-0'+day)   
                else:
                   seldays_sep.append(year+'-09-'+day)
                
    seldays_oct = []
    for year in ['2012','2013','2014','2015','2016']:
        for i in range(1,32,1): 
            for day in [str(i)]:
                if i in np.arange(1,10,1):
                   seldays_oct.append(year+'-10-0'+day)   
                else:
                   seldays_oct.append(year+'-10-'+day)
                
    seldays_all = seldays_sep + seldays_oct
    return seldays_all

def sexy_axes( ax, fs ):
    ax.spines['left'].set_bounds(5,15)
    #ax.spines['bottom'].set_bounds(-0.3,0.3)
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    #ax.set_xlabel(r'Net cloud-radiative heating [K day$^{-1}$]',fontsize=fs)
    ax.axvline(x=0, ymin=0.0, ymax=1,c='black', lw=1)
    #ax.set_xlim(-0.3, 0.3)
    ax.set_ylim(5,15)
    ax.yaxis.set_ticks([5,7,9,11,13,15])
    ax.yaxis.set_ticklabels([5,7,9,11,13,15], fontsize=fs)#,fontweight='bold')
    #ax.xaxis.set_ticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    #ax.xaxis.set_ticklabels([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3], fontsize=fs)#,fontweight='bold')
    ax.set_xticklabels(ax.get_xticks(), rotation=45)
    ax.set_yticklabels(ax.get_yticks(), rotation=45)


def process_sim( dataset_attribute ):
    # Options for ERA5 CRH 
    if dataset_attribute == 'ERA5_sep_1':
        pcolor = 'k'
        conv = 1
        mphys = 1 
        
    elif dataset_attribute == 'ERA5_oct_1':
        pcolor = 'k'
        conv = 1
        mphys = 1
         
    elif dataset_attribute == 'ERA5_sep_2':
        pcolor = 'k'
        conv = 1
        mphys = 2 
        
    elif dataset_attribute == 'ERA5_oct_2':
        pcolor = 'k'
        conv = 1
        mphys = 2
        
    # options for CloudSat data
    elif dataset_attribute == 'CC_1':
        pcolor = 'k'
        conv = 1
        mphys = 1
            
    elif dataset_attribute == 'CC_2':
        pcolor = 'k'
        conv = 1
        mphys = 2    
            
    # ICON_simulations
    else:
        # get plotting color according to ICON resolution
        pcolor = colordict[simdict[dataset_attribute]['res']]
        # get linestyle according to convection scheme
        conv = simdict[dataset_attribute]['conv']
        if conv==0:
           lstyle='--'
        elif conv==1:
           lstyle='-'
        elif conv==2:
           lstyle=':'
            
        mphys = simdict[dataset_attribute]['mphys']

    return pcolor, conv, mphys


def get_fulllevel_height():
    # define simulation
    resolution = '80km'
    sim = '0001'
    expid = 'nawdexnwp-' + resolution + '-mis-' + sim
    # read ocean mask
    ipath_oceanmask = '/xdisk/sylvia/nawdex-hackathon/domain-mean/'
    da_ocean = xr.open_dataset(ipath_oceanmask + '/openoceanmask/' + expid + \
                               '_openoceanmask.nc')['mask_openocean']
    index = np.where(da_ocean == 1)[0]
    del da_ocean, ipath_oceanmask
    # read z_ifc data
    ipath = '/xdisk/sylvia/nawdex-hackathon/domain-mean/'
    ds = xr.open_dataset(ipath + 'nawdexnwp-' + resolution + '-mis-' + \
                         sim + '_2016092200_fg_DOM01_ML_0036.nc')
    del ipath
    # apply ocean mask on z_ifc to make sure to look at ocean grid point
    ds = ds.isel(ncells=index)
    
    # calculate full levels based on z_ifc
    z_full = (ds.z_ifc[:,0] - (ds.z_ifc[:,0].diff('height_3')/2)).values
    
    # get half level height
    z_half = ds.z_ifc[:, 0].values
    
    del resolution, sim, expid, index, ds
    return z_full, z_half

    
def time_merge_list(sim_set,resolution):
    
    # empty list for diffrent setups
    ds_merge1 = []
    ds_merge2 = []
    ds_merge3 = []
    ds_merge4 = []
    ds_merge5 = []
    ds_merge6 = []
       
    for i in range(len(sim_set)):
        
        conv = simdict[sim_set[i].attrs['simulation']]['conv']
        mphys = simdict[sim_set[i].attrs['simulation']]['mphys']
        res = simdict[sim_set[i].attrs['simulation']]['res']
        
        if conv == 1 and mphys ==1 and res==resolution:
            
            ds_merge1.append(sim_set[i])
            ds1 = xr.concat(ds_merge1,dim='time')
            
        elif conv == 1 and mphys ==2 and res==resolution:
            
            ds_merge2.append(sim_set[i])
            ds2 = xr.concat(ds_merge2,dim='time')
            
        elif conv == 0 and mphys ==1 and res==resolution:
            
            ds_merge3.append(sim_set[i])
            ds3 = xr.concat(ds_merge3,dim='time')
         
        elif conv == 0 and mphys ==2 and res==resolution:
            
            ds_merge4.append(sim_set[i])
            ds4 = xr.concat(ds_merge4,dim='time')
            
        elif conv == 2 and mphys ==1 and res==resolution:
            
            ds_merge5.append(sim_set[i])
            ds5 = xr.concat(ds_merge5,dim='time')
            
        elif conv == 2 and mphys ==2 and res==resolution:
            
            ds_merge6.append(sim_set[i])
            ds6 = xr.concat(ds_merge6,dim='time')
            
    if resolution=='80km' or resolution=='40km' or resolution=='20km' or resolution=='10km' or resolution=='5km':
        
        return (ds1,ds2) # Only with convection on / 1 and 2 moment mphy
    
    else:
        
        return (ds3,ds4,ds5,ds6) #( for 2km resolution with conv:off/only shallow convection, 1 and 2 moment mphy)



# function for vertical profile plot
def plot_varmean(fig, _ds_list, _var, fs, lw):

    fig.tight_layout(pad=2.5)
    
    # empty lists for time average
    
    elist1 = []
    elist2 = []
    
    for ds in _ds_list:

        # make sure that _var is part of dataset
        if _var not in ds.keys():
            continue

            
        mphys = simdict[ds.attrs['simulation']]['mphys']
        
        if mphys==1:
            
            elist1.append(ds[_var].sel(height=slice(16,50)).mean('time').values)
        
        elif mphys==2:
            
            elist2.append(ds[_var].sel(height=slice(16,50)).mean('time').values)    

            
    ax = fig.add_subplot(1,2,1)
    
    # time average over 14 days
    line1 = np.mean(elist1,axis=0)
    zfull = get_fulllevel_height()
    
    if _var == 'lw_crh':
        plt.plot(line1*86400,zfull[15:50]/1e3,linewidth=lw,label='lw_crh',color=colordict['80km'])
    if _var == 'sw_crh':
        plt.plot(line1*86400,zfull[15:50]/1e3,linewidth=lw,label='sw_crh',color=colordict['80km'],linestyle='--')    
    if _var == 'ddt_temp_radlw':
        plt.plot(line1*86400,zfull[15:50]/1e3,linewidth=lw,label='lw_rad',color=colordict['40km'])
    if _var == 'ddt_temp_radsw':
        plt.plot(line1*86400,zfull[15:50]/1e3,linewidth=lw,label='sw_rad',color=colordict['40km'],linestyle='--'    )
    if _var == 'ddt_temp_dyn2':
        plt.plot(line1*86400,zfull[15:50]/1e3,linewidth=lw,label='dyn',color=colordict['20km'])
    if _var == 'ddt_temp_turb':
        plt.plot(line1*86400,zfull[15:50]/1e3,linewidth=lw,label='turb',color=colordict['10km'])
    if _var == 'ddt_temp_pconv':
        plt.plot(line1*86400,zfull[15:50]/1e3,linewidth=lw,label='conv',color=colordict['5km'])
    if _var == 'ddt_temp_mphy':
        plt.plot(line1*86400,zfull[15:50]/1e3,linewidth=lw,label='mphy',color=colordict['2km'])    
        
    plt.tick_params(labelsize=13)
    #plt.ylabel("Height (km)",fontsize=12)
            
    ax.spines['left'].set_bounds(5,15)
    ax.spines['bottom'].set_bounds(-2.5,2.5)#(-0.25,0.25)
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.axvline(x=0, ymin=0.0, ymax=1,c='black', lw=1)
    plt.xlim(-2.5,2.5)#(-0.25,0.25)
    plt.ylim(5,15)
    ax.yaxis.set_ticks([5,7,9,11,13,15])
    ax.yaxis.set_ticklabels([5,7,9,11,13,15], fontsize=fs)#,fontweight='bold')
    ax.xaxis.set_ticks([-2,-1,0,1,2])
    ax.xaxis.set_ticklabels([-2.0,-1.0,0.0,1.0,2.0], fontsize=fs)#,fontweight='bold')
    plt.xticks(rotation=45); plt.yticks(rotation=45)
     
    ax.text(0.02,1.05,'(a)',weight='bold',fontsize=fs+2,transform=ax.transAxes)
    ax.text(0.02,0.95,'1-moment',fontsize=fs,transform=ax.transAxes)
           
    #plt.title('Microphysics: One-moment scheme',fontsize=fs, pad=20)
    plt.xlabel('Heating rates (K/day)',fontsize=fs)
    plt.legend(fontsize=fs-3,frameon=False,loc='upper right',bbox_to_anchor=(1.05,0.98))
    plt.ylabel("Height (km)",fontsize=fs)
    
    
    # mphy = 2
            
    ax = fig.add_subplot(1,2,2)
    
    # time average over 14 days
    line2 = np.mean(elist2,axis=0)
    
    if _var == 'lw_crh':
        plt.plot(line2*86400,zfull[15:50]/1e3,linewidth=lw,label='lw_crh',color=colordict['80km'])
    if _var == 'sw_crh':
        plt.plot(line2*86400,zfull[15:50]/1e3,linewidth=lw,label='sw_crh',color=colordict['80km'],linestyle='--')    
    if _var == 'ddt_temp_radlw':
        plt.plot(line2*86400,zfull[15:50]/1e3,linewidth=lw,label='lw_rad',color=colordict['40km'])
    if _var == 'ddt_temp_radsw':
        plt.plot(line2*86400,zfull[15:50]/1e3,linewidth=lw,label='sw_rad',color=colordict['40km'],linestyle='--')    
    if _var == 'ddt_temp_dyn2':
        plt.plot(line2*86400,zfull[15:50]/1e3,linewidth=lw,label='dyn',color=colordict['20km'])
    if _var == 'ddt_temp_turb':
        plt.plot(line2*86400,zfull[15:50]/1e3,linewidth=lw,label='turb',color=colordict['10km'])
    if _var == 'ddt_temp_pconv':
        plt.plot(line2*86400,zfull[15:50]/1e3,linewidth=lw,label='conv',color=colordict['5km'])
    if _var == 'ddt_temp_mphy':
        plt.plot(line2*86400,zfull[15:50]/1e3,linewidth=lw,label='mphy',color=colordict['2km']) 
        
    plt.tick_params(labelsize=13)
    #plt.ylabel("Height (km)",fontsize=12)
            
    ax.spines['left'].set_bounds(5,15)
    ax.spines['bottom'].set_bounds(-2.5,2.5)#(-0.25,0.25)
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.axvline(x=0, ymin=0.0, ymax=1,c='black', lw=1)
    plt.xlim(-2.5,2.5)#(-0.25,0.25)
    plt.ylim(5,15)
    ax.yaxis.set_ticks([5,7,9,11,13,15])
    ax.yaxis.set_ticklabels([5,7,9,11,13,15], fontsize=fs)#,fontweight='bold')
    ax.xaxis.set_ticks([-2,-1,0,1,2])
    ax.xaxis.set_ticklabels([-2.0,-1.0,0.0,1.0,2.0], fontsize=fs)#,fontweight='bold')
    plt.xticks(rotation=45); plt.yticks(rotation=45)
                
    ax.text(0.02,1.05,'(b)',weight='bold',fontsize=fs+2,transform=ax.transAxes)
    ax.text(0.02,0.95,'2-moment',fontsize=fs,transform=ax.transAxes)   
 
    #plt.title('Microphysics: Two-moment scheme',fontsize=fs, pad=20)
    plt.xlabel('Heating rates (K/day)',fontsize=fs)

