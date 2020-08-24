# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 14:22:41 2020

@author: JiM
Processing Nick's probe data with DO included
First application Cape Cod Bay 2020
"""
#####  HARDCODES 
#fn=direct+'2002042_low_20200716_050955_DissolvedOxygen.csv'
#fnp=direct+'2002042_low_20200716_050955.gps'
testO=0.3# criteria to first difference test Oxygen in water
testT=21.0# criteria for max temperature 

inputdir='C:/Users/DELL/Downloads/lfa_dmf/'
#case='Kurt_Holmes_HOUND_DOG_'
#sn=list(range(2002042,2002047,1)) # probe serial #s
#testT=19.0
#case='Bill_Chaprales_RUEBY_'
#sn=[2006052,2006053,2006054,2006055,2006057] # probe serial #s
#testT=19.0
#case='Bill_Lister_MUSSEL_POINT_'
#sn=[2006066,2006067,2006068,2006069,2006072] # probe serial #s
#case='Mike_Rego_MISS_LILLY_'
#sn=[2006060,2006061,2006062,2006063,2006064] # probe serial #s
case='Willy_Ogg_Jr_HAPPY_TRAILS_'
sn=[2002047,2002048,2002049,2002050,2006051] # probe serial #s
#testT=21.0


#####  IMPORT modules
import numpy as np
import netCDF4
import pandas as pd
import matplotlib.pyplot as plt
import os

###  FUNCTIONS
def get_depth(loni,lati,mindist_allowed):
    # routine to get depth (meters) using vol1 from NGDC
    url='https://www.ngdc.noaa.gov/thredds/dodsC/crm/crm_vol1.nc'
    nc = netCDF4.Dataset(url).variables 
    lon=nc['x'][:]
    lat=nc['y'][:]
    xi,yi,min_dist= nearlonlat_zl(lon,lat,loni,lati) 
    if min_dist>mindist_allowed:
      depth=np.nan
    else:
      depth=nc['z'][yi,xi]
    return depth#,min_dist

def nearlonlat_zl(lon,lat,lonp,latp): # needed for the next function get_FVCOM_bottom_temp 
    """ 
    used in "get_depth" to find distance
    """ 
    # approximation for small distance 
    cp=np.cos(latp*np.pi/180.) 
    dx=(lon-lonp)*cp
    dy=lat-latp 
    xi=np.argmin(abs(dx)) 
    yi=np.argmin(abs(dy))
    min_dist=111*np.sqrt(dx[xi]**2+dy[yi]**2)
    return xi,yi,min_dist

### MAIN CODE
us=[i for i in range(len(case)) if case.startswith('_', i)]# gets index of underscores in "case" 
fman=case[0:us[1]+1] # directory with csv & gps files
if len(us)==5:
    vessel=case[us[0]+1:us[2]]
else:
    vessel=case[us[0]+1:us[1]]

for k in sn: # loop through probe serial numbers
  direct=inputdir+case+str(k)
  inc=0
  for file in os.listdir(direct):
    if (file.startswith(str(k))) and (file.endswith("DissolvedOxygen.csv")):
        print(os.path.join(direct, file))
        fn=os.path.join(direct, file)
        ##### READ & PLOT
        df=pd.read_csv(fn)# gets raw probe data
        df=df.drop('Dissolved Oxygen (%)',axis=1)# drop this since Nick says it is no good
        df=df.rename(columns={'DO Temperature (C)':'Temperature (C)'}) # rename column
        df['date']=pd.to_datetime(df['ISO 8601 Time'], format='%Y-%m-%dT%H:%M:%S') # form a datetime
        df=df.set_index('date') # make this the "index" column of the dataframe
        diffO=df['Dissolved Oxygen (mg/l)'].diff()# first differences of Oxygen
        diffT=df['Temperature (C)'].diff()# first differences of Temperature
        # find "in water" estimate where DO levels off in first 15 records 
        id=np.where(abs(diffO[0:15])>testO)
        if len(id[0])>1:
            df=df[id[-1][-1]+1:-1]
        df=df[df['Temperature (C)']<testT]    
        #df.plot(subplots=True) # plots both columns
        fnp=os.path.join(direct, file[0:27])+'.gps'
        #dfp=pd.read_csv(fnp,skiprows=3,header=None,sep='\s=')# gets raw probe data
        try:
            dfp=open(fnp)
            line=dfp.readline()
            if (line[0:3]=='RWS') & (line[5]=='4'):
                df['lat']=float(line[5:14])
                df['lon']=float(line[14:24])
            else:
                line=dfp.readline()
                line=dfp.readline()
                line=dfp.readline()
                if (line[0:3]=='SWS') & (line[5]=='4'):
                    df['lat']=float(line[5:14])
                    df['lon']=float(line[14:24])
                else:
                    df['lat']=np.nan
                    df['lon']=np.nan
            dfp.close()
        except:
            df['lat']=np.nan
            df['lon']=np.nan            
        inc=inc+1
        if inc==1:
            dfall=df
        else:
            dfall= pd.concat([dfall, df], axis=0) # adds to whole dataframe for this participant
  ax=dfall[['Dissolved Oxygen (mg/l)','Temperature (C)']].plot(subplots=True) # plots both columns
  FT=dfall['Temperature (C)'].values*1.8+32
  ax2=ax[1].twinx()
  ax2.set_ylim(ax[1].get_ylim()[0]*1.8+32,ax[1].get_ylim()[1]*1.8+32)
  ax2.plot(dfall.index,FT,color='r')
  ax2.set_ylabel('fahrenheit')
  depth=get_depth(dfall['lon'].mean(skipna=True),dfall['lat'].mean(skipna=True),0.4)
  #ax[0].set_title(fman+' SN='+str(k)+' @ '+"{0:.4g}".format(dfall['lat'].mean(skipna=True))+'N '+"{0:.4g}".format(dfall['lon'].mean(skipna=True))+'W', fontsize=12)
  ax[0].set_title(fman+' SN='+str(k)+' @ ~'+"{0:.4g}".format(dfall['lat'].mean(skipna=True))+'N '+"{0:.4g}".format(dfall['lon'].mean(skipna=True))+'W in ~'+"{0:.2g}".format(abs(depth))+'meters', fontsize=12)
  plt.savefig(str(k)+'_'+fman+'_time_series_plot_raw.png')

plt.show()

