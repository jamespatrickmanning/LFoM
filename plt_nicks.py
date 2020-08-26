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


inputdir='C:/Users/DELL/Downloads/lfa_dmf/'
case=['Kurt_Holmes_HOUND_DOG_','Bill_Chaprales_RUEBY_','Bill_Lister_MUSSEL_POINT_','Mike_Rego_MISS_LILLY_','Willy_Ogg_Jr_HAPPY_TRAILS_']
sn=[[2002042,2002043,2002044,2002045,2002046],[2006052,2006053,2006054,2006055,2006057],[2006066,2006067,2006068,2006069,2006072],[2006060,2006061,2006062,2006063,2006064],[2002047,2002048,2002049,2002050,2006051] ]
# testT criteria for max temperature 
testT=[19.0,19.0,20.0,22.0,19.0]
mindistfromharbor=0.2 # minimum dist from harbor acceptable
#####  IMPORT modules
import numpy as np
import netCDF4
import pandas as pd
import matplotlib.pyplot as plt
import os
#from conversions import dd2dm

###  FUNCTIONS
def dd2dm(lat,lon):
    """
    convert lat, lon from decimal degrees to degrees,minutes
    """
    lat_d = int(abs(lat))                #calculate latitude degrees
    lat_m = (abs(lat) - lat_d) * 60. #calculate latitude minutes

    lon_d = int(abs(lon))
    lon_m = (abs(lon) - lon_d) * 60.
    
    la=lat_d*100.+lat_m
    lo=lon_d*100.+lon_m
    return la,lo
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

def gps_compare_JiM(lat,lon,mindistfromharbor): #check to see if the boat is in the harbor derived from Huanxin's "wifipc.py" functions   
    # function returns yes if this position is within "mindistfromharbor" miles of a dock
    file='harborlist.txt' # has header line lat, lon, harbor
    df=pd.read_csv(file,sep=',')
    indice_lat=[i for i ,v in enumerate(abs(np.array(df['lat'])-lat)<mindistfromharbor) if v]
    indice_lon=[i for i ,v in enumerate(abs(np.array(df['lon'])-lon)<mindistfromharbor) if v]
    harbor_point_list=[i for i, j in zip(indice_lat,indice_lon) if i==j]
    if len(harbor_point_list)>0:
       near_harbor='yes'
    else:
       near_harbor='no'
    return near_harbor #yes or no
### MAIN CODE
incp=0
for f in range(len(case)):
 us=[i for i in range(len(case[f])) if case[f].startswith('_', i)]# gets index of underscores in "case" 
 fman=case[f][0:us[1]+1] # directory with csv & gps files
 if len(us)==5:
    vessel=case[f][us[0]+1:us[2]]
 else:
    vessel=case[f][us[0]+1:us[1]]
 for k in sn[f]: # loop through probe serial numbers
  plt.rcParams["figure.figsize"] = [16,9]
  direct=inputdir+case[f]+str(k)
  inc=0 #increments through sn
  for file in os.listdir(direct):
    if (file.startswith(str(k))) and (file.endswith("DissolvedOxygen.csv")):
        fnp=os.path.join(direct, file[0:27])+'.gps'
        try:
            dfp=open(fnp)
            line=dfp.readline()
            if (line[0:3]=='RWS') & (line[5]=='4'):
                lat=float(line[5:14])
                lon=float(line[14:24])
            else:
                line=dfp.readline()
                line=dfp.readline()
                line=dfp.readline()
                if (line[0:3]=='SWS') & (line[5]=='4'):
                    lat=float(line[5:14])
                    lon=float(line[14:24])
                else:
                    lat=np.nan
                    lon=np.nan
            dfp.close()
        except:
            lat=np.nan
            lon=np.nan
        if ~np.isnan(lat):
            [la,lo]=dd2dm(lat,lon)# converts to ddmm format change indexes 4/9/2019    
            yorn=gps_compare_JiM(la,lo,mindistfromharbor) # check to see if the data is from the dock (near harbor)
            if yorn=='y':
                    continue   
        incf=0 #increment through files with this sn
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
        df=df[df['Temperature (C)']<testT[f]] # remove data greater than "testT"   
        inc=inc+1
        df['lat']=lat
        df['lon']=lon        

        if inc==1: #increment through sn 
            dfall=df
        else:
            dfall= pd.concat([dfall, df], axis=0) # adds to whole dataframe for this participant
        incp=incp+1
        if incp==1:
            dfpout=pd.DataFrame([[fman,k,df['lat'].mean(),df['lon'].mean()]], columns=['fman','sn','lat','lon'])            
        else:
            dfp2out = pd.DataFrame([[fman,k,df['lat'].mean(),df['lon'].mean()]], columns=['fman','sn','lat','lon'])
            dfpout=pd.concat([dfpout, dfp2out])# saves all the position data
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
 plt.close('all')
dfpout.to_csv('LFoM_sites.csv')

