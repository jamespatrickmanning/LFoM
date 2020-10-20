# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 14:22:41 2020

@author: JiM
Processing Nick's probe data with DO included
First application Cape Cod Bay 2020 (old version in plt_nicks_preOct2020.py)
Modified in October 2020 to read Nick's new format
"""
#####  HARDCODES 
#fn=direct+'2002042_low_20200716_050955_DissolvedOxygen_GPS.csv'
#fnp=direct+'2002042_low_20200716_050955.gps'
testO=0.3# criteria to first difference test Oxygen in water


inputdir='C:/Users/Joann/Downloads/lfm_dmf/'
#inputdir='C:/Users/DELL/Downloads/lfa_dmf/'
case=['Kurt_Holmes_HOUND_DOG_','Bill_Chaprales_RUEBY_','Bill_Lister_MUSSEL_POINT_','Mike_Rego_MISS_LILLY_','Willy_Ogg_Jr_HAPPY_TRAILS_']
sn=[[2002042,2002043,2002044,2002045,2002046],[2006052,2006053,2006054,2006055,2006057],[2006066,2006067,2006068,2006069,2006072],[2006060,2006061,2006062,2006063,2006064],[2002047,2002048,2002049,2002050,2006051] ]
# testT criteria for max temperature 
testT=[19.0,19.0,20.0,22.0,19.0]
mindistfromharbor=.5# minimum dist from harbor acceptable
Host = '66.114.154.52' # ISOMEDIA machine where emolt.org is
UserName = 'mingchaossh'
Pswd = 'eMOLT1$'
remot_dir = '/httpdocs'
local_folder = ''
#####  IMPORT modules
import numpy as np
import netCDF4
import pandas as pd
import matplotlib.pyplot as plt
import ftplib
import os
import warnings
warnings.filterwarnings('ignore')
#import multiple_models as mm
#params = {'axes.labelsize': 18,
#          'axes.titlesize': 18}
#plt.rcParams.update(params)
plt.rcParams["font.size"] = 18
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
def upload_emolt(Host, UserName, Pswd,remot_dir, local_folder):
    '''upload result to cloud'''
    ftp = ftplib.FTP(Host)
    ftp.login(UserName, Pswd)
    ftp.cwd(remot_dir)
    command = 'STOR emolt_nicks.dat'
    sendFILEName = 'emolt_nicks.dat'
    ftp.storbinary(command, open(sendFILEName,'rb'))
    print('uploading emolt_nicks.dat to emolt.org ')
    ftp.quit()
### MAIN CODE
if os.path.exists('emolt_nicks.csv'): # refresh the daily average output file
    os.remove('emolt_nicks.csv')
incp=0
for f in range(len(case)):
#for f in [0,1]: # use this to test one case only    
 us=[i for i in range(len(case[f])) if case[f].startswith('_', i)]# gets index of underscores in "case" 
 fman=case[f][0:us[1]+1] # directory with csv  files has fishers name
 if len(us)==5: # Accounts for names like "willy_Ogg_Jr" having multiple parts
    vessel=case[f][us[0]+1:us[2]]
 else:
    vessel=case[f][us[0]+1:us[1]]
 for k in sn[f]: # loop through probe serial numbers
 #for k in [sn[f][3]]:    
  plt.rcParams["figure.figsize"] = [16,9]
  direct=inputdir+case[f]+str(k)
  inc=0 #increments through sn
  for file in os.listdir(direct):
    if (file.startswith(str(k))) and (file.endswith("DissolvedOxygen_GPS.csv")):
        incf=0 #increment through files with this sn
        print(os.path.join(direct, file))
        fn=os.path.join(direct, file)
        ##### READ & PLOT
        df=pd.read_csv(fn)# gets raw probe data
        #df=df.drop('Dissolved Oxygen (%)',axis=1)# drop this since Nick says it is no good
        df=df.rename(columns={'TEMP':'Temperature (C)'}) # rename column
        df=df.rename(columns={'DO':'Dissolved Oxygen (mg/l)'})
        df=df.rename(columns={'S_lat':'lat'}) # rename column
        df=df.rename(columns={'S_long':'lon'})
        for kk in range(len(df)):
            if np.isnan(df['lat'][kk]):
                df['lat'][kk]=df['E_lat'][kk]
                df['lon'][kk]=df['E_long'][kk]
        if ~np.isnan(df['lat'][0]):
            [la,lo]=dd2dm(df['lat'][0],df['lon'][0])# converts to ddmm format change indexes 4/9/2019  
            yorn=gps_compare_JiM(la,lo,mindistfromharbor) # check to see if the data is from the dock (near harbor)
            if yorn=='yes':
                continue  
        if np.isnan(df['lat'][0]):
            continue
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
        df=df[df['lat']!=np.nan]
        #if (f==0) & (inc==1): #first case and first serial number 
        if inc==1: # first case of this serial number 
            dfall=df
        else:
            dfall= pd.concat([dfall, df], axis=0) # adds to whole dataframe for this participant
        incp=incp+1
        if incp==1:
            dfpout=pd.DataFrame([[fman,k,df['lat'].mean(),df['lon'].mean()]], columns=['fman','sn','lat','lon'])            
        else:
            dfp2out = pd.DataFrame([[fman,k,df['lat'].mean(),df['lon'].mean()]], columns=['fman','sn','lat','lon'])
            dfpout=pd.concat([dfpout, dfp2out])# saves all the position data
  # here's where we save an "emolt.dat" -like file which includes:
  # vessel_# sn mth day hr mn yd lon lat 1.0 nan depth range_depth days temp std_temp year
  #'''
  dfall['mth']=0
  dfall['day']=0
  dfall['yrday']=0.0#initializes as float  
  #dfall['doppio']=[]#model estimate
  for jj in range(len(dfall)):
      dfall['mth'][jj]=dfall.index.month[jj].astype(np.int64)
      dfall['day'][jj]=dfall.index.day[jj].astype(np.int64)
      dfall['yrday'][jj]=dfall.index[jj].timetuple().tm_yday+(dfall.index[jj].hour+dfall.index[jj].minute/60.)/24.
  dfall['std_temp']=dfall['Temperature (C)'].std()

  dfallda=dfall.resample('D').mean()
  dfallda['vessel_num']=50+f
  dfallda['serial_num']=k
  dfallda['hour']=[12]*len(dfallda)
  dfallda['min']=[0]*len(dfallda)
  dfallda['days']=[1.0]*len(dfallda)
  dfallda['year']=dfallda.index.year 
  dfallda['dum1']=[0]*len(dfallda)
  #depth=get_depth(dfallda['lon'],dfallda['lat'],0.4)
  try:
      depth=get_depth(dfall['lon'].mean(skipna=True),dfall['lat'].mean(skipna=True),0.4)
  except:
      depth=np.nan
  dfallda['depth']=[depth]*len(dfallda)
  dfallda['range_depth']=[0]*len(dfallda)
  use_col=['vessel_num','serial_num', 'mth', 'day', 'hour', 'min', 'yrday','lon','lat','Dissolved Oxygen (mg/l)','dum1','depth', 'range_depth', 'days', 'Temperature (C)','std_temp', 'year']
  #dfallda=dfallda[dfallda['lat']!=np.nan]
  dfallda=dfallda.dropna(axis=0,subset=['lat'])
  dfallda.to_csv('emolt_nicks.dat',sep=' ',columns=use_col,header=False,mode='a',index=False,float_format='%.3f')  
  #dfallda[[use_col]].to_csv('emolt_nicks.csv',sep=' ',header=False,mode='a',index=False,float_format='%.3f')  
  #'''    
  ax=dfall[['Dissolved Oxygen (mg/l)','Temperature (C)']].plot(subplots=True) # plots both columns
  FT=dfall['Temperature (C)'].values*1.8+32
  ax2=ax[1].twinx()
  ax2.set_ylim(ax[1].get_ylim()[0]*1.8+32,ax[1].get_ylim()[1]*1.8+32)
  ax2.plot(dfall.index,FT,color='r')
  ax2.set_ylabel('fahrenheit')
  #ax[0].set_title(fman+' SN='+str(k)+' @ '+"{0:.4g}".format(dfall['lat'].mean(skipna=True))+'N '+"{0:.4g}".format(dfall['lon'].mean(skipna=True))+'W', fontsize=12)
  #ax[0].set_title(fman+' SN='+str(k)+' @ ~'+"{0:.4g}".format(dfall['lat'].mean(skipna=True))+'N '+"{0:.4g}".format(dfall['lon'].mean(skipna=True))+'W in ~'+"{0:.2g}".format(abs(depth))+'meters', fontsize=12)
  ax[0].set_title(fman+' SN='+str(k)+' in ~'+"{0:.2g}".format(abs(3.28*depth))+' feet', fontsize=24)
  plt.savefig(str(k)+'_'+fman+'_time_series_plot_raw.png')
 plt.show()
 plt.close('all')
dfpout.to_csv('LFoM_sites.csv')
#upload_emolt(Host, UserName, Pswd,remot_dir, local_folder)

