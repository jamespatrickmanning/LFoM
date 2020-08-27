# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 11:26:50 2020

@author: JiM
routine to plot positions on a basemap
specifically for the LFoM sites in Summer 2020
borrowed code from "getsst.py"
May be problems with "basemap" if not installed on your machine.
"""
###### IMPORT MODULES
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#NOTE:  JiM NEEDED THE FOLLOWING LINE TO POINT TO his PROJ LIBRARY
import os
os.environ['PROJ_LIB'] = 'c:\\Users\\Joann\\anaconda3\\pkgs\\proj4-5.2.0-ha925a31_1\\Library\share'
from mpl_toolkits.basemap import Basemap
###### HARDCODES
area='inside_CCBAY'
mode='hide_detail'#detail
def getgbox(area):
  # gets geographic box based on area
  if area=='SNE':
    gbox=[-71.,-66.,39.,42.] # for SNE
  elif area=='OOI':
    gbox=[-72.,-69.5,39.5,41.5] # for OOI
  elif area=='GBANK':
    gbox=[-70.,-64.,39.,42.] # for GBANK
  elif area=='GS':           
    gbox=[-71.,-63.,38.,42.5] # for Gulf Stream
  elif area=='NorthShore':
    gbox=[-71.,-69.5,41.5,43.] # for north shore
  elif area=='WNERR':
    gbox=[-71.,-70.,42.5,43.3] # for WNERR deployment
  elif area=='DESPASEATO':
    gbox=[-71.,-69.5,42.6,43.25] # for miniboat Despaseato deployment
  elif area=='CCBAY':
    gbox=[-70.75,-69.8,41.5,42.23] # CCBAY
  elif area=='inside_CCBAY':
    gbox=[-70.7,-70.,41.7,42.1] # inside_CCBAY
  elif area=='NEC':
    gbox=[-69.,-64.,39.,43.5] # NE Channel
  elif area=='NE':
    gbox=[-76.,-66.,35.,44.5] # NE Shelf 
  return gbox
gbox=getgbox(area) # uses the getgbox function to define lat/lon boundary
latsize=[gbox[2],gbox[3]]
lonsize=[gbox[0],gbox[1]]
tick_int=(gbox[3]-gbox[2])/4. # allow for 3-4 tick axis label intervals
if tick_int>2:
    tick_int=int(tick_int)   # make the tick_interval integer increments
if tick_int<=2:
    tick_int=.3
fig = plt.figure(figsize=(12,10))
m = Basemap(projection='merc',llcrnrlat=min(latsize),urcrnrlat=max(latsize),\
            llcrnrlon=min(lonsize),urcrnrlon=max(lonsize),resolution='f')
m.fillcontinents(color='gray')
df=pd.read_csv('LFoM_sites.csv').dropna()
fman=np.unique(df['fman'])
col=['r','g','b','c','m']
sym=['+','*','o','^','v']
kinc=0#increment of fishermen
for k in fman:# loop through fishermen
    kinc=kinc+1
    dffman=df[df['fman']==k]
    if mode=='detail':# this includes actual positions
        xfm,yfm=m(dffman['lon'].max(),dffman['lat'].mean())
        plt.annotate(k[0:-1],xy=(xfm,yfm),color=col[kinc-1],va="bottom",ha="left",fontsize=20,fontweight='bold')
        jinc=0
        for j in np.unique(dffman['sn']):#loop through probe sn
            jinc=jinc+1
            dfsn=dffman[dffman['sn']==j]# make dataframe for just this probe sn
            x,y=m(dfsn['lon'].values,dfsn['lat'].values) # convert to basemap coordinates
            m.scatter(x,y,color=col[kinc-1],marker=sym[jinc-1])
            xsn,ysn=m(dfsn['lon'].min(),dfsn['lat'].mean())# find mean pos for this probe
            l2d=str(dfsn['sn'].values[0]%100).zfill(2)
            plt.annotate(l2d,xy=(xsn,ysn),color=col[kinc-1],ha="right",va='bottom',fontsize=16)
    else:#'hide_detail'
        xfm,yfm=m(dffman['lon'].mean(),dffman['lat'].mean())  
        plt.annotate(k[0:-1],xy=(xfm,yfm),color=col[kinc-1],va="bottom",ha="left",fontsize=20,fontweight='bold')

m.drawparallels(np.arange(min(latsize),max(latsize)+1,tick_int),labels=[1,0,0,0])
m.drawmeridians(np.arange(min(lonsize),max(lonsize)+1,tick_int),labels=[0,0,0,1])
#m.drawcoastlines()
m.drawmapboundary()
plt.title('LFoM sites 2020',fontsize=24)
plt.savefig('LFoM_sites_'+mode+'.png')
plt.show()