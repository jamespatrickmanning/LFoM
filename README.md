# LFoM
JiM Aug 2020

Code to plot and process Nick's temp & dissolved oxygen files.

There are two routines: one plot time series plots (plt_nicks.py) and one to make a map (plt_LFoM_sites.py).
The first one creates a file "LFoM_sites.csv" which is read by the second.

Assumes you have:

-the .csv and .gps files in same directory

-installed Python with netCDF4 module to get NGDC depths

-have all the docks included in the "harborlist.txt" file

-examined the editting criteria like "testT", for example,  which filters out temp > that value
