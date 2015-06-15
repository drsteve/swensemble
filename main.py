import os
import sys
sys.path.append('/home/ehab/MyFiles/Softex/spacePy/spacepy-0.1.5')
import math
import numpy
import scipy
import bisect
import scipy.stats
import datetime
from datetime import date
from spacepy import pycdf
from spacepy import seapy
from itertools import chain
import spacepy.time as spt
import matplotlib.pyplot as plt
from getswdata import getOMNIfiles, dataClean

OMNIfnames = getOMNIfiles(['2005'],'/home/ehab/SWData','hourly')
cdfData = []
for i in range(len(OMNIfnames)):
 cdf = pycdf.CDF(OMNIfnames[i])
 cdfData.append(cdf.copy())
cdf.close()

epochs = []
SWV = []; SWN = []; SWT = []
SWBx= []; SWBy= []; SWBz= []
for i in range(len(OMNIfnames)): 
 epochs = epochs + list(cdfData[i]['Epoch'])
 SWV    = SWV + list(cdfData[i]['V'])
 SWN    = SWN + list(cdfData[i]['N'])
 SWT    = SWT + list(cdfData[i]['T'])
 SWBx   = SWBx+ list(cdfData[i]['BX_GSE'])
 SWBy   = SWBy+ list(cdfData[i]['BY_GSE'])
 SWBz   = SWBz+ list(cdfData[i]['BZ_GSE'])

SWV = dataClean(SWV,2500.0)
SWN = dataClean(SWN,999.0)
SWT = dataClean(SWT,1.0e6)

sEpoch = datetime.datetime(2005, 1, 1, 0, 0)
eEpoch = datetime.datetime(2005,12,31,23,59)
sEpochYYID = bisect.bisect_left(epochs, sEpoch)
eEpochYYID = bisect.bisect_left(epochs, eEpoch)

sEpoch = datetime.datetime(2005, 1, 1, 0, 0)
eEpoch = datetime.datetime(2005, 6,30,23,59)
sEpochJJID = bisect.bisect_left(epochs, sEpoch)
eEpochJJID = bisect.bisect_left(epochs, eEpoch)

sEpoch = datetime.datetime(2005, 7, 1, 0, 0)
eEpoch = datetime.datetime(2005,12,31,23,59)
sEpochJDID = bisect.bisect_left(epochs, sEpoch)
eEpochJDID = bisect.bisect_left(epochs, eEpoch)

plt.figure(1)
time = epochs[sEpochYYID:eEpochYYID]
data = SWV[sEpochYYID:eEpochYYID]
plt.plot(time,data)
time = epochs[sEpochJJID:eEpochJJID]
data = SWV[sEpochJJID:eEpochJJID]
plt.plot(time,data)
time = epochs[sEpochJDID:eEpochJDID]
data = SWV[sEpochJDID:eEpochJDID]
plt.plot(time,data)
plt.show()

plt.figure(2)
plt.subplot(2, 3, 1)
stride = (max(SWV[sEpochYYID:eEpochYYID])-min(SWV[sEpochYYID:eEpochYYID]))/200
bins=numpy.arange(min(SWV[sEpochYYID:eEpochYYID])-stride,max(SWV[sEpochYYID:eEpochYYID])+stride,stride)
SWVHST = scipy.stats.histogram2(SWV[sEpochYYID:eEpochYYID],bins)
plt.plot(bins,SWVHST,label='Jan-Dec')
SWVHST = scipy.stats.histogram2(SWV[sEpochJJID:eEpochJJID],bins)
plt.plot(bins,SWVHST,label='Jan-Jun')
SWVHST = scipy.stats.histogram2(SWV[sEpochJDID:eEpochJDID],bins)
plt.plot(bins,SWVHST,label='Jul-Dec')
plt.title('Solar Wind Speed PDF')
plt.xlabel('Solar Wind Speed (Km/s) (stride = '+ str(stride) + ' Km/s)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 2)
stride = (max(SWN[sEpochYYID:eEpochYYID])-min(SWN[sEpochYYID:eEpochYYID]))/200
bins=numpy.arange(min(SWN[sEpochYYID:eEpochYYID])-stride,max(SWN[sEpochYYID:eEpochYYID])+stride,stride)
SWNHST = scipy.stats.histogram2(SWN[sEpochYYID:eEpochYYID],bins)
plt.plot(bins,SWNHST,label='Jan-Dec')
SWNHST = scipy.stats.histogram2(SWN[sEpochJJID:eEpochJJID],bins)
plt.plot(bins,SWNHST,label='Jan-Jun')
SWNHST = scipy.stats.histogram2(SWN[sEpochJDID:eEpochJDID],bins)
plt.plot(bins,SWNHST,label='Jul-Dec')
plt.title('Solar Wind Density PDF')
plt.xlabel('Solar Wind Density (N/cc) (stride = '+ str(stride) + ' N/cc)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 3)
stride = (max(SWT[sEpochYYID:eEpochYYID])-min(SWT[sEpochYYID:eEpochYYID]))/1000
bins=numpy.arange(min(SWT[sEpochYYID:eEpochYYID])-stride,max(SWT[sEpochYYID:eEpochYYID])+stride,stride)
SWTHST = scipy.stats.histogram2(SWT[sEpochYYID:eEpochYYID],bins)
plt.plot(bins,SWTHST,label='Jan-Dec')
SWTHST = scipy.stats.histogram2(SWT[sEpochJJID:eEpochJJID],bins)
plt.plot(bins,SWTHST,label='Jan-Jun')
SWTHST = scipy.stats.histogram2(SWT[sEpochJDID:eEpochJDID],bins)
plt.plot(bins,SWTHST,label='Jul-Dec')
plt.title('Solar Wind Temperature PDF')
plt.xlabel('Solar Wind Temperature ($K^o$) (stride = '+ str(stride) + ' $K^o$)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 4)
stride = (max(SWBx[sEpochYYID:eEpochYYID])-min(SWBx[sEpochYYID:eEpochYYID]))/100
bins=numpy.arange(min(SWBx[sEpochYYID:eEpochYYID])-stride,max(SWBx[sEpochYYID:eEpochYYID])+stride,stride)
SWBHST = scipy.stats.histogram2(SWBx[sEpochYYID:eEpochYYID],bins)
plt.plot(bins,SWBHST,label='Jan-Dec')
SWBHST = scipy.stats.histogram2(SWBx[sEpochJJID:eEpochJJID],bins)
plt.plot(bins,SWBHST,label='Jan-Jun')
SWBHST = scipy.stats.histogram2(SWBx[sEpochJDID:eEpochJDID],bins)
plt.plot(bins,SWBHST,label='Jul-Dec')
plt.title('Solar Wind Bx_GSE PDF')
plt.xlabel('Solar Wind Bx_GSE (nT) (stride = '+ str(stride) + ' nT)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 5)
stride = (max(SWBy[sEpochYYID:eEpochYYID])-min(SWBy[sEpochYYID:eEpochYYID]))/100
bins=numpy.arange(min(SWBy[sEpochYYID:eEpochYYID])-stride,max(SWBy[sEpochYYID:eEpochYYID])+stride,stride)
SWBHST = scipy.stats.histogram2(SWBy[sEpochYYID:eEpochYYID],bins)
plt.plot(bins,SWBHST,label='Jan-Dec')
SWBHST = scipy.stats.histogram2(SWBy[sEpochJJID:eEpochJJID],bins)
plt.plot(bins,SWBHST,label='Jan-Jun')
SWBHST = scipy.stats.histogram2(SWBy[sEpochJDID:eEpochJDID],bins)
plt.plot(bins,SWBHST,label='Jul-Dec')
plt.title('Solar Wind By_GSE PDF')
plt.xlabel('Solar Wind By_GSE (nT) (stride = '+ str(stride) + ' nT)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 6)
stride = (max(SWBz[sEpochYYID:eEpochYYID])-min(SWBz[sEpochYYID:eEpochYYID]))/100
bins=numpy.arange(min(SWBz[sEpochYYID:eEpochYYID])-stride,max(SWBz[sEpochYYID:eEpochYYID])+stride,stride)
SWBHST = scipy.stats.histogram2(SWBz[sEpochYYID:eEpochYYID],bins)
plt.plot(bins,SWBHST,label='Jan-Dec')
SWBHST = scipy.stats.histogram2(SWBz[sEpochJJID:eEpochJJID],bins)
plt.plot(bins,SWBHST,label='Jan-Jun')
SWBHST = scipy.stats.histogram2(SWBz[sEpochJDID:eEpochJDID],bins)
plt.plot(bins,SWBHST,label='Jul-Dec')
plt.title('Solar Wind Bz_GSE PDF')
plt.xlabel('Solar Wind Bz_GSE (nT) (stride = '+ str(stride) + ' nT)')
plt.ylabel('Number of Occurance')
plt.legend()


plt.suptitle('OMNI Hourly Data for 2014')
plt.show()


sys.exit()

