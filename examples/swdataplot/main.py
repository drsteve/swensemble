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
from scipy.integrate import simps, trapz
from swdatanal import getDistrib, omniDataCorr
from getswdata import getOMNIfiles, dataClean, dateShift, dateList

startDate   = datetime.datetime(2014,1,1,0,0,0)
endDate     = datetime.datetime(2014,12,31,23,59,59)
locDateList = dateList(startDate, endDate, shift = 'month')
OMNIfnames  = getOMNIfiles(locDateList,'/home/ehab/SWData','hourly')

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

SWV = dataClean(SWV,[2500.0],['>='])
SWN = dataClean(SWN,[999.0,0.0],['>=','<'])
SWT = dataClean(SWT,[1.0e6],['>='])

sEpoch = datetime.datetime(2014, 1, 1, 0, 0)
eEpoch = datetime.datetime(2014,12,31,23,59)
sEpochYYID = bisect.bisect_left(epochs, sEpoch)
eEpochYYID = bisect.bisect_left(epochs, eEpoch)

sEpoch = datetime.datetime(2014, 1, 1, 0, 0)
eEpoch = datetime.datetime(2014, 6,30,23,59)
sEpochJJID = bisect.bisect_left(epochs, sEpoch)
eEpochJJID = bisect.bisect_left(epochs, eEpoch)

sEpoch = datetime.datetime(2014, 7, 1, 0, 0)
eEpoch = datetime.datetime(2014,12,31,23,59)
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

SWV = filter(lambda v: v==v, SWV)
SWN = filter(lambda v: v==v, SWN)
SWT = filter(lambda v: v==v, SWT)


plt.figure(2)
Vstride = 5
plt.subplot(2, 3, 1)
SWVHST, bins = getDistrib(SWV[sEpochYYID:eEpochYYID], stride = Vstride, norm=True)
plt.plot(bins,SWVHST,label='Jan-Dec')
SWVHST = getDistrib(SWV[sEpochJJID:eEpochJJID], bins=bins, norm=True)
plt.plot(bins,SWVHST,label='Jan-Jun')
SWVHST = getDistrib(SWV[sEpochJDID:eEpochJDID], bins=bins, norm=True)
plt.plot(bins,SWVHST,label='Jul-Dec')
plt.title('Solar Wind Speed PDF')
plt.xlabel('Solar Wind Speed (Km/s) (stride = '+ str(Vstride) + ' Km/s)')
plt.ylabel('PDF')
plt.legend()

Nstride = 5
plt.subplot(2, 3, 2)
SWNHST, bins = getDistrib(SWN[sEpochYYID:eEpochYYID], stride = Nstride, norm=True)
plt.plot(bins,SWNHST,label='Jan-Dec')
SWNHST = getDistrib(SWN[sEpochJJID:eEpochJJID], bins=bins, norm=True)
plt.plot(bins,SWNHST,label='Jan-Jun')
SWNHST = getDistrib(SWN[sEpochJDID:eEpochJDID], bins=bins, norm=True)
plt.plot(bins,SWNHST,label='Jul-Dec')
plt.title('Solar Wind Density PDF')
plt.xlabel('Solar Wind Density (N/cc) (stride = '+ str(Nstride) + ' N/cc)')
plt.ylabel('PDF')
plt.legend()

Tstride = 2000
plt.subplot(2, 3, 3)
SWTHST, bins = getDistrib(SWT[sEpochYYID:eEpochYYID], stride = Tstride, norm=True)
plt.plot(bins,SWTHST,label='Jan-Dec')
SWTHST = getDistrib(SWT[sEpochJJID:eEpochJJID], bins=bins, norm=True)
plt.plot(bins,SWTHST,label='Jan-Jun')
SWTHST = getDistrib(SWT[sEpochJDID:eEpochJDID], bins=bins, norm=True)
plt.plot(bins,SWTHST,label='Jul-Dec')
plt.title('Solar Wind Temperature PDF')
plt.xlabel('Solar Wind Temperature ($K^o$) (stride = '+ str(Tstride) + ' $K^o$)')
plt.ylabel('PDF')
plt.legend()

plt.subplot(2, 3, 4)
stride = (max(SWBx[sEpochYYID:eEpochYYID])-min(SWBx[sEpochYYID:eEpochYYID]))/100
bins=numpy.arange(min(SWBx[sEpochYYID:eEpochYYID])-stride,max(SWBx[sEpochYYID:eEpochYYID])+stride,stride)
SWBHST = map(float, scipy.stats.histogram2(SWBx[sEpochYYID:eEpochYYID],bins))
SWBHST = [SWBHST[i]/sum(SWBHST) for i in range(len(SWBHST))]
plt.plot(bins,SWBHST,label='Jan-Dec')
SWBHST = map(float, scipy.stats.histogram2(SWBx[sEpochJJID:eEpochJJID],bins))
SWBHST = [SWBHST[i]/sum(SWBHST) for i in range(len(SWBHST))]
plt.plot(bins,SWBHST,label='Jan-Jun')
SWBHST = map(float, scipy.stats.histogram2(SWBx[sEpochJDID:eEpochJDID],bins))
SWBHST = [SWBHST[i]/sum(SWBHST) for i in range(len(SWBHST))]
plt.plot(bins,SWBHST,label='Jul-Dec')
plt.title('Solar Wind Bx_GSE PDF')
plt.xlabel('Solar Wind Bx_GSE (nT) (stride = '+ str(stride) + ' nT)')
plt.ylabel('PDF')
plt.legend()

plt.subplot(2, 3, 5)
stride = (max(SWBy[sEpochYYID:eEpochYYID])-min(SWBy[sEpochYYID:eEpochYYID]))/100
bins=numpy.arange(min(SWBy[sEpochYYID:eEpochYYID])-stride,max(SWBy[sEpochYYID:eEpochYYID])+stride,stride)
SWBHST = map(float, scipy.stats.histogram2(SWBy[sEpochYYID:eEpochYYID],bins))
SWBHST = [SWBHST[i]/sum(SWBHST) for i in range(len(SWBHST))]
plt.plot(bins,SWBHST,label='Jan-Dec')
SWBHST = map(float, scipy.stats.histogram2(SWBy[sEpochJJID:eEpochJJID],bins))
SWBHST = [SWBHST[i]/sum(SWBHST) for i in range(len(SWBHST))]
plt.plot(bins,SWBHST,label='Jan-Jun')
SWBHST = map(float, scipy.stats.histogram2(SWBy[sEpochJDID:eEpochJDID],bins))
SWBHST = [SWBHST[i]/sum(SWBHST) for i in range(len(SWBHST))]
plt.plot(bins,SWBHST,label='Jul-Dec')
plt.title('Solar Wind By_GSE PDF')
plt.xlabel('Solar Wind By_GSE (nT) (stride = '+ str(stride) + ' nT)')
plt.ylabel('PDF')
plt.legend()

plt.subplot(2, 3, 6)
stride = (max(SWBz[sEpochYYID:eEpochYYID])-min(SWBz[sEpochYYID:eEpochYYID]))/100
bins = numpy.arange(math.floor(min(SWBz[sEpochYYID:eEpochYYID])-stride),math.ceil(max(SWBz[sEpochYYID:eEpochYYID])+stride),math.ceil(stride))
SWBHST = map(float, scipy.stats.histogram2(SWBz[sEpochYYID:eEpochYYID],bins))
SWBHST = [SWBHST[i]/sum(SWBHST) for i in range(len(SWBHST))]
area = trapz(SWBHST, dx=1)
plt.plot(bins,SWBHST,label='Jan-Dec')
SWBHST = map(float, scipy.stats.histogram2(SWBz[sEpochJJID:eEpochJJID],bins))
SWBHST = [SWBHST[i]/sum(SWBHST) for i in range(len(SWBHST))]
area = trapz(SWBHST, dx=1)
plt.plot(bins,SWBHST,label='Jan-Jun')
SWBHST = map(float, scipy.stats.histogram2(SWBz[sEpochJDID:eEpochJDID],bins))
SWBHST = [SWBHST[i]/sum(SWBHST) for i in range(len(SWBHST))]
area = trapz(SWBHST, dx=1)
plt.plot(bins,SWBHST,label='Jul-Dec')
plt.title('Solar Wind Bz_GSE PDF')
plt.xlabel('Solar Wind Bz_GSE (nT) (stride = '+ str(stride) + ' nT)')
plt.ylabel('PDF')
plt.legend()


plt.suptitle('OMNI Hourly Data for 2014')
plt.show()


sys.exit()

