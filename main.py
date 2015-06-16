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

startDate = datetime.datetime(2008,7,1)
endDate   = datetime.datetime(2008,9,1 )
refDate   = datetime.datetime(2014,7,1)
locDateList = dateList(startDate, endDate, shift = 'month')
OMNIfnames = getOMNIfiles(locDateList,'/home/ehab/SWData','hourly')
if refDate in locDateList:
 pass
else:
 OMNIfnames = OMNIfnames + getOMNIfiles([refDate],'/home/ehab/SWData','hourly')

cdfData = []
for i in range(len(OMNIfnames)):
 cdf = pycdf.CDF(OMNIfnames[i])
 cdfData.append(cdf.copy())
cdf.close()

epochs = []; SWP = []
for i in range(len(OMNIfnames)): 
 epochs = epochs + list(cdfData[i]['Epoch'])
 SWP    = SWP + list(cdfData[i]['V'])
SWP = dataClean(SWP,[2500],['>='])

SWPDateRng, cepochs, KSVals, KSDist, aepochs = omniDataCorr(refDate, startDate, endDate, epochs, SWP, CorrTime = 'Day')

fig = plt.figure()
ax = plt.subplot2grid((2,2), (0,0), colspan=2)
plt.plot(cepochs,SWPDateRng)
plt.title('Solar Wind Speed')
plt.ylabel('Speed (km/s')

ax = plt.subplot2grid((2,2),(1, 0))
plt.plot(aepochs,KSVals[:,0],label='KS-Stat')
plt.plot(aepochs,KSVals[:,1],label='P-Value')
plt.title('Kolmogorov-Smirnov test between Velocity Samples')
plt.xlabel('Month')
plt.ylabel('KS-Stat and P-Value')
ax.set_yscale('log')

ax = plt.subplot2grid((2,2),(1, 1))
plt.plot(aepochs,KSDist[:,0],label='KS-Stat')
plt.plot(aepochs,KSDist[:,1],label='P-Value')
plt.title('Kolmogorov-Smirnov test between Velcotity Distributions')
plt.xlabel('Month')
plt.ylabel('KS-Stat and P-Value')
plt.legend(loc='best')

plt.suptitle('Compare the Velocity Daily Variation Relative to First Date')
plt.show()
sys.exit()





















epochs = []
SWV = []; SWN = []; SWT = []
SWBx= []; SWBy= []; SWBz= []
for i in range(len(OMNIfnames)): 
 epochs = epochs + list(cdfData[i]['Epoch'])
 SWV    = SWV + list(cdfData[i]['flow_speed'])
#SWV    = SWV + list(cdfData[i]['proton_density'])
#SWV    = SWV + list(cdfData[i]['V'])
#SWN    = SWN + list(cdfData[i]['N'])
#SWT    = SWT + list(cdfData[i]['T'])
#SWBx   = SWBx+ list(cdfData[i]['BX_GSE'])
#SWBy   = SWBy+ list(cdfData[i]['BY_GSE'])
#SWBz   = SWBz+ list(cdfData[i]['BZ_GSE'])

SWV = dataClean(SWV,[2500.0],['>='])
#SWN = dataClean(SWN,[999.0,0.0],['>=','<'])
#SWT = dataClean(SWT,[1.0e6],['>='])

sys.exit()

OMNIfnames = getOMNIfiles(['2014'],'/home/ehab/SWData','hourly')
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

#plt.figure(1)
#time = epochs[sEpochYYID:eEpochYYID]
#data = SWV[sEpochYYID:eEpochYYID]
#plt.plot(time,data)
#time = epochs[sEpochJJID:eEpochJJID]
#data = SWV[sEpochJJID:eEpochJJID]
#plt.plot(time,data)
#time = epochs[sEpochJDID:eEpochJDID]
#data = SWV[sEpochJDID:eEpochJDID]
#plt.plot(time,data)
#plt.show()

filter(lambda v: v==v, SWV)
filter(lambda v: v==v, SWN)
filter(lambda v: v==v, SWT)

plt.figure(2)
plt.subplot(2, 3, 1)
SWVHST, bins, stride = getDistrib(SWV[sEpochYYID:eEpochYYID], nbins=200, norm=True)
plt.plot(bins,SWVHST,label='Jan-Dec')
SWVHST = getDistrib(SWV[sEpochJJID:eEpochJJID], bins=bins, norm=True)
plt.plot(bins,SWVHST,label='Jan-Jun')
SWVHST = getDistrib(SWV[sEpochJDID:eEpochJDID], bins=bins, norm=True)
plt.plot(bins,SWVHST,label='Jul-Dec')
plt.title('Solar Wind Speed PDF')
plt.xlabel('Solar Wind Speed (Km/s) (stride = '+ str(stride) + ' Km/s)')
plt.ylabel('PDF')
plt.legend()

plt.subplot(2, 3, 2)
SWNHST, bins, stride = getDistrib(SWN[sEpochYYID:eEpochYYID], nbins=200, norm=True)
plt.plot(bins,SWNHST,label='Jan-Dec')
SWNHST = getDistrib(SWN[sEpochJJID:eEpochJJID], bins=bins, norm=True)
plt.plot(bins,SWNHST,label='Jan-Jun')
SWNHST = getDistrib(SWN[sEpochJDID:eEpochJDID], bins=bins, norm=True)
plt.plot(bins,SWNHST,label='Jul-Dec')
plt.title('Solar Wind Density PDF')
plt.xlabel('Solar Wind Density (N/cc) (stride = '+ str(stride) + ' N/cc)')
plt.ylabel('PDF')
plt.legend()

plt.subplot(2, 3, 3)
SWTHST, bins, stride = getDistrib(SWT[sEpochYYID:eEpochYYID], nbins=200, norm=True)
plt.plot(bins,SWTHST,label='Jan-Dec')
SWTHST = getDistrib(SWT[sEpochJJID:eEpochJJID], nbins=200, norm=True)
plt.plot(bins,SWTHST,label='Jan-Jun')
SWTHST = getDistrib(SWT[sEpochJDID:eEpochJDID], nbins=200, norm=True)
plt.plot(bins,SWTHST,label='Jul-Dec')
plt.title('Solar Wind Temperature PDF')
plt.xlabel('Solar Wind Temperature ($K^o$) (stride = '+ str(stride) + ' $K^o$)')
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

