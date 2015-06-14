import os
import sys
sys.path.append('/home/ehab/MyFiles/Softex/spacePy/spacepy-0.1.5')
import math
import numpy
import scipy
import scipy.stats
import datetime
from spacepy import pycdf
from itertools import chain
import spacepy.time as spt
import matplotlib.pyplot as plt
from getswdata import getOMNIfiles

OMNIfnames = getOMNIfiles(2012,'/home/ehab/SWData','hourly')
print OMNIfnames[0]
print OMNIfnames[1]
sys.exit()

cdf = pycdf.CDF(OMNIfnames[0])
cdfData01 = cdf.copy()
cdf = pycdf.CDF(OMNIfnames[1])
cdfData02 = cdf.copy()
cdf.close()

def dataClean(listName,threshold):
    for i in listName:
     if i >= threshold: listName.remove(i)
    return listName

#print cdfData01.keys()

plt.subplot(2, 3, 1)
SWV = list(chain(*zip(cdfData01['V'],cdfData02['V'])))
SWV = dataClean(SWV,2500)
stride = (max(SWV)-min(SWV))/100
bins=numpy.arange(min(SWV)-stride,max(SWV)+stride,stride)
SWVHST = scipy.stats.histogram2(SWV,bins)
plt.plot(bins,SWVHST,label='Jan-Dec')
SWV = cdfData01['V']
SWVHST = scipy.stats.histogram2(SWV,bins)
plt.plot(bins,SWVHST,label='Jan-Jun')
SWV = cdfData02['V']
SWVHST = scipy.stats.histogram2(SWV,bins)
plt.plot(bins,SWVHST,label='Jul-Dec')
plt.title('Solar Wind Speed PDF')
plt.xlabel('Solar Wind Speed (Km/s) (stride = '+ str(stride) + ' Km/s)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 2)
SWN = list(chain(*zip(cdfData01['N'],cdfData02['N'])))
SWN = dataClean(SWN,999)
stride = (max(SWN)-min(SWN))/100
bins=numpy.arange(min(SWN)-stride,max(SWN)+stride,stride)
SWNHST = scipy.stats.histogram2(SWN,bins)
plt.plot(bins,SWNHST,label='Jan-Dec')
SWN = cdfData01['N']
SWNHST = scipy.stats.histogram2(SWN,bins)
plt.plot(bins,SWNHST,label='Jan-Jun')
SWN = cdfData02['N']
SWNHST = scipy.stats.histogram2(SWN,bins)
plt.plot(bins,SWNHST,label='Jul-Dec')
plt.title('Solar Wind Density PDF')
plt.xlabel('Solar Wind Density (N/cc) (stride = '+ str(stride) + ' N/cc)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 3)
SWT = list(chain(*zip(cdfData01['T'],cdfData02['T'])))
stride = (max(SWT)-min(SWT))/100
bins=numpy.arange(min(SWT)-stride,max(SWT)+stride,stride)
SWTHST = scipy.stats.histogram2(SWT,bins)
plt.plot(bins,SWTHST,label='Jan-Dec')
SWT = cdfData01['T']
SWTHST = scipy.stats.histogram2(SWT,bins)
plt.plot(bins,SWTHST,label='Jan-Jun')
SWT = cdfData02['T']
SWTHST = scipy.stats.histogram2(SWT,bins)
plt.plot(bins,SWTHST,label='Jul-Dec')
plt.title('Solar Wind Temperature PDF')
plt.xlabel('Solar Wind Temperature ($K^o$) (stride = '+ str(stride) + ' $K^o$)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 4)
SWB = list(chain(*zip(cdfData01['BX_GSE'],cdfData02['BX_GSE'])))
stride = (max(SWB)-min(SWB))/100
bins=numpy.arange(min(SWB)-stride,max(SWB)+stride,stride)
SWBHST = scipy.stats.histogram2(SWB,bins)
plt.plot(bins,SWBHST,label='Jan-Dec')
SWB = cdfData01['BX_GSE']
SWBHST = scipy.stats.histogram2(SWB,bins)
plt.plot(bins,SWBHST,label='Jan-Jun')
SWB = cdfData02['BX_GSE']
SWBHST = scipy.stats.histogram2(SWB,bins)
plt.plot(bins,SWBHST,label='Jul-Dec')
plt.title('Solar Wind Bx PDF')
plt.xlabel('Solar Wind Bx (nT) (stride = '+ str(stride) + ' nT)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 5)
SWB = list(chain(*zip(cdfData01['BY_GSE'],cdfData02['BY_GSE'])))
stride = (max(SWB)-min(SWB))/100
bins=numpy.arange(min(SWB)-stride,max(SWB)+stride,stride)
SWBHST = scipy.stats.histogram2(SWB,bins)
plt.plot(bins,SWBHST,label='Jan-Dec')
SWB = cdfData01['BY_GSE']
SWBHST = scipy.stats.histogram2(SWB,bins)
plt.plot(bins,SWBHST,label='Jan-Jun')
SWB = cdfData02['BY_GSE']
SWBHST = scipy.stats.histogram2(SWB,bins)
plt.plot(bins,SWBHST,label='Jul-Dec')
plt.title('Solar Wind By PDF')
plt.xlabel('Solar Wind By (nT) (stride = '+ str(stride) + ' nT)')
plt.ylabel('Number of Occurance')
plt.legend()

plt.subplot(2, 3, 6)
SWB = list(chain(*zip(cdfData01['BZ_GSE'],cdfData02['BZ_GSE'])))
stride = (max(SWB)-min(SWB))/100
bins=numpy.arange(min(SWB)-stride,max(SWB)+stride,stride)
SWBHST = scipy.stats.histogram2(SWB,bins)
plt.plot(bins,SWBHST,label='Jan-Dec')
SWB = cdfData01['BZ_GSE']
SWBHST = scipy.stats.histogram2(SWB,bins)
plt.plot(bins,SWBHST,label='Jan-Jun')
SWB = cdfData02['BZ_GSE']
SWBHST = scipy.stats.histogram2(SWB,bins)
plt.plot(bins,SWBHST,label='Jul-Dec')
plt.title('Solar Wind Bz PDF')
plt.xlabel('Solar Wind Bz (nT) (stride = '+ str(stride) + ' nT)')
plt.ylabel('Number of Occurance')
plt.legend()


plt.suptitle('OMNI Hourly Data for 2014')
plt.show()



sys.exit()





#HUR = cdfData['HR']
#cdfDates=[0 for i in range((len(cdfData['YR'])))]
#cdfDates=list(range(len(cdfData['YR'])))
#for i in range(len(cdfData['YR'])):
# d = datetime.date(int(cdfData['YR'][i]),1,1) + datetime.timedelta(days=int(cdfData['Day'][i]-1))
# t = datetime.time(int(cdfData['HR'][i]))
# cdfDates[i] = datetime.combine(d,t)
#print datetime.combine(d,t), type(datetime.combine(d,t))
#print cdfDates[i]
#plt.plot(HUR,SWV)
#plt.show()
