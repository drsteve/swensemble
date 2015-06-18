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

startDate   = datetime.datetime(2014,3,1,0,0,0)
endDate     = datetime.datetime(2014,5,31,23,59,59)
locDateList = dateList(startDate, endDate, shift = 'month')
srefDate    = datetime.datetime(2008,1,1,0,0,0)
erefDate    = datetime.datetime(2008,1,1,23,59,59)
refDateList = dateList(srefDate, erefDate, shift = 'month')
OMNIfnames  = getOMNIfiles(locDateList,'/home/ehab/SWData','hourly')
if refDateList in locDateList:
 pass
else:
 OMNIfnames = sorted(OMNIfnames + getOMNIfiles(refDateList,'/home/ehab/SWData','hourly'))

cdfData = []
for i in range(len(OMNIfnames)):
 cdf = pycdf.CDF(OMNIfnames[i])
 cdfData.append(cdf.copy())
cdf.close()

epochs = []; SWP = []
for i in range(len(OMNIfnames)): 
 epochs = epochs + list(cdfData[i]['Epoch'])
#SWP    = SWP + list(cdfData[i]['proton_density'])
 SWP    = SWP + list(cdfData[i]['N'])
SWP = dataClean(SWP,[999.0,0.0],['>=','<'])
#SWP = dataClean(SWP,[2500],['>='])
stride = 10
SWPDateRng, cepochs, KSVals, KSDist, aepochs = omniDataCorr(srefDate, erefDate, startDate, endDate, epochs, SWP, binStride = stride, CorrTime = 'Day', CorrType = 'pearson')

if SWPDateRng > 0: 
 fig = plt.figure()
 ax = plt.subplot2grid((2,2), (0,0), colspan=2)
 plt.plot(cepochs,SWPDateRng)
 plt.title('Solar Wind Density')
 plt.ylabel('Density (N/cc)')

 ax = plt.subplot2grid((2,2),(1, 0))
 plt.plot(aepochs,KSVals[:,0],label='KS-Stat')
 plt.plot(aepochs,KSVals[:,1],label='P-Value')
 plt.title('Kolmogorov-Smirnov test between Density Samples')
 plt.xlabel('Month')
 plt.ylabel('KS-Stat and P-Value')
 plt.xticks(rotation='30')
 plt.legend(loc='best')

 ax = plt.subplot2grid((2,2),(1, 1))
 plt.plot(aepochs,KSDist[:,0],label='KS-Stat')
 plt.plot(aepochs,KSDist[:,1],label='P-Value')
 plt.title('Kolmogorov-Smirnov test between Density Distributions')
 plt.xlabel('Month')
 plt.ylabel('KS-Stat and P-Value')
 plt.xticks(rotation='30')
 plt.legend(loc='best')

 plt.suptitle('Compare the Density Daily Variation Relative to First Date (stride = ' + str(stride) + ')')
 plt.show()

sys.exit()

