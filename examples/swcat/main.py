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
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import simps, trapz
from swdatanal import getDistrib, omniDataCorr, getSolarWindType
from getswdata import getOMNIfiles, dataClean, dateShift, dateList

#omniRecTime = 'hourly'
omniRecTime = '5min'

startDate   = datetime.datetime(2014, 1, 1, 0, 0, 0)
endDate     = datetime.datetime(2014, 2,28,23,59,59)
locDateList = dateList(startDate, endDate, shift = 'month')
OMNIfnames  = getOMNIfiles(locDateList,'/home/ehab/SWData',omniRecTime)

cdfData = []
for i in range(len(OMNIfnames)):
 cdf = pycdf.CDF(OMNIfnames[i])
 cdfData.append(cdf.copy())
cdf.close()

epochs = []; SWV = []; SWN = []; SWT = []; SWB = []
for i in range(len(OMNIfnames)): 
 epochs = epochs + list(cdfData[i]['Epoch'])
 if omniRecTime == 'hourly':
  SWN    = SWN + list(cdfData[i]['N'])
  SWT    = SWT + list(cdfData[i]['T'])
  SWV    = SWV + list(cdfData[i]['V'])
  SWB    = SWB + list(cdfData[i]['ABS_B'])
 else:
  SWN    = SWN + list(cdfData[i]['proton_density'])
  SWT    = SWT + list(cdfData[i]['T'])
  SWV    = SWV + list(cdfData[i]['flow_speed'])
  for j in range(len(cdfData[i]['BX_GSE'])):
   B     = numpy.linalg.norm([cdfData[i]['BX_GSE'][j],cdfData[i]['BX_GSE'][j],cdfData[i]['BX_GSE'][j]])
   SWB   = SWB + [B]

SWN = dataClean(SWN,[999.0,0.0],['>=','<='])
SWV = dataClean(SWV,[2500],['>='])
SWT = dataClean(SWT,[1.0e7],['>='])
SWB = dataClean(SWB,[0.0],['<='])

SWPList = {'V': SWV, 'N': SWN, 'T': SWT, 'B': SWB}

Sp, Tr, VA, SWClass = getSolarWindType(SWPList)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

EJTSp = [Sp[item] for item in range(len(SWClass)) if SWClass[item] == 'EJT']
EJTVA = [VA[item] for item in range(len(SWClass)) if SWClass[item] == 'EJT']
EJTTr = [Tr[item] for item in range(len(SWClass)) if SWClass[item] == 'EJT']
if EJTSp != []:
 ax.scatter(EJTSp, EJTVA, EJTTr, color = 'b', marker = 'o', label = 'Ejecta')

CHOSp = [Sp[item] for item in range(len(SWClass)) if SWClass[item] == 'CHO']
CHOVA = [VA[item] for item in range(len(SWClass)) if SWClass[item] == 'CHO']
CHOTr = [Tr[item] for item in range(len(SWClass)) if SWClass[item] == 'CHO']
if CHOSp != []:
 ax.scatter(CHOSp, CHOVA, CHOTr, color = 'r', marker = 'o', label = 'Coronal-Hole')

SRRSp = [Sp[item] for item in range(len(SWClass)) if SWClass[item] == 'SRR']
SRRVA = [VA[item] for item in range(len(SWClass)) if SWClass[item] == 'SRR']
SRRTr = [Tr[item] for item in range(len(SWClass)) if SWClass[item] == 'SRR']
if SRRSp != []:
 ax.scatter(SRRSp, SRRVA, SRRTr, color = 'm', marker = 'o', label = 'Sector-Reversal')

SBOSp = [Sp[item] for item in range(len(SWClass)) if SWClass[item] == 'SBO']
SBOVA = [VA[item] for item in range(len(SWClass)) if SWClass[item] == 'SBO']
SBOTr = [Tr[item] for item in range(len(SWClass)) if SWClass[item] == 'SBO']
if SBOSp != []:
 ax.scatter(SBOSp, SBOVA, SBOTr, color = 'g', marker = 'o', label = 'Streamer-Belt')

plt.legend()
ax.set_xlabel('$S_p = \\frac{T_p}{n_p^{2/3}}$ (eV.$cm^2$)')
ax.set_ylabel('$V_A = \\frac{B}{(4\pi m_pn_p)^{1/2}}$ (Km/s)')
ax.set_zlabel('$T_r=\\frac{T_{exp}}{T_p}$')
plt.show()

sys.exit()

sEpoch = datetime.datetime(2014, 1, 1, 0, 0)
eEpoch = datetime.datetime(2014,12,31,23,59)
sEpochID = bisect.bisect_left(epochs, sEpoch)
eEpochID = bisect.bisect_left(epochs, eEpoch)

SWV = filter(lambda v: v==v, SWV)
SWN = filter(lambda v: v==v, SWN)
SWT = filter(lambda v: v==v, SWT)
SWB = filter(lambda v: v==v, SWB)

