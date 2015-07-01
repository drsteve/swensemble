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
from getswdata import getOMNIfiles, getACEfiles, getWINDfiles, getIMP8files
from getswdata import getOMNIdata
from getswdata import dataClean, dateShift, dateList

dataSRC = 'omni'
dataRecTime = 'hourly'

startDate   = datetime.datetime(1999, 1, 1, 0, 0, 0)
endDate     = datetime.datetime(1999,12,31,23,59,59)
locDateList = dateList(startDate, endDate, shift = 'day')

swParms = ['V','N','T','BX_GSE','BY_GSE','BZ_GSE','ABS_B']
omniData  = getOMNIdata(locDateList,'/home/ehab/SWData',swParms,dataRecTime,dataStat='clean')

Sp, Tr, VA, SWClass = getSolarWindType(omniData)

fig = plt.figure(1)
plt.suptitle('PDE Speed and Density in Solar Wind Categories for ' + str(startDate.year))

SWVEJT = [omniData['V'][item] for item in range(len(SWClass)) if SWClass[item] == 'EJT']
SWNEJT = [omniData['N'][item] for item in range(len(SWClass)) if SWClass[item] == 'EJT']
plt.subplot(2, 4, 1)
plt.title('Ejecta')
SWVHST, bins = getDistrib(filter(lambda v: v==v, SWVEJT), stride = 5, norm=True)
plt.plot(bins,SWVHST,label='Ejecta')
plt.xlabel('Speed (Km/s)')
plt.ylabel('Solar Wind Speed PDE')
plt.subplot(2, 4, 5)
SWNHST, bins = getDistrib(filter(lambda v: v==v, SWNEJT), stride = 1, norm=True)
plt.plot(bins,SWNHST,label='Ejecta')
plt.xlabel('Density (N/cc)')
plt.ylabel('Solar Wind Density PDE')

SWVEJT = [omniData['V'][item] for item in range(len(SWClass)) if SWClass[item] == 'CHO']
SWNEJT = [omniData['N'][item] for item in range(len(SWClass)) if SWClass[item] == 'CHO']
fig = plt.figure(1)
plt.subplot(2, 4, 2)
plt.title('Coronal Hole Origin')
SWVHST, bins = getDistrib(filter(lambda v: v==v, SWVEJT), stride = 5, norm=True)
plt.plot(bins,SWVHST,label='Coronal-Hole')
plt.xlabel('Speed (Km/s)')
plt.subplot(2, 4, 6)
SWNHST, bins = getDistrib(filter(lambda v: v==v, SWNEJT), stride = 1, norm=True)
plt.plot(bins,SWNHST,label='Coronal-Hole')
plt.xlabel('Density (N/cc)')

SWVEJT = [omniData['V'][item] for item in range(len(SWClass)) if SWClass[item] == 'SRR']
SWNEJT = [omniData['N'][item] for item in range(len(SWClass)) if SWClass[item] == 'SRR']
fig = plt.figure(1)
plt.subplot(2, 4, 3)
plt.title('Sector Reverse Origin')
SWVHST, bins = getDistrib(filter(lambda v: v==v, SWVEJT), stride = 5, norm=True)
plt.plot(bins,SWVHST,label='Sector-Reverse')
plt.xlabel('Speed (Km/s)')
plt.subplot(2, 4, 7)
SWNHST, bins = getDistrib(filter(lambda v: v==v, SWNEJT), stride = 1, norm=True)
plt.plot(bins,SWNHST,label='Sector-Reverse')
plt.xlabel('Density (N/cc)')

SWVEJT = [omniData['V'][item] for item in range(len(SWClass)) if SWClass[item] == 'SBO']
SWNEJT = [omniData['N'][item] for item in range(len(SWClass)) if SWClass[item] == 'SBO']
fig = plt.figure(1)
plt.subplot(2, 4, 4)
plt.title('Streamer Belt Origin')
SWVHST, bins = getDistrib(filter(lambda v: v==v, SWVEJT), stride = 5, norm=True)
plt.plot(bins,SWVHST,label='Streamer-Belt')
plt.xlabel('Speed (Km/s)')
plt.subplot(2, 4, 8)
SWNHST, bins = getDistrib(filter(lambda v: v==v, SWNEJT), stride = 1, norm=True)
plt.plot(bins,SWNHST,label='Streamer-Belt')
plt.xlabel('Density (N/cc)')


fig = plt.figure(2)
plt.suptitle('Solar Wind Categories for ' + str(startDate.year))
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


