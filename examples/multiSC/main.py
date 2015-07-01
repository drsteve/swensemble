import os
import sys
sys.path.append('/home/ehab/MyFiles/Softex/spacePy/spacepy-0.1.5')
import math
import numpy
import scipy
import bisect
import scipy.stats
import datetime
from numpy import transpose
from numpy.linalg import norm
from datetime import date
from spacepy import pycdf
from spacepy import seapy
from itertools import chain
import spacepy.time as spt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import simps, trapz
from swdatanal import getDistrib, getIndices, omniDataCorr, getSolarWindType
from getswdata import getOMNIfiles, getOMNIdata, getOMNIparams
from getswdata import getACEfiles, getACEdata, getACEparams
from getswdata import getIMP8files,getIMP8data, getIMP8params
from getswdata import getWINDfiles,getWINDdata, getWINDparams
from getswdata import dataClean, dateShift, dateList


startDate   = datetime.datetime(2003,3, 1, 0, 0, 0)
endDate     = datetime.datetime(2003,3,31,23,59,59)
locDateList = dateList(startDate, endDate, shift = 'day')

dataSRC = ""

dataSRC = 'wind'
if dataSRC == 'wind':
 dataRecTime = '1min'
 windParams = getWINDparams(locDateList,'/home/ehab/SWData',dataRecTime,windDat='all')
 if windParams != '--empty--':
  swParams = ['BZ','Proton_Np_nonlin','Proton_V_nonlin']
  windData = getWINDdata(locDateList,'/home/ehab/SWData',swParams,dataRecTime,dataStat='clean')

#dataSRC = 'imp8'
if dataSRC == 'imp8':
 dataRecTime = '4min'
 imp8Params = getIMP8params(locDateList,'/home/ehab/SWData',dataRecTime,imp8Dat='plasma')
 dataRecTime = '15sec'
 imp8Params = getIMP8params(locDateList,'/home/ehab/SWData',dataRecTime,imp8Dat='magnetic')
 if imp8Params != '--empty--':
  swParams = ['B_Vector_GSE','proton_density_fit','V_fit','protonV_thermal_fit']
  imp8Data = getIMP8data(locDateList,'/home/ehab/SWData',swParams,dataRecTime,dataStat='clean')


#dataSRC = 'ace'
if dataSRC == 'ace':
 dataRecTime = 'hourly'
 aceParams = getACEparams(locDateList,'/home/ehab/SWData',dataRecTime,aceDat='magnetic')
 aceParams = getACEparams(locDateList,'/home/ehab/SWData',dataRecTime,aceDat='plasma')
 if aceParams != '--empty--':
  swParams = ['Magnitude','Np','Vp','Tpr','BGSEc','V_GSE']
  aceData  = getACEdata(locDateList,'/home/ehab/SWData',swParams,dataRecTime,dataStat='clean')
  swCatParams = {'V':aceData['Vp'], 'T':aceData['Tpr'], 'N':aceData['Np'], 'B':aceData['Magnitude']}

  _, _, _, SWClass = getSolarWindType(swCatParams)
  aceVEJT=[];aceVCHO=[];aceVSRR=[];aceVSBO=[]
  aceNEJT=[];aceNCHO=[];aceNSRR=[];aceNSBO=[]
  aceBEJT=[];aceBCHO=[];aceBSRR=[];aceBSBO=[]
  acePEEJT=[];acePECHO=[];acePESRR=[];acePESBO=[]
  aceBEEJT=[];aceBECHO=[];aceBESRR=[];aceBESBO=[]
  for item in range(len(SWClass)):
   if SWClass[item] == 'EJT':
    aceVEJT.extend([aceData['Vp'][item]])
    aceNEJT.extend([aceData['Np'][item]])
    acePEEJT.extend([aceData['plasmaEpoch'][item]])
    aceBEEJT.extend([aceData['magneticEpoch'][item]])
   elif SWClass[item] == 'CHO':
    aceVCHO.extend([aceData['Vp'][item]])
    aceNCHO.extend([aceData['Np'][item]])
    aceBCHO.extend([aceData['Magnitude'][item]])
    acePECHO.extend([aceData['plasmaEpoch'][item]])
    aceBECHO.extend([aceData['magneticEpoch'][item]])
   elif SWClass[item] == 'SRR':
    aceVSRR.extend([aceData['Vp'][item]])
    aceNSRR.extend([aceData['Np'][item]])
    aceBSRR.extend([aceData['Magnitude'][item]])
    acePESRR.extend([aceData['plasmaEpoch'][item]])
    aceBESRR.extend([aceData['magneticEpoch'][item]])
   elif SWClass[item] == 'SBO':
    aceVSBO.extend([aceData['Vp'][item]])
    aceNSBO.extend([aceData['Np'][item]])
    aceBSBO.extend([aceData['Magnitude'][item]])
    acePESBO.extend([aceData['plasmaEpoch'][item]])
    aceBESBO.extend([aceData['magneticEpoch'][item]])

  plt.figure(1)
  plt.subplot(3,1,1)
  plt.plot(acePEEJT,aceVEJT, color = 'b', marker = '*', label = 'Ejecta')
  plt.plot(acePECHO,aceVCHO, color = 'r', marker = '*', label = 'Conronal-Hole')
  plt.plot(acePESRR,aceVSRR, color = 'm', marker = '*', label = 'Sector-Reversal')
  plt.plot(acePESBO,aceVSBO, color = 'g', marker = '*', label = 'Streamer-Belt')
  plt.ylabel('Flow Velocity (km/s)')
  plt.title('ACE Data Solar Wind Categorization March 2003')
  plt.legend(loc='best')
  plt.subplot(3,1,2)
  plt.plot(acePEEJT,aceNEJT, color = 'b', marker = '*', label = 'Ejecta')
  plt.plot(acePECHO,aceNCHO, color = 'r', marker = '*', label = 'Conronal-Hole')
  plt.plot(acePESRR,aceNSRR, color = 'm', marker = '*', label = 'Sector-Reversal')
  plt.plot(acePESBO,aceNSBO, color = 'g', marker = '*', label = 'Streamer-Belt')
  plt.ylabel('Flow Density (N/cc)')
  plt.legend(loc='best')
  plt.subplot(3,1,3)
 #plt.plot(aceBEEJT,aceBEJT, color = 'b', marker = '*', label = 'Ejecta')
  plt.plot(aceBECHO,aceBCHO, color = 'r', marker = '*', label = 'Conronal-Hole')
  plt.plot(aceBESRR,aceBSRR, color = 'm', marker = '*', label = 'Sector-Reversal')
  plt.plot(aceBESBO,aceBSBO, color = 'g', marker = '*', label = 'Streamer-Belt')
  plt.ylabel('$B_z$ (nT)')
  plt.xlabel('Date')
  plt.legend(loc='best')
  plt.show()


#dataSRC = 'omni'
if dataSRC == 'omni':
 dataRecTime = 'hourly'
 omniParams = getOMNIparams(locDateList,'/home/ehab/SWData',dataRecTime)
 if omniParams != '--empty--':
  swParams = ['V','N','T','BX_GSE','BY_GSE','BZ_GSE','ABS_B']
  omniData  = getOMNIdata(locDateList,'/home/ehab/SWData',swParams,dataRecTime,dataStat='clean')
  swCatParams = {'V':omniData['V'], 'T':omniData['T'], 'N':omniData['N'], 'B':omniData['ABS_B']}

  _, _, _, SWClass = getSolarWindType(swCatParams)
  omniVEJT=[];omniVCHO=[];omniVSRR=[];omniVSBO=[]
  omniNEJT=[];omniNCHO=[];omniNSRR=[];omniNSBO=[]
  omniBEJT=[];omniBCHO=[];omniBSRR=[];omniBSBO=[]
  omniEEJT=[];omniECHO=[];omniESRR=[];omniESBO=[]
  for item in range(len(SWClass)):
   if SWClass[item] == 'EJT':
    omniVEJT.extend([omniData['V'][item]])
    omniNEJT.extend([omniData['N'][item]])
    omniBEJT.extend([omniData['ABS_B'][item]])
    omniEEJT.extend([omniData['Epoch'][item]])
   elif SWClass[item] == 'CHO':
    omniVCHO.extend([omniData['V'][item]])
    omniNCHO.extend([omniData['N'][item]])
    omniBCHO.extend([omniData['ABS_B'][item]])
    omniECHO.extend([omniData['Epoch'][item]])
   elif SWClass[item] == 'SRR':
    omniVSRR.extend([omniData['V'][item]])
    omniNSRR.extend([omniData['N'][item]])
    omniBSRR.extend([omniData['ABS_B'][item]])
    omniESRR.extend([omniData['Epoch'][item]])
   elif SWClass[item] == 'SBO':
    omniVSBO.extend([omniData['V'][item]])
    omniNSBO.extend([omniData['N'][item]])
    omniBSBO.extend([omniData['ABS_B'][item]])
    omniESBO.extend([omniData['Epoch'][item]])

  plt.figure(1)
  plt.subplot(3,1,1)
  plt.plot(omniEEJT,omniVEJT, color = 'b', marker = 'o', label = 'Ejecta')
  plt.plot(omniECHO,omniVCHO, color = 'r', marker = 'o', label = 'Conronal-Hole')
  plt.plot(omniESRR,omniVSRR, color = 'm', marker = 'o', label = 'Sector-Reversal')
  plt.plot(omniESBO,omniVSBO, color = 'g', marker = 'o', label = 'Streamer-Belt')
  plt.ylabel('Flow Velocity (km/s)')
  plt.legend(loc='best')
  plt.subplot(3,1,2)
  plt.plot(omniEEJT,omniNEJT, color = 'b', marker = 'o', label = 'Ejecta')
  plt.plot(omniECHO,omniNCHO, color = 'r', marker = 'o', label = 'Conronal-Hole')
  plt.plot(omniESRR,omniNSRR, color = 'm', marker = 'o', label = 'Sector-Reversal')
  plt.plot(omniESBO,omniNSBO, color = 'g', marker = 'o', label = 'Streamer-Belt')
  plt.ylabel('Flow Density (N/cc)')
  plt.legend(loc='best')
  plt.subplot(3,1,3)
  plt.plot(omniEEJT,omniBEJT, color = 'b', marker = 'o', label = 'Ejecta')
  plt.plot(omniECHO,omniBCHO, color = 'r', marker = 'o', label = 'Conronal-Hole')
  plt.plot(omniESRR,omniBSRR, color = 'm', marker = 'o', label = 'Sector-Reversal')
  plt.plot(omniESBO,omniBSBO, color = 'g', marker = 'o', label = 'Streamer-Belt')
  plt.ylabel('$B_z$ (nT)')
  plt.xlabel('Date')
  plt.legend(loc='best')
  plt.title('OMNI Data Solar Wind Categorization')
  plt.show()

sys.exit()

plt.figure(2)
plt.subplot(3,1,1)
plt.plot(windData['Epoch'],windData['Proton_Np_nonlin'],label='WIND')
plt.plot(imp8Data['plasmaEpoch'],imp8Data['proton_density_fit'],label='IMP8')
plt.plot(aceData['plasmaEpoch'],aceData['Np'],label='ACE')
plt.plot(omniData['Epoch'],omniData['N'],label='OMNI')
plt.title('Proton Number Density')
plt.ylabel('Proton Number Density (N/cc)')
plt.xlim(startDate,endDate)
plt.legend(loc='best')

plt.subplot(3,1,2)
plt.plot(windData['Epoch'],windData['Proton_V_nonlin'],label='WIND')
plt.plot(imp8Data['plasmaEpoch'],imp8Data['V_fit'],label='IMP8')
plt.plot(aceData['plasmaEpoch'],aceData['Vp'],label='ACE')
plt.plot(omniData['Epoch'],omniData['V'],label='OMNI')
plt.title('Proton Flow Velocity')
plt.ylabel('Proton Flow Velocity (Km/s)')
plt.xlim(startDate,endDate)
plt.legend(loc='best')

plt.subplot(3,1,3)
plt.plot(windData['Epoch'],windData['BZ'],label='WIND')
plt.plot(imp8Data['magneticEpoch'],imp8Data['B_Vector_GSE'][2],label='IMP8')
plt.plot(aceData['magneticEpoch'],aceData['BGSEc'][2],label='ACE')
plt.plot(omniData['Epoch'],omniData['BZ_GSE'],label='OMNI')
plt.title('$B_z$')
plt.xlabel('Date')
plt.ylabel('$B_z$ (nT)')
plt.xlim(startDate,endDate)
plt.legend(loc='best')

plt.suptitle('Measurements from Different Spacecrafts')
plt.show()

sys.exit()

swCatParams = {'V':omniData['V'], 'T':omniData['T'], 'N':omniData['N'], 'B':omniData['ABS_B']}
Sp, Tr, VA, SWClass = getSolarWindType(swCatParams)

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

if dataSRC == 'omni':
 OMNIfnames  = getOMNIfiles(locDateList,'/home/ehab/SWData',dataRecTime)
 print OMNIfnames
elif dataSRC == 'ace':
 ACEfnames  = getACEfiles(locDateList,'/home/ehab/SWData',dataRecTime,'plasma')
 print ACEfnames
elif dataSRC == 'wind':
 WINDfnames  = getWINDfiles(locDateList,'/home/ehab/SWData',dataRecTime,'iplasma')
 print WINDfnames
elif dataSRC == 'imp8':
 WINDfnames  = getIMP8files(locDateList,'/home/ehab/SWData',dataRecTime,'plasma')
 print WINDfnames

cdfData = []
for i in range(len(OMNIfnames)):
 cdf = pycdf.CDF(OMNIfnames[i])
 cdfData.append(cdf.copy())
cdf.close()

#print cdfData[0].keys()

sys.exit()

epochs = []; SWV = []; SWN = []; SWT = []; SWB = []
for i in range(len(OMNIfnames)): 
 epochs = epochs + list(cdfData[i]['Epoch'])
 if dataSRC == 'omni':
  if datRecTime == 'hourly':
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
 elif dataSRC == 'ace':
  print 'coming soon'

SWN = dataClean(SWN,[999.0,0.0],['>=','<='])
SWV = dataClean(SWV,[2500],['>='])
SWT = dataClean(SWT,[1.0e7],['>='])
SWB = dataClean(SWB,[0.0],['<='])

SWPList = {'V': SWV, 'N': SWN, 'T': SWT, 'B': SWB}

Sp, Tr, VA, SWClass = getSolarWindType(SWPList)

fig = plt.figure(1)
plt.suptitle('PDE Speed and Density in Solar Wind Categories for ' + str(startDate.year))

SWVEJT = [SWV[item] for item in range(len(SWClass)) if SWClass[item] == 'EJT']
SWNEJT = [SWN[item] for item in range(len(SWClass)) if SWClass[item] == 'EJT']
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

SWVEJT = [SWV[item] for item in range(len(SWClass)) if SWClass[item] == 'CHO']
SWNEJT = [SWN[item] for item in range(len(SWClass)) if SWClass[item] == 'CHO']
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

SWVEJT = [SWV[item] for item in range(len(SWClass)) if SWClass[item] == 'SRR']
SWNEJT = [SWN[item] for item in range(len(SWClass)) if SWClass[item] == 'SRR']
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

SWVEJT = [SWV[item] for item in range(len(SWClass)) if SWClass[item] == 'SBO']
SWNEJT = [SWN[item] for item in range(len(SWClass)) if SWClass[item] == 'SBO']
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

