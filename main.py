import os
import sys
sys.path.append('/home/ehab/MyFiles/Softex/spacePy/spacepy-0.1.5')
import math
import numpy
import scipy
import bisect
import scipy.stats
import datetime, time
from numpy import transpose
from numpy.linalg import norm
from datetime import date
from spacepy import pycdf
from spacepy import seapy
from itertools import chain
from scipy.stats import pearsonr, spearmanr
import spacepy.time as spt
import matplotlib.pyplot as plt
from scipy.signal import medfilt
from mpl_toolkits.mplot3d import Axes3D

from scipy.integrate import simps, trapz
from scipy.stats import pearsonr, gaussian_kde

from swdatanal import getDistrib, getIndices, omniDataCorr, ccorr, xcorr
from swdatanal import kdeBW, getDesKDE, epochShift, findContiguousData
from swdatanal import search, epochBlock, findCorrEpoch
from swdatanal import getSolarWindType, getTimeLag, swMedFilter
from getswdata import getOMNIfiles, getOMNIdata, getOMNIparams
from getswdata import getACEfiles, getACEdata, getACEparams, aceDataAdjust
from getswdata import getIMP8files,getIMP8data, getIMP8params, imp8DataAdjust
from getswdata import getWINDfiles,getWINDdata, getWINDparams
from getswdata import getGeotailfiles,getGeotaildata, getGeotailparams
from getswdata import dataClean, dateShift, dateList, epochMatch, removeNaN

#startDate   = datetime.datetime(1998, 1,14, 0, 0, 0)
#endDate     = datetime.datetime(1998, 1,16, 0, 0, 0)
#startDate   = datetime.datetime(2000, 1,14,18, 0, 0)
#endDate     = datetime.datetime(2000, 1,15,12, 0, 0)
#startDate   = datetime.datetime(1999, 7,23, 0, 0, 0)
#endDate     = datetime.datetime(1999, 7,24, 0, 0, 0)
startDate   = datetime.datetime(1998, 1, 1, 0, 0, 0)
endDate     = datetime.datetime(2000,12,31,23,59,59)
locDateList = dateList(startDate, endDate, shift = 'hour')

#timeShift = 4
#iCounter = 0
#while 1:
# sTime.extend([dateShift(startDate, hours=iCounter*timeShift)])
# eTime.extend([dateShift(startDate, hours=(iCounter+1)*timeShift)])
# if eTime[iCounter] >= endDate: break
# iCounter = iCounter + 1

dataSRC = ""
swCats = False

#dataSRC = 'geotail'
if dataSRC == 'geotail':
 geotailParams = getGeotailparams(locDateList,'/home/ehab/SWData',geotailSet='1min')
 if geotailParams != '--empty--':
  swParams = ['V','N','T','BX_GSE','BY_GSE','BZ_GSE','ABS_B','X','Y','Z','VX_GSE','VY_GSE','VZ_GSE']
  geotailData = getGeotaildata(locDateList,'/home/ehab/SWData',swParams,geotailSet='1min',dataStat='clean')
  geotailData['B']   = geotailData['ABS_B']
 #geotailData['N']   = geotailData['N']
 #geotailData['V']   = geotailData['V']
 #geotailData['T']   = geotailData['T']
  geotailData['Bx']  = geotailData['BX_GSE']
  geotailData['By']  = geotailData['BY_GSE']
  geotailData['Bz']  = geotailData['BZ_GSE']
  geotailData['Vx']  = geotailData['VX_GSE']
  geotailData['Vy']  = geotailData['VY_GSE']
  geotailData['Vz']  = geotailData['VZ_GSE']
  geotailData['SCxGSE'] = geotailData['X']
  geotailData['SCyGSE'] = geotailData['Y']
  geotailData['SCzGSE'] = geotailData['Z']

#dataSRC = 'wind'
if dataSRC == 'wind':
 windParams = getWINDparams(locDateList,'/home/ehab/SWData',windSet='1min',windDat='merged')
#windParams = getWINDparams(locDateList,'/home/ehab/SWData',windSet='1min',windDat='magnetic')
#windParams = getWINDparams(locDateList,'/home/ehab/SWData',windSet='1min',windDat='plasma')
 if windParams != '--empty--':
  RE = 6371.0
  swParams = ['BX','BY','BZ','Proton_Np_nonlin','Proton_V_nonlin','Proton_VX_nonlin','Proton_VY_nonlin','Proton_VZ_nonlin','xgse','ygse','zgse']
  windData = getWINDdata(locDateList,'/home/ehab/SWData',swParams,windSet='1min',windDat='merged',dataStat='clean')
 #windData['B']   = numpy.array(windData['Magnitude'])
  windData['N']   = numpy.array(windData['Proton_Np_nonlin'])
  windData['V']   = numpy.array(windData['Proton_V_nonlin'])
 #windData['T']   = numpy.array(windData['protonV_thermal_fit'])
  windData['Bx']  = numpy.array(windData['BX'])
  windData['By']  = numpy.array(windData['BY'])
  windData['Bz']  = numpy.array(windData['BZ'])
  windData['Vx']  = numpy.array(windData['Proton_VX_nonlin'])
  windData['Vy']  = numpy.array(windData['Proton_VY_nonlin'])
  windData['Vz']  = numpy.array(windData['Proton_VZ_nonlin'])
  windData['SCxGSE'] = numpy.array(windData['xgse'])*RE
  windData['SCyGSE'] = numpy.array(windData['ygse'])*RE
  windData['SCzGSE'] = numpy.array(windData['zgse'])*RE

dataSRC = 'imp8'
if dataSRC == 'imp8':
 imp8Params = getIMP8params(locDateList,'/home/ehab/SWData','15sec',imp8Dat='magnetic')
 imp8Params = getIMP8params(locDateList,'/home/ehab/SWData','1min',imp8Dat='plasma')
 if imp8Params != '--empty--':
  RE = 6371.0
  swParams = ['B_Vector_GSE','proton_density_fit','V_fit','protonV_thermal_fit','SC_Pos_GSE','SC_Pos_GSM']
  imp8Data = getIMP8data(locDateList,'/home/ehab/SWData',swParams,imp8Set='15sec',dataStat='clean')
  imp8Data = imp8DataAdjust(imp8Data)

dataSRC = 'ace'
if dataSRC == 'ace':
 aceParams = getACEparams(locDateList,'/home/ehab/SWData',aceSet='16sec',aceDat='magnetic')
 aceParams = getACEparams(locDateList,'/home/ehab/SWData',aceSet='1min',aceDat='plasma')
 if aceParams != '--empty--':
  swParams = ['Magnitude','Np','Vp','Tpr','BGSEc','V_GSE','SC_pos_GSE']
  aceData  = getACEdata(locDateList,'/home/ehab/SWData',swParams,['16sec','1min'],dataStat='clean')
  aceData  = aceDataAdjust(aceData)

  if swCats:
   swCatParams = {'V':aceData['V'], 'T':aceData['T'], 'N':aceData['N'], 'B':aceData['B']}
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

  if swCats:
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



RE = 6371.0
plt.figure(200)
plt.subplot(3,1,1)
plt.plot(imp8Data['magneticEpoch'],imp8Data['SCxGSE']/RE,label='IMP8')
plt.plot( aceData['magneticEpoch'], aceData['SCxGSE']/RE,label='ACE')
plt.ylabel('X-Position ($R_E$)')
plt.title('ACE and IMP8 Location')
plt.legend(loc='best')
plt.subplot(3,1,2)
plt.plot(imp8Data['magneticEpoch'],imp8Data['SCyGSE']/RE,label='Imp8')
plt.plot( aceData['magneticEpoch'], aceData['SCyGSE']/RE,label='ACE')
plt.ylabel('Y-Position ($R_E$)')
plt.legend(loc='best')
plt.subplot(3,1,3)
plt.plot(imp8Data['magneticEpoch'],imp8Data['SCzGSE']/RE,label='IMP8')
plt.plot( aceData['magneticEpoch'], aceData['SCzGSE']/RE,label='ACE')
plt.ylabel('Z-Position ($R_E$)')
plt.xlabel('Epoch Time')
plt.legend(loc='best')
plt.show()
sys.exit()

timeShift = 4
timeGap = 900

pBlockStart = epochBlock(imp8Data['plasmaEpoch'], imp8Data['V'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(pBlockStart,pBlockStart)
pBlockStart = epochBlock(imp8Data['plasmaEpoch'], imp8Data['N'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(imp8Data['magneticEpoch'], imp8Data['Bz'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)

pBlockStart = epochBlock(aceData['plasmaEpoch'], aceData['Vx'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(aceData['plasmaEpoch'], aceData['V'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(aceData['plasmaEpoch'], aceData['N'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(aceData['magneticEpoch'], aceData['Bz'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)

sTime = []; eTime = []
for ind in range(len(refBlockStart)):
 sTime.extend([refBlockStart[ind]])
 eTime.extend([dateShift(refBlockStart[ind], hours = timeShift)])

accSrcParam = []
accDesParam = []

for i in range(len(sTime)):
 print sTime[i]
 uniDateList = dateList(sTime[i], eTime[i], shift = 'minute')

 sEpochBIDI  = bisect.bisect_left(imp8Data['magneticEpoch'], sTime[i])
 eEpochBIDI  = bisect.bisect_left(imp8Data['magneticEpoch'], eTime[i])
 sEpochPIDI  = bisect.bisect_left(imp8Data['plasmaEpoch'], sTime[i])
 eEpochPIDI  = bisect.bisect_left(imp8Data['plasmaEpoch'], eTime[i])
 impPEpoch   = numpy.array(imp8Data['plasmaEpoch'][sEpochPIDI:eEpochPIDI])
 impBEpoch   = numpy.array(imp8Data['magneticEpoch'][sEpochBIDI:eEpochBIDI])
 
 sEpochBIDA  = bisect.bisect_left(aceData['magneticEpoch'], sTime[i])
 eEpochBIDA  = bisect.bisect_left(aceData['magneticEpoch'], eTime[i])
 sEpochPIDA  = bisect.bisect_left(aceData['plasmaEpoch'], sTime[i])
 eEpochPIDA  = bisect.bisect_left(aceData['plasmaEpoch'], eTime[i])
 acePEpoch   = numpy.array(aceData['plasmaEpoch'][sEpochPIDA:eEpochPIDA])
 aceBEpoch   = numpy.array(aceData['magneticEpoch'][sEpochBIDA:eEpochBIDA])

 TT,VV   = removeNaN(acePEpoch,aceData['V'][sEpochPIDA:eEpochPIDA])
 aceVmod = epochMatch(uniDateList,TT,VV,interpKind='linear')
 TT,VV   = removeNaN(impPEpoch,imp8Data['V'][sEpochPIDI:eEpochPIDI])
 impVmod = epochMatch(uniDateList,TT,VV,interpKind='linear')

 TT,NN   = removeNaN(acePEpoch,aceData['N'][sEpochPIDA:eEpochPIDA])
 aceNmod = epochMatch(uniDateList,TT,NN,interpKind='linear')
 TT,NN   = removeNaN(impPEpoch,imp8Data['N'][sEpochPIDI:eEpochPIDI])
 impNmod = epochMatch(uniDateList,TT,NN,interpKind='linear')

#TT,BB   = removeNaN(aceBEpoch,aceData['Bx'][sEpochBIDA:eEpochBIDA])
#aceBxmod= epochMatch(uniDateList,TT,BB,interpKind='linear')
#TT,BB   = removeNaN(impBEpoch,imp8Data['Bx'][sEpochBIDI:eEpochBIDI])
#impBxmod= epochMatch(uniDateList,TT,BB,interpKind='linear')

#TT,BB   = removeNaN(aceBEpoch,aceData['By'][sEpochBIDA:eEpochBIDA])
#aceBymod= epochMatch(uniDateList,TT,BB,interpKind='linear')
#TT,BB   = removeNaN(impPEpoch,imp8Data['By'][sEpochBIDI:eEpochBIDI])
#impBymod= epochMatch(uniDateList,TT,BB,interpKind='linear')

 TT,BB   = removeNaN(aceBEpoch,aceData['Bz'][sEpochBIDA:eEpochBIDA])
 aceBzmod= epochMatch(uniDateList,TT,BB,interpKind='linear')
 TT,BB   = removeNaN(impBEpoch,imp8Data['Bz'][sEpochBIDI:eEpochBIDI])
 impBzmod= epochMatch(uniDateList,TT,BB,interpKind='linear')

 TT,VV   = removeNaN(acePEpoch,aceData['Vx'])
 aceVxmod= epochMatch(uniDateList,TT,VV,interpKind='linear')

 aceXmod = epochMatch(uniDateList,aceBEpoch,aceData['SCxGSE'][sEpochBIDA:eEpochBIDA],interpKind='linear')
 impXmod = epochMatch(uniDateList,impBEpoch,imp8Data['SCxGSE'][sEpochBIDI:eEpochBIDI],interpKind='linear')

 aceYmod = epochMatch(uniDateList,aceBEpoch,aceData['SCyGSE'][sEpochBIDA:eEpochBIDA],interpKind='linear')
 impYmod = epochMatch(uniDateList,impBEpoch,imp8Data['SCyGSE'][sEpochBIDI:eEpochBIDI],interpKind='linear')

 aceZmod = epochMatch(uniDateList,aceBEpoch,aceData['SCzGSE'][sEpochBIDA:eEpochBIDA],interpKind='linear')
 impZmod = epochMatch(uniDateList,impBEpoch,imp8Data['SCzGSE'][sEpochBIDI:eEpochBIDI],interpKind='linear')

 for ind in range(len(aceVxmod)):
  if abs(aceVxmod[ind]) > 0.0: break

 imp8DataMD = {'epoch':uniDateList[ind:-ind]}
 imp8DataMD['V'] = swMedFilter(imp8DataMD['epoch'],impVmod[ind:-ind],15*60)
 imp8DataMD['N'] = swMedFilter(imp8DataMD['epoch'],impNmod[ind:-ind],15*60)
 imp8DataMD['Bz'] = swMedFilter(imp8DataMD['epoch'],impBzmod[ind:-ind],15*60)
 imp8DataMD['SCxGSE'] = impXmod[ind:-ind]
 imp8DataMD['SCyGSE'] = impYmod[ind:-ind]
 imp8DataMD['SCzGSE'] = impZmod[ind:-ind]

 aceDataMD = {'epoch':uniDateList[ind:-ind]}
 aceDataMD['V'] = swMedFilter(aceDataMD['epoch'],aceVmod[ind:-ind],15*60)
 aceDataMD['Vx'] = swMedFilter(aceDataMD['epoch'],aceVxmod[ind:-ind],15*60)
 aceDataMD['N'] = swMedFilter(aceDataMD['epoch'],aceNmod[ind:-ind],15*60)
#aceDataMD['Bx'] = swMedFilter(aceDataMD['epoch'],aceBxmod[ind:-ind],15*60)
#aceDataMD['By'] = swMedFilter(aceDataMD['epoch'],aceBymod[ind:-ind],15*60)
 aceDataMD['Bz'] = swMedFilter(aceDataMD['epoch'],aceBzmod[ind:-ind],15*60)
 aceDataMD['SCxGSE'] = aceXmod[ind:-ind]
 aceDataMD['SCyGSE'] = aceYmod[ind:-ind]
 aceDataMD['SCzGSE'] = aceZmod[ind:-ind]

#plt.figure(100+i+1)
#plt.subplot(3,1,1)
#plt.plot(imp8DataMD['epoch'],imp8DataMD['V'])
#plt.plot( aceDataMD['epoch'], aceDataMD['V'])
#plt.subplot(3,1,2)
#plt.plot(imp8DataMD['epoch'],imp8DataMD['N'])
#plt.plot( aceDataMD['epoch'], aceDataMD['N'])
#plt.subplot(3,1,3)
#plt.plot(imp8DataMD['epoch'],imp8DataMD['Bz'])
#plt.plot( aceDataMD['epoch'], aceDataMD['Bz'])

 destPos = {'X':imp8DataMD['SCxGSE'],'Y':imp8DataMD['SCyGSE'],'Z':imp8DataMD['SCzGSE']}
 lagging, corEpoch = getTimeLag(aceDataMD,destPos,method='standard')
 lagging = dataClean(lagging,[1e-12],['<='])

#plt.figure(300+i+1)
#plt.plot(corEpoch,lagging)
#plt.title('Time Lag between ACE and IMP8 as a function of Epoch time')
#plt.xlabel('Epoch')
#plt.ylabel('Time Lag (seconds)')

#sEpoch1ID  = bisect.bisect_left(imp8DataMD['epoch'], corEpoch[0])
#eEpoch1ID  = bisect.bisect_left(imp8DataMD['epoch'],imp8DataMD['epoch'][-1])
#sEpoch2ID  = 0
#eEpoch2ID  = eEpoch1ID-sEpoch1ID

#SCMatchData={'aceEpoch':corEpoch[sEpoch2ID:eEpoch2ID],'impEpoch':imp8DataMD['epoch'][sEpoch1ID:eEpoch1ID]}
#SCMatchData['aceV']  = aceVmod[sEpoch2ID:eEpoch2ID]
#SCMatchData['impV']  = impVmod[sEpoch1ID:eEpoch1ID]
#SCMatchData['aceN']  = aceNmod[sEpoch2ID:eEpoch2ID]
#SCMatchData['impN']  = impNmod[sEpoch1ID:eEpoch1ID]
#SCMatchData['aceBz'] = aceBzmod[sEpoch2ID:eEpoch2ID]
#SCMatchData['impBz'] = impBzmod[sEpoch1ID:eEpoch1ID]

 cEpoch, Param = epochShift(imp8DataMD['epoch'][ind:-ind],impVmod[ind:-ind],lagging)
 SCMatchData          = {'impEpoch':cEpoch}
 SCMatchData['impV']  = Param
 cEpoch, Param = epochShift(imp8DataMD['epoch'][ind:-ind],impNmod[ind:-ind],lagging)
 SCMatchData['impN']  = Param
 cEpoch, Param = epochShift(imp8DataMD['epoch'][ind:-ind],impBzmod[ind:-ind],lagging)
 SCMatchData['impBz'] = Param

 SCMatchData['aceEpoch'] = cEpoch
 SCMatchData['aceV']     = epochMatch(cEpoch,aceDataMD['epoch'],aceVmod[ind:-ind],interpKind='linear')
 SCMatchData['aceN']     = epochMatch(cEpoch,aceDataMD['epoch'],aceNmod[ind:-ind],interpKind='linear')
 SCMatchData['aceBz']    = epochMatch(cEpoch,aceDataMD['epoch'],aceBzmod[ind:-ind],interpKind='linear')

#plt.figure(400+i+1)
#plt.subplot(3,1,1)
#plt.plot(SCMatchData['aceEpoch'],SCMatchData['aceV'],label='propagated ACE')
#plt.plot(SCMatchData['impEpoch'],SCMatchData['impV'],label='IMP8')
#plt.legend(loc='best')
#plt.ylabel('Solar Wind Speed (km/s)')
#plt.title('Solar Wind - Time Lag')
#plt.subplot(3,1,2)
#plt.plot(SCMatchData['aceEpoch'],SCMatchData['aceN'],label='Propagated ACE')
#plt.plot(SCMatchData['impEpoch'],SCMatchData['impN'],label='IMP8')
#plt.legend(loc='best')
#plt.ylabel('Solar Wind Density (N/cc)')
#plt.subplot(3,1,3)
#plt.plot(SCMatchData['aceEpoch'],SCMatchData['aceBz'],label='Propagated ACE')
#plt.plot(SCMatchData['impEpoch'],SCMatchData['impBz'],label='IMP8')
#plt.legend(loc='best')
#plt.ylabel('$B_z$ (nT)')
#plt.suptitle('Solar Wind Data in' + str(sTime[i]) + '-' + str(eTime[i]))

#plt.figure(500+i+1)
#plt.subplot(3,1,1)
#vCorr, vLag = xcorr(SCMatchData['impV'],SCMatchData['aceV'])
#plt.plot(vLag, vCorr)
#plt.ylabel('Velocity Correlation')
#plt.title('Correlation between Solar Wind Parameters at ACE (after lagging) and IMP8')
#plt.subplot(3,1,2)
#nCorr, nLag = xcorr(SCMatchData['impN'],SCMatchData['aceN'])
#plt.plot(nLag, nCorr)
#plt.ylabel('Density Correlation')
#plt.subplot(3,1,3)
#bCorr, bLag = xcorr(SCMatchData['impBz'],SCMatchData['aceBz'])
#plt.plot(bLag, bCorr)
#plt.ylabel('Bz Correlation')
#plt.xlabel('Index')

 accSrcParam.extend(SCMatchData['aceV'])
 accDesParam.extend(SCMatchData['impV'])

nPins = 100
vRanges = [[251,300],[301,350],[351,400],[401,450],[451,500],[501,550],[551,600]]
vRanges = vRanges + [[601,650],[651,700],[701,750],[751,800],[801,850],[851,900]]
vRanges = vRanges + [[901,950],[951,1000],[1001,1050],[1051,1100],[1101,1150],[1151,1200]]
DesRanges, DesKDE = getDesKDE(accSrcParam,accDesParam,vRanges,nPins)

for j in range(len(DesKDE)):
 if DesKDE[j] != []:
  xVals = numpy.linspace(min(DesRanges[j]),max(DesRanges[j]),nPins)
  plt.figure(600+j+1)
  plt.plot(xVals,DesKDE[j],label=str(sTime[i]))
  plt.suptitle('Source Speed Ranges = [' + str(vRanges[j][0]) + ',' + str(vRanges[j][1]) + '] (km/s)')
  plt.title('KDE of Propagated Solar Wind for 6 hours Slots')
  plt.xlabel('Solar Wind Speed at Destination (km/s)')
  plt.ylabel('Kernel Density Estimation')
  plt.legend(loc='best')

plt.show()

sys.exit()


