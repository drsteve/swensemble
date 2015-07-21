import os
import sys
sys.path.append('/home/ehab/MyFiles/Softex/spacePy/spacepy-0.1.5')
import json
import math
import numpy
import scipy
import bisect
from bisect import bisect_left
import scipy.stats
import datetime, time
from datetime import timedelta
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
from scipy.stats import ks_2samp, pearsonr, gaussian_kde

from swdatanal import getDistrib, getIndices, omniDataCorr, ccorr, xcorr
from swdatanal import kdeBW, getDesKDE, findContiguousData
from swdatanal import search, epochBlock, findCorrEpoch, dataFilter
from swdatanal import getSolarWindType, getTimeLag, swMedFilter, rejectOutliers
from getswdata import getOMNIfiles, getOMNIdata, getOMNIparams, omniDataAdjust
from getswdata import getACEfiles, getACEdata, getACEparams, aceDataAdjust
from getswdata import getIMP8files,getIMP8data, getIMP8params, imp8DataAdjust
from getswdata import getWINDfiles,getWINDdata, getWINDparams, windDataAdjust
from getswdata import getGeotailfiles,getGeotaildata, getGeotailparams, geotailDataAdjust
from getswdata import dataClean, dateShift, dateList, mapDataToEpoch, epochShift, commonEpoch, removeNaN

startDate   = datetime.datetime(1999, 1, 1, 0, 0, 0)
endDate     = datetime.datetime(1999, 1,31,23,59,59)
locDateList = dateList(startDate, endDate, shift = 'hour')

dataSRC = ""

#dataSRC = 'geotail'
if dataSRC == 'geotail':
 geotailParams = getGeotailparams(locDateList,'/home/ehab/SWData',geotailSet='1min')
#if geotailParams != '--empty--':
# swParams = ['V','N','T','BX_GSE','BY_GSE','BZ_GSE','ABS_B','X','Y','Z','VX_GSE','VY_GSE','VZ_GSE']
# geotailData = getGeotaildata(locDateList,'/home/ehab/SWData',swParams,geotailSet='1min',dataStat='clean')
# geotailData = geotailDataAdjust(geotailData)

#dataSRC = 'wind'
if dataSRC == 'wind':
 windParams = getWINDparams(locDateList,'/home/ehab/SWData',windSet='1min',windDat='merged')
#if windParams != '--empty--':
# swParams = ['BX','BY','BZ','Proton_Np_nonlin','Proton_V_nonlin','Proton_VX_nonlin','Proton_VY_nonlin','Proton_VZ_nonlin','Proton_W_nonlin','xgse','ygse','zgse']
# windData = getWINDdata(locDateList,'/home/ehab/SWData',swParams,windSet='1min',windDat='merged',dataStat='clean')
# windData = windDataAdjust(windData)

dataSRC = 'imp8'
if dataSRC == 'imp8':
 imp8Params = getIMP8params(locDateList,'/home/ehab/SWData','15sec',imp8Dat='magnetic')
 imp8Params = getIMP8params(locDateList,'/home/ehab/SWData','1min',imp8Dat='plasma')
 if imp8Params != '--empty--':
  swParams = ['B_Vector_GSE','proton_density_fit','V_fit','protonV_thermal_fit','protonV_thermal_mom','SC_Pos_GSE','SC_Pos_GSM']
  imp8Data = getIMP8data(locDateList,'/home/ehab/SWData',swParams,imp8Set='15sec',dataStat='clean')
  imp8Data = imp8DataAdjust(imp8Data)
 #imp8FH = open('imp8DF.json', 'w+')
 #json.dump(imp8Data, imp8FH)

dataSRC = 'ace'
if dataSRC == 'ace':
 aceParams = getACEparams(locDateList,'/home/ehab/SWData',aceSet='16sec',aceDat='magnetic')
 aceParams = getACEparams(locDateList,'/home/ehab/SWData',aceSet='1min',aceDat='plasma')
 if aceParams != '--empty--':
  swParams = ['Magnitude','Np','Vp','Tpr','BGSEc','V_GSE','SC_pos_GSE']
  aceData  = getACEdata(locDateList,'/home/ehab/SWData',swParams,['16sec','1min'],dataStat='clean')
  aceData  = aceDataAdjust(aceData)
 #swClass = getSolarWindType(aceData)
 #aceFH = open('aceDF.json', 'w+')
 #json.dump(aceData, aceFH)

#dataSRC = 'omni'
if dataSRC == 'omni':
 omniParams = getOMNIparams(locDateList,'/home/ehab/SWData',omniSet='1min')
#if omniParams != '--empty--':
# swParams = ['proton_density','flow_speed','T','BX_GSE','BY_GSE','BZ_GSE','Vx','Vy','Vz']
# omniData  = getOMNIdata(locDateList,'/home/ehab/SWData',swParams,omniSet='1min',dataStat='clean')
# omniData  = omniDataAdjust(omniData)
# swClass = getSolarWindType(omniData)


#RE = 6371.0
#plt.figure(200)
#plt.subplot(3,1,1)
#plt.plot(imp8Data['magneticEpoch'],imp8Data['SCxGSE']/RE,label='IMP8')
#plt.plot( aceData['magneticEpoch'], aceData['SCxGSE']/RE,label='ACE')
#plt.ylabel('X-Position ($R_E$)')
#plt.title('ACE and IMP8 Location')
#plt.legend(loc='best')
#plt.subplot(3,1,2)
#plt.plot(imp8Data['magneticEpoch'],imp8Data['SCyGSE']/RE,label='Imp8')
#plt.plot( aceData['magneticEpoch'], aceData['SCyGSE']/RE,label='ACE')
#plt.ylabel('Y-Position ($R_E$)')
#plt.legend(loc='best')
#plt.subplot(3,1,3)
#plt.plot(imp8Data['magneticEpoch'],imp8Data['SCzGSE']/RE,label='IMP8')
#plt.plot( aceData['magneticEpoch'], aceData['SCzGSE']/RE,label='ACE')
#plt.ylabel('Z-Position ($R_E$)')
#plt.xlabel('Epoch Time')
#plt.legend(loc='best')

timeShift = 4
timeGap = 900

pBlockStart = epochBlock(imp8Data['plasmaEpoch'], imp8Data['V'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(pBlockStart,pBlockStart)
pBlockStart = epochBlock(imp8Data['plasmaEpoch'], imp8Data['N'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(imp8Data['magneticEpoch'], imp8Data['Bz'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(imp8Data['plasmaEpoch'], imp8Data['T'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)

pBlockStart = epochBlock(aceData['plasmaEpoch'], aceData['Vx'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(aceData['plasmaEpoch'], aceData['V'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(aceData['plasmaEpoch'], aceData['N'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(aceData['plasmaEpoch'], aceData['T'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(aceData['magneticEpoch'], aceData['Bz'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)

for iBlock in refBlockStart:
 sEpochID  = bisect_left(imp8Data['plasmaEpoch'], iBlock + timedelta(0,0))
 eEpochID  = bisect_left(imp8Data['plasmaEpoch'], iBlock + timedelta(0,timeShift*3600))
 if len(dataFilter(imp8Data['SCxGSE'][sEpochID:eEpochID],0.0,'<')) >= 1: refBlockStart.remove(iBlock)

sTime = []; eTime = []
for ind in range(len(refBlockStart)):
 sTime.extend([refBlockStart[ind]])
 eTime.extend([dateShift(refBlockStart[ind], hours = timeShift)])

accSrcE = []; accSrcV = []; accSrcN = []; accSrcT = []; accSrcB = []
accDesE = []; accDesV = []; accDesN = []; accDesT = []; accDesB = []

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

 try:
  EE,VV   = removeNaN(acePEpoch,aceData['V'][sEpochPIDA:eEpochPIDA])
  aceVmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(impPEpoch,imp8Data['V'][sEpochPIDI:eEpochPIDI])
  impVmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(acePEpoch,aceData['N'][sEpochPIDA:eEpochPIDA])
  aceNmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(impPEpoch,imp8Data['N'][sEpochPIDI:eEpochPIDI])
  impNmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(acePEpoch,aceData['T'][sEpochPIDA:eEpochPIDA])
  aceTmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(impPEpoch,imp8Data['T'][sEpochPIDI:eEpochPIDI])
  impTmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(aceBEpoch,aceData['Bz'][sEpochBIDA:eEpochBIDA])
  aceBzmod= mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(impBEpoch,imp8Data['Bz'][sEpochBIDI:eEpochBIDI])
  impBzmod= mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(aceBEpoch,aceData['B'][sEpochBIDA:eEpochBIDA])
  aceBmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(impBEpoch,imp8Data['B'][sEpochBIDI:eEpochBIDI])
  impBmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(acePEpoch,aceData['Vx'][sEpochPIDA:eEpochPIDA])
  aceVxmod= mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  aceXmod = mapDataToEpoch(uniDateList,aceBEpoch,aceData['SCxGSE'][sEpochBIDA:eEpochBIDA],interpKind='linear')
  impXmod = mapDataToEpoch(uniDateList,impBEpoch,imp8Data['SCxGSE'][sEpochBIDI:eEpochBIDI],interpKind='linear')

  aceYmod = mapDataToEpoch(uniDateList,aceBEpoch,aceData['SCyGSE'][sEpochBIDA:eEpochBIDA],interpKind='linear')
  impYmod = mapDataToEpoch(uniDateList,impBEpoch,imp8Data['SCyGSE'][sEpochBIDI:eEpochBIDI],interpKind='linear')

  aceZmod = mapDataToEpoch(uniDateList,aceBEpoch,aceData['SCzGSE'][sEpochBIDA:eEpochBIDA],interpKind='linear')
  impZmod = mapDataToEpoch(uniDateList,impBEpoch,imp8Data['SCzGSE'][sEpochBIDI:eEpochBIDI],interpKind='linear')
 except:
  print 'Error'
  continue

 for ind in range(len(aceVxmod)):
  if abs(aceVxmod[ind]) > 0.0: break

 uniDateList = uniDateList[ind:-ind:1]

 impVmod     = impVmod[ind:-ind:1]
 impNmod     = impNmod[ind:-ind:1]
 impTmod     = impTmod[ind:-ind:1]
 impBmod     = impBmod[ind:-ind:1]
 impXmod     = impXmod[ind:-ind:1]
 impYmod     = impYmod[ind:-ind:1]
 impZmod     = impZmod[ind:-ind:1]
 impBzmod    = impBzmod[ind:-ind:1]

 aceVmod     = aceVmod[ind:-ind:1]
 aceNmod     = aceNmod[ind:-ind:1]
 aceTmod     = aceTmod[ind:-ind:1]
 aceBmod     = aceBmod[ind:-ind:1]
 aceXmod     = aceXmod[ind:-ind:1]
 aceYmod     = aceYmod[ind:-ind:1]
 aceZmod     = aceZmod[ind:-ind:1]
 aceVxmod    = aceVxmod[ind:-ind:1]
 aceBzmod    = aceBzmod[ind:-ind:1]

 imp8DataMD = {'epoch':uniDateList}
 imp8DataMD['V'] = swMedFilter(imp8DataMD['epoch'],impVmod,15*60)
 imp8DataMD['N'] = swMedFilter(imp8DataMD['epoch'],impNmod,15*60)
 imp8DataMD['Bz'] = swMedFilter(imp8DataMD['epoch'],impBzmod,15*60)
 imp8DataMD['SCxGSE'] = impXmod
 imp8DataMD['SCyGSE'] = impYmod
 imp8DataMD['SCzGSE'] = impZmod

 aceDataMD = {'epoch':uniDateList}
 aceDataMD['V'] = swMedFilter(aceDataMD['epoch'],aceVmod,15*60)
 aceDataMD['Vx'] = swMedFilter(aceDataMD['epoch'],aceVxmod,15*60)
 aceDataMD['N'] = swMedFilter(aceDataMD['epoch'],aceNmod,15*60)
 aceDataMD['Bz'] = swMedFilter(aceDataMD['epoch'],aceBzmod,15*60)
 aceDataMD['SCxGSE'] = aceXmod
 aceDataMD['SCyGSE'] = aceYmod
 aceDataMD['SCzGSE'] = aceZmod

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
 lagging, sEpoch = getTimeLag(uniDateList,aceDataMD,destPos,method='standard')
 lagging = dataClean(lagging,[1e-12],['<='])

 cmnEpoch, indEpoch = commonEpoch(uniDateList,sEpoch)

#plt.figure(300+i+1)
#plt.plot(corEpoch,lagging)
#plt.title('Time Lag between ACE and IMP8 as a function of Epoch time')
#plt.xlabel('Epoch')
#plt.ylabel('Time Lag (seconds)')

 SCMatchData          = {'impEpoch':cmnEpoch,'aceEpoch':uniDateList[0:len(cmnEpoch)]}
 SCMatchData['impV']  = [impVmod[ind] for ind in indEpoch]
 SCMatchData['impN']  = [impNmod[ind] for ind in indEpoch]
 SCMatchData['impT']  = [impTmod[ind] for ind in indEpoch]
 SCMatchData['impBz'] = [impBzmod[ind] for ind in indEpoch]
 SCMatchData['aceV']  = aceVmod[0:len(cmnEpoch)] 
 SCMatchData['aceN']  = aceNmod[0:len(cmnEpoch)] 
 SCMatchData['aceT']  = aceTmod[0:len(cmnEpoch)] 
 SCMatchData['aceBz'] = aceBzmod[0:len(cmnEpoch)] 

 accDesE.extend(cmnEpoch)
 accDesV.extend([impVmod[ind] for ind in indEpoch])
 accDesN.extend([impNmod[ind] for ind in indEpoch])
 accDesT.extend([impTmod[ind] for ind in indEpoch])
 accDesB.extend([impBmod[ind] for ind in indEpoch])

 accSrcE.extend(uniDateList[0:len(cmnEpoch)])
 accSrcV.extend(aceVmod[0:len(cmnEpoch)])
 accSrcN.extend(aceNmod[0:len(cmnEpoch)])
 accSrcT.extend(aceTmod[0:len(cmnEpoch)])
 accSrcB.extend(aceBmod[0:len(cmnEpoch)])

#plt.figure(400+i+1)
#plt.subplot(4,1,1)
#plt.plot(SCMatchData['impEpoch'],SCMatchData['aceV'],label='propagated ACE')
#plt.plot(SCMatchData['impEpoch'],SCMatchData['impV'],label='IMP8')
#plt.legend(loc='best')
#plt.ylabel('Solar Wind Speed (km/s)')
#plt.title('Solar Wind - Time Lag')
#plt.subplot(4,1,2)
#plt.plot(SCMatchData['impEpoch'],SCMatchData['aceN'],label='Propagated ACE')
#plt.plot(SCMatchData['impEpoch'],SCMatchData['impN'],label='IMP8')
#plt.legend(loc='best')
#plt.ylabel('Solar Wind Density (N/cc)')
#plt.subplot(4,1,3)
#plt.plot(SCMatchData['impEpoch'],SCMatchData['aceBz'],label='Propagated ACE')
#plt.plot(SCMatchData['impEpoch'],SCMatchData['impBz'],label='IMP8')
#plt.legend(loc='best')
#plt.ylabel('$B_z$ (nT)')
#plt.suptitle('Solar Wind Data in' + str(sTime[i]) + '-' + str(eTime[i]))
#plt.subplot(4,1,4)
#plt.plot(SCMatchData['impEpoch'],lagging[0:len(cmnEpoch)])
#plt.title('Time Lag between ACE and IMP8 as a function of Epoch time')
#plt.ylabel('Time Lag (seconds)')
#plt.xlabel('IMP8 Epoch')

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

#srcData = {'epoch':accSrcE,'V':accSrcV,'N':accSrcN,'T':accSrcT,'B':accSrcB}
srcData = {'epoch':accDesE,'V':accSrcV,'N':accSrcN,'T':accSrcT,'B':accSrcB}
desData = {'epoch':accDesE,'V':accDesV,'N':accDesN,'T':accDesT,'B':accDesB}

accSrcE = []; accSrcV  = []; accSrcN = []; accSrcT = []; accSrcB = []
accDesE = []; accDesV  = []; accDesN = []; accDesT = []; accDesB = []

a,b = ccorr(accSrcV,accDesV)
plt.plot(b,a)
plt.show()

sys.exit()

swClassSrc = getSolarWindType(srcData,gplot=False)
swClassDes = getSolarWindType(desData,gplot=False)

plt.figure(1)
plt.plot(swClassSrc['EEJT'],swClassSrc['VEJT'],'bo')
plt.plot(swClassDes['EEJT'],swClassDes['VEJT'],'b^')

plt.figure(2)
plt.plot(swClassSrc['ECHO'],swClassSrc['VCHO'],'ro')
plt.plot(swClassDes['ECHO'],swClassDes['VCHO'],'r^')

plt.figure(3)
plt.plot(swClassSrc['ESRR'],swClassSrc['VSRR'],'mo')
plt.plot(swClassDes['ESRR'],swClassDes['VSRR'],'m^')

plt.figure(4)
plt.plot(swClassSrc['ESBO'],swClassSrc['VSBO'],'go')
plt.plot(swClassDes['ESBO'],swClassDes['VSBO'],'g^')

nPins = 100
vRanges = [[201,250],[251,300],[301,350],[351,400],[401,450],[451,500],[501,550],[551,600]]
vRanges = vRanges + [[601,650],[651,700],[701,750],[751,800],[801,850],[851,900]]
vRanges = vRanges + [[901,950],[951,1000],[1001,1050],[1051,1100],[1101,1150],[1151,1200]]

DesRanges, DesKDE = getDesKDE(srcData['V'],desData['V'],vRanges,nPins)
for j in range(len(DesKDE)):
 if DesKDE[j] != []:
  xVals = numpy.linspace(min(DesRanges[j]),max(DesRanges[j]),nPins)
  plt.figure(600+j+1)
  plt.plot(xVals,DesKDE[j],color='k')
  plt.suptitle('Source Speed Ranges = [' + str(vRanges[j][0]) + ',' + str(vRanges[j][1]) + '] (km/s)')
  plt.title('KDE of Propagated Solar Wind for 6 hours Slots')
  plt.xlabel('Solar Wind Speed at Destination (km/s)')
  plt.ylabel('Kernel Density Estimation')
 if DesKDE[j] != []:
  xVals = numpy.linspace(min(DesRanges[j]),max(DesRanges[j]),nPins)
  plt.figure(700+j+1)
  plt.plot(xVals,DesKDE[j],color='k',label='Uncategorized')
  plt.suptitle('Source Speed Ranges = [' + str(vRanges[j][0]) + ',' + str(vRanges[j][1]) + '] (km/s)')
  plt.title('KDE of Propagated Solar Wind for 6 hours Slots')
  plt.xlabel('Solar Wind Speed at Destination (km/s)')
  plt.ylabel('Kernel Density Estimation')

try:
 DesRanges, DesKDE = getDesKDE(swClassSrc['VEJT'],swClassDes['VEJT'],vRanges,nPins)
 for j in range(len(DesKDE)):
  if DesKDE[j] != []:
   xVals = numpy.linspace(min(DesRanges[j]),max(DesRanges[j]),nPins)
   plt.figure(600+j+1)
   plt.plot(xVals,DesKDE[j],color='b',label='Ejecta')
except:
 pass

try:
 DesRanges, DesKDE = getDesKDE(swClassSrc['VCHO'],swClassDes['VCHO'],vRanges,nPins)
 for j in range(len(DesKDE)):
  if DesKDE[j] != []:
   xVals = numpy.linspace(min(DesRanges[j]),max(DesRanges[j]),nPins)
   plt.figure(600+j+1)
   plt.plot(xVals,DesKDE[j],color='r',label='Coronal-Hole')
except:
 pass

try:
 DesRanges, DesKDE = getDesKDE(swClassSrc['VSRR'],swClassDes['VSRR'],vRanges,nPins)
 for j in range(len(DesKDE)):
  if DesKDE[j] != []:
   xVals = numpy.linspace(min(DesRanges[j]),max(DesRanges[j]),nPins)
   plt.figure(600+j+1)
   plt.plot(xVals,DesKDE[j],color='m',label='Sector-Reversal')
except:
 pass

try:
 DesRanges, DesKDE = getDesKDE(swClassSrc['VSBO'],swClassDes['VSBO'],vRanges,nPins)
 for j in range(len(DesKDE)):
  if DesKDE[j] != []:
   xVals = numpy.linspace(min(DesRanges[j]),max(DesRanges[j]),nPins)
   plt.figure(600+j+1)
   plt.plot(xVals,DesKDE[j],color='g',label='Streamer-Belt')
except:
 pass

plt.legend(loc = 'best')
#srcDataEJT = {'epoch':swClassSrc['EEJT'],'V':swClassSrc['VEJT'],'N':swClassSrc['NEJT'],'T':swClassSrc['TEJT'],'B':swClassSrc['BEJT']}
#srcDataSBO = {'epoch':swClassSrc['ESBO'],'V':swClassSrc['VSBO'],'N':swClassSrc['NSBO'],'T':swClassSrc['TSBO'],'B':swClassSrc['BSBO']}
#srcDataSRR = {'epoch':swClassSrc['ESRR'],'V':swClassSrc['VSRR'],'N':swClassSrc['NSRR'],'T':swClassSrc['TSRR'],'B':swClassSrc['BSRR']}
#srcDataCHO = {'epoch':swClassSrc['ECHO'],'V':swClassSrc['VCHO'],'N':swClassSrc['NCHO'],'T':swClassSrc['TCHO'],'B':swClassSrc['BCHO']}

#desDataEJT = {'epoch':swClassDes['EEJT'],'V':swClassDes['VEJT'],'N':swClassDes['NEJT'],'T':swClassDes['TEJT'],'B':swClassDes['BEJT']}
#desDataSBO = {'epoch':swClassDes['ESBO'],'V':swClassDes['VSBO'],'N':swClassDes['NSBO'],'T':swClassDes['TSBO'],'B':swClassDes['BSBO']}
#desDataSRR = {'epoch':swClassDes['ESRR'],'V':swClassDes['VSRR'],'N':swClassDes['NSRR'],'T':swClassDes['TSRR'],'B':swClassDes['BSRR']}
#desDataCHO = {'epoch':swClassDes['ECHO'],'V':swClassDes['VCHO'],'N':swClassDes['NCHO'],'T':swClassDes['TCHO'],'B':swClassDes['BCHO']}


plt.show()

sys.exit()


