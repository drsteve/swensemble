import os
import sys
sys.path.append('/home/ehab/MyFiles/Softex/spacePy/spacepy-0.1.5')
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
from swdatanal import kdeBW, getDesKDE, epochShift, findContiguousData
from swdatanal import search, epochBlock, findCorrEpoch, dataFilter
from swdatanal import getSolarWindType, getTimeLag, swMedFilter, rejectOutliers
from getswdata import getOMNIfiles, getOMNIdata, getOMNIparams, omniDataAdjust
from getswdata import getACEfiles, getACEdata, getACEparams, aceDataAdjust
from getswdata import getIMP8files,getIMP8data, getIMP8params, imp8DataAdjust
from getswdata import getWINDfiles,getWINDdata, getWINDparams, windDataAdjust
from getswdata import getGeotailfiles,getGeotaildata, getGeotailparams, geotailDataAdjust
from getswdata import dataClean, dateShift, dateList, epochMatch, removeNaN

startDate   = datetime.datetime(1999, 1, 1, 0, 0, 0)
endDate     = datetime.datetime(1999, 1,31,23,59,59)
locDateList = dateList(startDate, endDate, shift = 'hour')


dataSRC = ""
swCats = False

#dataSRC = 'geotail'
if dataSRC == 'geotail':
 geotailParams = getGeotailparams(locDateList,'/home/ehab/SWData',geotailSet='1min')
 if geotailParams != '--empty--':
  swParams = ['V','N','T','BX_GSE','BY_GSE','BZ_GSE','ABS_B','X','Y','Z','VX_GSE','VY_GSE','VZ_GSE']
  geotailData = getGeotaildata(locDateList,'/home/ehab/SWData',swParams,geotailSet='1min',dataStat='clean')
  geotailData = geotailDataAdjust(geotailData)

#dataSRC = 'wind'
if dataSRC == 'wind':
 windParams = getWINDparams(locDateList,'/home/ehab/SWData',windSet='1min',windDat='merged')
#windParams = getWINDparams(locDateList,'/home/ehab/SWData',windSet='1min',windDat='magnetic')
#windParams = getWINDparams(locDateList,'/home/ehab/SWData',windSet='1min',windDat='plasma')
 if windParams != '--empty--':
  RE = 6371.0
  swParams = ['BX','BY','BZ','Proton_Np_nonlin','Proton_V_nonlin','Proton_VX_nonlin','Proton_VY_nonlin','Proton_VZ_nonlin','Proton_W_nonlin','xgse','ygse','zgse']
  windData = getWINDdata(locDateList,'/home/ehab/SWData',swParams,windSet='1min',windDat='merged',dataStat='clean')
  windData = windDataAdjust(windData)

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
 #aceData  = getACEdata(locDateList,'/home/ehab/SWData',swParams,['hourly','hourly'],dataStat='clean')
  aceData  = aceDataAdjust(aceData)

  if swCats:
   _, _, _, SWClass = getSolarWindType(aceData)
   aceVEJT=[];aceVCHO=[];aceVSRR=[];aceVSBO=[]
   aceNEJT=[];aceNCHO=[];aceNSRR=[];aceNSBO=[]
   aceBEJT=[];aceBCHO=[];aceBSRR=[];aceBSBO=[]
   acePEEJT=[];acePECHO=[];acePESRR=[];acePESBO=[]
   aceBEEJT=[];aceBECHO=[];aceBESRR=[];aceBESBO=[]
   for item in range(len(SWClass)):
    if SWClass[item] == 'EJT':
     aceVEJT.extend([aceData['V'][item]])
     aceNEJT.extend([aceData['N'][item]])
     aceBEJT.extend([aceData['B'][item]])
     acePEEJT.extend([aceData['plasmaEpoch'][item]])
     aceBEEJT.extend([aceData['magneticEpoch'][item]])
    elif SWClass[item] == 'CHO':
     aceVCHO.extend([aceData['V'][item]])
     aceNCHO.extend([aceData['N'][item]])
     aceBCHO.extend([aceData['B'][item]])
     acePECHO.extend([aceData['plasmaEpoch'][item]])
     aceBECHO.extend([aceData['magneticEpoch'][item]])
    elif SWClass[item] == 'SRR':
     aceVSRR.extend([aceData['V'][item]])
     aceNSRR.extend([aceData['N'][item]])
     aceBSRR.extend([aceData['B'][item]])
     acePESRR.extend([aceData['plasmaEpoch'][item]])
     aceBESRR.extend([aceData['magneticEpoch'][item]])
    elif SWClass[item] == 'SBO':
     aceVSBO.extend([aceData['V'][item]])
     aceNSBO.extend([aceData['N'][item]])
     aceBSBO.extend([aceData['B'][item]])
     acePESBO.extend([aceData['plasmaEpoch'][item]])
     aceBESBO.extend([aceData['magneticEpoch'][item]])
 
   plt.figure(1)
   plt.subplot(3,1,1)
   plt.plot(acePEEJT,aceVEJT, 'bo', label = 'Ejecta')
   plt.plot(acePECHO,aceVCHO, 'ro', label = 'Conronal-Hole')
   plt.plot(acePESRR,aceVSRR, 'mo', label = 'Sector-Reversal')
   plt.plot(acePESBO,aceVSBO, 'go', label = 'Streamer-Belt')
   plt.ylabel('Flow Velocity (km/s)')
   plt.title('ACE Data Solar Wind Categorization March 2003')
   plt.legend(loc='best')
   plt.subplot(3,1,2)
   plt.plot(acePEEJT,aceNEJT, 'bo', label = 'Ejecta')
   plt.plot(acePECHO,aceNCHO, 'ro', label = 'Conronal-Hole')
   plt.plot(acePESRR,aceNSRR, 'mo', label = 'Sector-Reversal')
   plt.plot(acePESBO,aceNSBO, 'go', label = 'Streamer-Belt')
   plt.ylabel('Flow Density (N/cc)')
   plt.legend(loc='best')
   plt.subplot(3,1,3)
   plt.plot(aceBEEJT,aceBEJT, 'bo', label = 'Ejecta')
   plt.plot(aceBECHO,aceBCHO, 'ro', label = 'Conronal-Hole')
   plt.plot(aceBESRR,aceBSRR, 'mo', label = 'Sector-Reversal')
   plt.plot(aceBESBO,aceBSBO, 'go', label = 'Streamer-Belt')
   plt.ylabel('$B_z$ (nT)')
   plt.xlabel('Date')
   plt.legend(loc='best')
   plt.show()

#dataSRC = 'omni'
if dataSRC == 'omni':
 omniParams = getOMNIparams(locDateList,'/home/ehab/SWData',omniSet='hourly')
 if omniParams != '--empty--':
  swParams = ['V','N','T','BX_GSE','BY_GSE','BZ_GSE','ABS_B']
 #swParams = ['proton_density','flow_speed','T','BX_GSE','BY_GSE','BZ_GSE','Vx','Vy','Vz']
  omniData  = getOMNIdata(locDateList,'/home/ehab/SWData',swParams,omniSet='hourly',dataStat='clean')
  omniData  = omniDataAdjust(omniData)

  if swCats:
   _, _, _, SWClass = getSolarWindType(omniData)
   omniVEJT=[];omniVCHO=[];omniVSRR=[];omniVSBO=[]
   omniNEJT=[];omniNCHO=[];omniNSRR=[];omniNSBO=[]
   omniBEJT=[];omniBCHO=[];omniBSRR=[];omniBSBO=[]
   omniEEJT=[];omniECHO=[];omniESRR=[];omniESBO=[]
   for item in range(len(SWClass)):
    if SWClass[item] == 'EJT':
     omniVEJT.extend([omniData['V'][item]])
     omniNEJT.extend([omniData['N'][item]])
     omniBEJT.extend([omniData['B'][item]])
     omniEEJT.extend([omniData['epoch'][item]])
    elif SWClass[item] == 'CHO':
     omniVCHO.extend([omniData['V'][item]])
     omniNCHO.extend([omniData['N'][item]])
     omniBCHO.extend([omniData['B'][item]])
     omniECHO.extend([omniData['epoch'][item]])
    elif SWClass[item] == 'SRR':
     omniVSRR.extend([omniData['V'][item]])
     omniNSRR.extend([omniData['N'][item]])
     omniBSRR.extend([omniData['B'][item]])
     omniESRR.extend([omniData['epoch'][item]])
    elif SWClass[item] == 'SBO':
     omniVSBO.extend([omniData['V'][item]])
     omniNSBO.extend([omniData['N'][item]])
     omniBSBO.extend([omniData['B'][item]])
     omniESBO.extend([omniData['epoch'][item]])
 
   plt.figure(1)
   plt.subplot(3,1,1)
   plt.plot(omniEEJT,omniVEJT, 'bo', label = 'Ejecta')
   plt.plot(omniECHO,omniVCHO, 'ro', label = 'Conronal-Hole')
   plt.plot(omniESRR,omniVSRR, 'mo', label = 'Sector-Reversal')
   plt.plot(omniESBO,omniVSBO, 'go', label = 'Streamer-Belt')
   plt.ylabel('Flow Velocity (km/s)')
   plt.legend(loc='best')
   plt.subplot(3,1,2)
   plt.plot(omniEEJT,omniNEJT, 'bo', label = 'Ejecta')
   plt.plot(omniECHO,omniNCHO, 'ro', label = 'Conronal-Hole')
   plt.plot(omniESRR,omniNSRR, 'mo', label = 'Sector-Reversal')
   plt.plot(omniESBO,omniNSBO, 'go', label = 'Streamer-Belt')
   plt.ylabel('Flow Density (N/cc)')
   plt.legend(loc='best')
   plt.subplot(3,1,3)
   plt.plot(omniEEJT,omniBEJT, 'bo', label = 'Ejecta')
   plt.plot(omniECHO,omniBCHO, 'ro', label = 'Conronal-Hole')
   plt.plot(omniESRR,omniBSRR, 'mo', label = 'Sector-Reversal')
   plt.plot(omniESBO,omniBSBO, 'go', label = 'Streamer-Belt')
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

for iBlock in refBlockStart:
 sEpochID  = bisect_left(imp8Data['plasmaEpoch'], iBlock + timedelta(0,0))
 eEpochID  = bisect_left(imp8Data['plasmaEpoch'], iBlock + timedelta(0,timeShift*3600))
 if len(dataFilter(imp8Data['SCxGSE'][sEpochID:eEpochID],0.0,'<')) >= 1: refBlockStart.remove(iBlock)

sTime = []; eTime = []
for ind in range(len(refBlockStart)):
 sTime.extend([refBlockStart[ind]])
 eTime.extend([dateShift(refBlockStart[ind], hours = timeShift)])

accEpoch = []

accSrcParam = []
accDesParam = []

accSrcV = []; accSrcN = []; accSrcT = []; accSrcB = []
accDesV = []; accDesN = []; accDesT = []; accDesB = []

for i in range(len(sTime)):
 print sTime[i]
 uniDateList = dateList(sTime[i], eTime[i], shift = 'minute')
 accEpoch.extend(uniDateList)

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
 accSrcV.extend(aceVmod)
 TT,VV   = removeNaN(impPEpoch,imp8Data['V'][sEpochPIDI:eEpochPIDI])
 impVmod = epochMatch(uniDateList,TT,VV,interpKind='linear')
 accDesV.extend(impVmod)

 TT,NN   = removeNaN(acePEpoch,aceData['N'][sEpochPIDA:eEpochPIDA])
 aceNmod = epochMatch(uniDateList,TT,NN,interpKind='linear')
 accSrcN.extend(aceNmod)
 TT,NN   = removeNaN(impPEpoch,imp8Data['N'][sEpochPIDI:eEpochPIDI])
 impNmod = epochMatch(uniDateList,TT,NN,interpKind='linear')
 accDesN.extend(impNmod)

 TT,VV   = removeNaN(acePEpoch,aceData['T'][sEpochPIDA:eEpochPIDA])
 aceTmod = epochMatch(uniDateList,TT,VV,interpKind='linear')
 accSrcT.extend(aceTmod)
 TT,VV   = removeNaN(impPEpoch,imp8Data['T'][sEpochPIDI:eEpochPIDI])
 impTmod = epochMatch(uniDateList,TT,VV,interpKind='linear')
 accDesT.extend(impTmod)

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

 TT,BB   = removeNaN(aceBEpoch,aceData['B'][sEpochBIDA:eEpochBIDA])
 aceBmod= epochMatch(uniDateList,TT,BB,interpKind='linear')
 accSrcB.extend(aceBmod)
 TT,BB   = removeNaN(impBEpoch,imp8Data['B'][sEpochBIDI:eEpochBIDI])
 impBmod= epochMatch(uniDateList,TT,BB,interpKind='linear')
 accDesB.extend(impBmod)

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

swCatData = {'epoch':accEpoch,'V':accSrcV,'N':accSrcN,'T':accSrcT,'B':accSrcB}
swCats = True
if swCats:
 swClass = getSolarWindType(swCatData)

sys.exit()

nPins = 100
vRanges = [[251,300],[301,350],[351,400],[401,450],[451,500],[501,550],[551,600]]
vRanges = vRanges + [[601,650],[651,700],[701,750],[751,800],[801,850],[851,900]]
vRanges = vRanges + [[901,950],[951,1000],[1001,1050],[1051,1100],[1101,1150],[1151,1200]]
DesRanges, DesKDE = getDesKDE(accSrcParam,accDesParam,vRanges,nPins)

for j in range(len(DesKDE)):
 if DesKDE[j] != []:
  xVals = numpy.linspace(min(DesRanges[j]),max(DesRanges[j]),nPins)
  plt.figure(600+j+1)
  plt.plot(xVals,DesKDE[j],label='Raw KDE')
  plt.suptitle('Source Speed Ranges = [' + str(vRanges[j][0]) + ',' + str(vRanges[j][1]) + '] (km/s)')
  plt.title('KDE of Propagated Solar Wind for 6 hours Slots')
  plt.xlabel('Solar Wind Speed at Destination (km/s)')
  plt.ylabel('Kernel Density Estimation')

  xvals = numpy.arange(nPins)
  dist_names = ['cauchy']
  for dist_name in dist_names:
    dist = getattr(scipy.stats, dist_name)
    param = dist.fit(DesKDE[j])
    pdf_fitted = dist.pdf(xvals, *param[:-2], loc=param[-2], scale=param[-1])
   #pdf_fitted = dist.pdf(xvals, *param[:-2], loc=param[-2], scale=param[-1]) * len(xVals)
    plt.plot(xVals,pdf_fitted, label=dist_name)

  plt.legend(loc='best')

plt.figure(600)
plt.plot(accEpoch,accSrcV,'bo')
plt.plot(accEpoch,accDesV,'ro')

plt.show()

sys.exit()


