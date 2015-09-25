import os
import sys
sys.path.append('/home/ehab/MyFiles/Softex/spacePy/spacepy-0.1.5')

from matplotlib.backends.backend_pdf import *

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
from scipy.signal import medfilt

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.integrate import simps, trapz
from scipy.stats import ks_2samp, pearsonr, gaussian_kde

from swdatanal import getDistrib, getIndices, omniDataCorr, ccorr, xcorr
from swdatanal import kdeBW, getDesKDE, getSWPRange, getSurrogate
from swdatanal import search, epochBlock, findCorrEpoch, dataFilter
from swdatanal import getSolarWindType, getTimeLag, swMedFilter, rejectOutliers
from getswdata import getOMNIfiles, getOMNIdata, getOMNIparams, omniDataAdjust
from getswdata import getACEfiles, getACEdata, getACEparams, aceDataAdjust
from getswdata import getIMP8files,getIMP8data, getIMP8params, imp8DataAdjust
from getswdata import getWINDfiles,getWINDdata, getWINDparams, windDataAdjust
from getswdata import getGeotailfiles,getGeotaildata, getGeotailparams, geotailDataAdjust
from getswdata import dataClean, dateShift, dateList, mapDataToEpoch, epochShift, commonEpoch, removeNaN


startDate   = datetime.datetime(1998, 1, 1, 0, 0, 0)
endDate     = datetime.datetime(2000,12,31,23,59,59)
locDateList = dateList(startDate, endDate, shift = 'hour')

GeotailDataFlag = False
WindDataFlag = False
IMP8DataFlag = True
ACEDataFlag = True
OMNIDataFlag = False

if GeotailDataFlag:
 dataSRC = 'geotail'

#swParams = ['V','N','T','BX_GSE','BY_GSE','BZ_GSE','ABS_B','X','Y','Z','VX_GSE','VY_GSE','VZ_GSE']
#geotailData = getGeotaildata(locDateList,'/home/ehab/SWData',swParams,geotailSet='1min',dataStat='clean')
#geotailData = geotailDataAdjust(geotailData)

if WindDataFlag:
 dataSRC = 'wind'
 windParams = getWINDparams(locDateList,'/home/ehab/SWData',windSet='1min',windDat='merged')

#swParams = ['BX','BY','BZ','Proton_Np_nonlin','Proton_V_nonlin','Proton_VX_nonlin','Proton_VY_nonlin','Proton_VZ_nonlin','Proton_W_nonlin','xgse','ygse','zgse']
#windData = getWINDdata(locDateList,'/home/ehab/SWData',swParams,windSet='1min',windDat='merged',dataStat='clean')
#windData = windDataAdjust(windData)

if IMP8DataFlag:
 dataSRC = 'imp8'
#imp8Params = getIMP8params(locDateList,'/home/ehab/SWData','15sec',imp8Dat='magnetic')
#imp8Params = getIMP8params(locDateList,'/home/ehab/SWData','1min',imp8Dat='plasma')

 swParams = ['B_Vector_GSE','proton_density_fit','V_fit','protonV_thermal_fit','protonV_thermal_mom','SC_Pos_GSE','SC_Pos_GSM']
 print('Reading IMP8 Dataset')
 imp8Data = getIMP8data(locDateList,'/home/ehab/SWData',swParams,imp8Set='15sec',dataStat='clean')
 print('Adjusting IMP8 Dataset')
 imp8Data = imp8DataAdjust(imp8Data)

if ACEDataFlag:
 dataSRC = 'ace'
#aceParams = getACEparams(locDateList,'/home/ehab/SWData',aceSet='16sec',aceDat='magnetic')
#aceParams = getACEparams(locDateList,'/home/ehab/SWData',aceSet='1min',aceDat='plasma')

 swParams = ['Magnitude','Np','Vp','Tpr','BGSEc','V_GSE','SC_pos_GSE']
 print('Reading ACE Dataset')
 aceData  = getACEdata(locDateList,'/home/ehab/SWData',swParams,['16sec','1min'],dataStat='clean')
 print('Adjusting ACE Dataset')
 aceData  = aceDataAdjust(aceData)
#swClass = getSolarWindType(aceData)

if OMNIDataFlag:
 dataSRC = 'omni'
#omniParams = getOMNIparams(locDateList,'/home/ehab/SWData',omniSet='hourly')

#swParams = ['proton_density','flow_speed','T','ABS_B','BX_GSE','BY_GSE','BZ_GSE','Vx','Vy','Vz']
 swParams = ['V','N','T','ABS_B','BX_GSE','BY_GSE','BZ_GSE']
 omniData  = getOMNIdata(locDateList,'/home/ehab/SWData',swParams,omniSet='hourly',dataStat='clean')
 omniData  = omniDataAdjust(omniData)
 swClass = getSolarWindType(omniData,nCats=3)

'''
pp = PdfPages('./figures/TMP/plot.pdf')
fig = plt.figure(100)
ax=plt.subplot(4,1,1)
plt.plot(imp8Data['plasmaEpoch'],imp8Data['V'])
plt.plot( aceData['plasmaEpoch'], aceData['V'])
plt.title('Solar Wind Parameters at IMP8')
plt.ylabel('$\upsilon_p$ (km/s)')
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax=plt.subplot(4,1,2)
plt.plot(imp8Data['plasmaEpoch'],imp8Data['N'])
plt.plot( aceData['plasmaEpoch'], aceData['N'])
plt.ylabel('$n_p$ (N/cc)')
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax=plt.subplot(4,1,3)
plt.plot(imp8Data['plasmaEpoch'],imp8Data['T'])
plt.plot( aceData['plasmaEpoch'], aceData['T'])
plt.ylabel('$T_p$ (K)')
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax=plt.subplot(4,1,4)
plt.plot(imp8Data['magneticEpoch'],imp8Data['Bz'])
plt.plot( aceData['magneticEpoch'], aceData['Bz'])
plt.ylabel('$B_z$ (nT)')
pp.savefig(fig)
pp.close()

RE = 6371.0
plt.figure(200,figsize=(20,10))
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
'''

timeShift = 4
timeGap = 900

refBlockStart = []
print('Filtering IMP8 Valid Epoch')
pBlockStart = epochBlock(imp8Data['plasmaEpoch'], imp8Data['V'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(pBlockStart,pBlockStart)
pBlockStart = epochBlock(imp8Data['plasmaEpoch'], imp8Data['N'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(imp8Data['plasmaEpoch'], imp8Data['T'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)
pBlockStart = epochBlock(imp8Data['magneticEpoch'], imp8Data['Bz'], blockLen = timeShift, gapLen = timeGap, fixStep = True)
refBlockStart = findCorrEpoch(refBlockStart,pBlockStart)

print('Filtering ACE Valid Epoch')
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
  EE,VV   = removeNaN(aceData['V'][sEpochPIDA:eEpochPIDA],acePEpoch)
  aceVmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(imp8Data['V'][sEpochPIDI:eEpochPIDI],impPEpoch)
  impVmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(aceData['N'][sEpochPIDA:eEpochPIDA],acePEpoch)
  aceNmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(imp8Data['N'][sEpochPIDI:eEpochPIDI],impPEpoch)
  impNmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(aceData['T'][sEpochPIDA:eEpochPIDA],acePEpoch)
  aceTmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(imp8Data['T'][sEpochPIDI:eEpochPIDI],impPEpoch)
  impTmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(aceData['Bz'][sEpochBIDA:eEpochBIDA],aceBEpoch)
  aceBzmod= mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(imp8Data['Bz'][sEpochBIDI:eEpochBIDI],impBEpoch)
  impBzmod= mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(aceData['B'][sEpochBIDA:eEpochBIDA],aceBEpoch)
  aceBmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
  EE,VV   = removeNaN(imp8Data['B'][sEpochBIDI:eEpochBIDI],impBEpoch)
  impBmod = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  EE,VV   = removeNaN(aceData['Vx'][sEpochPIDA:eEpochPIDA],acePEpoch)
  aceVxmod= mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')

  aceXmod = mapDataToEpoch(uniDateList,aceBEpoch,aceData['SCxGSE'][sEpochBIDA:eEpochBIDA],interpKind='linear')
  impXmod = mapDataToEpoch(uniDateList,impBEpoch,imp8Data['SCxGSE'][sEpochBIDI:eEpochBIDI],interpKind='linear')

  aceYmod = mapDataToEpoch(uniDateList,aceBEpoch,aceData['SCyGSE'][sEpochBIDA:eEpochBIDA],interpKind='linear')
  impYmod = mapDataToEpoch(uniDateList,impBEpoch,imp8Data['SCyGSE'][sEpochBIDI:eEpochBIDI],interpKind='linear')

  aceZmod = mapDataToEpoch(uniDateList,aceBEpoch,aceData['SCzGSE'][sEpochBIDA:eEpochBIDA],interpKind='linear')
  impZmod = mapDataToEpoch(uniDateList,impBEpoch,imp8Data['SCzGSE'][sEpochBIDI:eEpochBIDI],interpKind='linear')

  for ind in range(len(aceVxmod)):
   if abs(aceVxmod[ind]) > 0.0: break

  uniDateList = uniDateList[ind:-ind-1:1]

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

  MedFilterTS = 30
  imp8DataMD = {'epoch':uniDateList}
  imp8DataMD['V'] = swMedFilter(imp8DataMD['epoch'],impVmod,MedFilterTS*60)
  imp8DataMD['N'] = swMedFilter(imp8DataMD['epoch'],impNmod,MedFilterTS*60)
  imp8DataMD['Bz'] = swMedFilter(imp8DataMD['epoch'],impBzmod,MedFilterTS*60)
  imp8DataMD['SCxGSE'] = impXmod
  imp8DataMD['SCyGSE'] = impYmod
  imp8DataMD['SCzGSE'] = impZmod

  aceDataMD = {'epoch':uniDateList}
  aceDataMD['V'] = swMedFilter(aceDataMD['epoch'],aceVmod,MedFilterTS*60)
  aceDataMD['Vx'] = swMedFilter(aceDataMD['epoch'],aceVxmod,MedFilterTS*60)
  aceDataMD['N'] = swMedFilter(aceDataMD['epoch'],aceNmod,MedFilterTS*60)
  aceDataMD['Bz'] = swMedFilter(aceDataMD['epoch'],aceBzmod,MedFilterTS*60)
  aceDataMD['SCxGSE'] = aceXmod
  aceDataMD['SCyGSE'] = aceYmod
  aceDataMD['SCzGSE'] = aceZmod
 except:
  print 'Error'
  continue

 '''
 plt.figure(100+i+1,figsize=(20,10))
 plt.subplot(3,1,1)
 plt.plot(imp8DataMD['epoch'],imp8DataMD['V'])
 plt.plot( aceDataMD['epoch'], aceDataMD['V'])
 plt.subplot(3,1,2)
 plt.plot(imp8DataMD['epoch'],imp8DataMD['N'])
 plt.plot( aceDataMD['epoch'], aceDataMD['N'])
 plt.subplot(3,1,3)
 plt.plot(imp8DataMD['epoch'],imp8DataMD['Bz'])
 plt.plot( aceDataMD['epoch'], aceDataMD['Bz'])
 '''

 destPos = {'X':imp8DataMD['SCxGSE'],'Y':imp8DataMD['SCyGSE'],'Z':imp8DataMD['SCzGSE']}
 lagging, sEpoch = getTimeLag(uniDateList,aceDataMD,destPos,method='flat')
 lagging = dataClean(lagging,[1e-12],['<='])

 cmnEpoch, indEpoch = commonEpoch(uniDateList,sEpoch)

 '''
 fig = plt.figure(300+i+1,figsize=(20,10))
 plt.plot(sEpoch,lagging)
 plt.title('Time Lag between ACE and IMP8 as a function of Epoch time')
 plt.xlabel('Epoch')
 plt.ylabel('Time Lag (seconds)')
 plt.savefig('./figures/lag' + str(301+i))
 plt.close(fig)
 '''

 SCMatchData          = {'impEpoch':cmnEpoch,'aceEpoch':uniDateList[0:len(cmnEpoch)]}
 SCMatchData['impV']  = [impVmod[item] for item in indEpoch]
 SCMatchData['impN']  = [impNmod[item] for item in indEpoch]
 SCMatchData['impT']  = [impTmod[item] for item in indEpoch]
 SCMatchData['impBz'] = [impBzmod[item] for item in indEpoch]
 SCMatchData['aceV']  = aceVmod[0:len(cmnEpoch)] 
 SCMatchData['aceN']  = aceNmod[0:len(cmnEpoch)] 
 SCMatchData['aceT']  = aceTmod[0:len(cmnEpoch)] 
 SCMatchData['aceBz'] = aceBzmod[0:len(cmnEpoch)] 

 accDesE.extend(cmnEpoch)
 accDesV.extend([impVmod[item] for item in indEpoch])
 accDesN.extend([impNmod[item] for item in indEpoch])
 accDesT.extend([impTmod[item] for item in indEpoch])
 accDesB.extend([impBmod[item] for item in indEpoch])

 accSrcE.extend(uniDateList[0:len(cmnEpoch)])
 accSrcV.extend(aceVmod[0:len(cmnEpoch)])
 accSrcN.extend(aceNmod[0:len(cmnEpoch)])
 accSrcT.extend(aceTmod[0:len(cmnEpoch)])
 accSrcB.extend(aceBmod[0:len(cmnEpoch)])

 pp400 = PdfPages('./figures/TMP/MatchACEIMP8.pdf')
 fig = plt.figure(400+i+1,figsize=(20,10))
 plt.subplot(4,1,1)
 plt.plot(SCMatchData['impEpoch'],SCMatchData['aceV'],label='ACE')
 plt.plot(SCMatchData['impEpoch'],SCMatchData['impV'],label='IMP8')
 plt.legend(loc='best')
 plt.ylabel('$\upsilon_p$ (km/s)')
 plt.title('Solar Wind - Time Lag')
 plt.subplot(4,1,2)
 plt.plot(SCMatchData['impEpoch'],SCMatchData['aceN'],label='ACE')
 plt.plot(SCMatchData['impEpoch'],SCMatchData['impN'],label='IMP8')
 plt.legend(loc='best')
 plt.ylabel('$n_p$ (N/cc)')
 plt.subplot(4,1,3)
 plt.plot(SCMatchData['impEpoch'],SCMatchData['aceBz'],label='ACE')
 plt.plot(SCMatchData['impEpoch'],SCMatchData['impBz'],label='IMP8')
 plt.legend(loc='best')
 plt.ylabel('$B_z$ (nT)')
 plt.suptitle('Solar Wind Data in' + str(sTime[i]) + '-' + str(eTime[i]))
 plt.subplot(4,1,4)
 plt.plot(SCMatchData['impEpoch'],lagging[0:len(cmnEpoch)])
 plt.ylabel('Lag (seconds)')
 plt.xlabel('IMP8 Epoch')
 pp400.savefig(fig)
 plt.close(fig)
pp400.close()

'''
 nShift=20
 vCorr, vLag = xcorr(SCMatchData['impV'],SCMatchData['aceV'],shift=nShift)
 nCorr, nLag = xcorr(SCMatchData['impN'],SCMatchData['aceN'],shift=nShift)
 bCorr, bLag = xcorr(SCMatchData['impBz'],SCMatchData['aceBz'],shift=nShift)
 vMaxCorr=[];nMaxCorr=[];bMaxCorr=[]
 vMaxLag=[]; nMaxLag=[]; bMaxLag=[]
 for iCounter in range(len(vCorr)):
  vMaxCorr.extend([max(abs(vCorr[iCounter]))])
  nMaxCorr.extend([max(abs(nCorr[iCounter]))])
  bMaxCorr.extend([max(abs(bCorr[iCounter]))])

  vMaxLag.extend(search(abs(vCorr[iCounter]),vMaxCorr[iCounter]))
  nMaxLag.extend(search(abs(nCorr[iCounter]),nMaxCorr[iCounter]))
  bMaxLag.extend(search(abs(bCorr[iCounter]),bMaxCorr[iCounter]))

 plt.figure(500+i+1,figsize=(20,10))
 plt.subplot(3,1,1)
 plt.plot(vMaxLag,vMaxCorr)
#for iCounter in range(len(vCorr)):
# plt.plot(vLag[iCounter], abs(vCorr[iCounter]),label=str(iCounter*nShift))
 plt.legend(loc='best',fontsize=9)
 plt.ylabel('Velocity Correlation')
 plt.title('Correlation between Solar Wind Parameters at ACE (after lagging) and IMP8')
 plt.subplot(3,1,2)
 plt.plot(nMaxLag,nMaxCorr)
#for iCounter in range(len(nCorr)):
# plt.plot(nLag[iCounter], abs(nCorr[iCounter]),label=str(iCounter*nShift))
 plt.legend(loc='best',fontsize=9)
 plt.ylabel('Density Correlation')
 plt.subplot(3,1,3)
 plt.plot(bMaxLag,bMaxCorr)
#for iCounter in range(len(nCorr)):
# plt.plot(bLag[iCounter], abs(bCorr[iCounter]),label=str(iCounter*nShift))
 plt.legend(loc='best',fontsize=9)
 plt.ylabel('Bz Correlation')
 plt.xlabel('Index')
 plt.suptitle('Time: ' + str(sTime[i]) + '-' + str(eTime[i]))
#plt.savefig('./figures/corr' + str(501+i))
#plt.close(fig)
 plt.show()
'''

srcData = {'epoch':accDesE,'V':accSrcV,'N':accSrcN,'T':accSrcT,'B':accSrcB}
desData = {'epoch':accDesE,'V':accDesV,'N':accDesN,'T':accDesT,'B':accDesB}

#SrcData = {'epoch':accSrcE,'V':accSrcV,'N':accSrcN,'T':accSrcT,'B':accSrcB}
#aceFH = open('aceDF.json', 'w+')
#json.dump(SrcData, aceFH)
#DesData = {'epoch':accDesE,'V':accDesV,'N':accDesN,'T':accDesT,'B':accDesB}
#impFH = open('impDF.json', 'w+')
#json.dump(DesData, impFH)

accSrcE = []; accSrcV  = []; accSrcN = []; accSrcT = []; accSrcB = []
accDesE = []; accDesV  = []; accDesN = []; accDesT = []; accDesB = []

swClassSrc = getSolarWindType(srcData,gplot=False)
swClassDes = getSolarWindType(desData,gplot=False)

'''
fig = plt.figure(1,figsize=(20,10))
plt.plot(swClassSrc['EEJT'],swClassSrc['VEJT'],'bo',label='ACE')
plt.plot(swClassDes['EEJT'],swClassDes['VEJT'],'b^',label='Target')
plt.title('Solar Wind Speed (Ejecta)')
plt.xlabel('epoch time')
plt.ylabel('speed (km/s)')
plt.savefig('./figures/swCatEJT')
plt.close(fig)

fig = plt.figure(2,figsize=(20,10))
plt.plot(swClassSrc['ECHO'],swClassSrc['VCHO'],'ro',label='ACE')
plt.plot(swClassDes['ECHO'],swClassDes['VCHO'],'r^',label='Target')
plt.title('Solar Wind Speed (Coronal-Hole Origin)')
plt.xlabel('epoch time')
plt.ylabel('speed (km/s)')
plt.savefig('./figures/swCatCHO')
plt.close(fig)

fig = plt.figure(3,figsize=(20,10))
plt.plot(swClassSrc['ESRR'],swClassSrc['VSRR'],'mo',label='ACE')
plt.plot(swClassDes['ESRR'],swClassDes['VSRR'],'m^',label='Target')
plt.title('Solar Wind Speed (Sector-Reversal Region)')
plt.xlabel('epoch time')
plt.ylabel('speed (km/s)')
plt.savefig('./figures/swCatSRR')
plt.close(fig)

fig = plt.figure(4,figsize=(20,10))
plt.plot(swClassSrc['ESBO'],swClassSrc['VSBO'],'go',label='ACE')
plt.plot(swClassDes['ESBO'],swClassDes['VSBO'],'g^',label='Target')
plt.title('Solar Wind Speed (Streamer-Belt Origin)')
plt.xlabel('epoch time')
plt.ylabel('speed (km/s)')
plt.savefig('./figures/swCatSBO')
plt.close(fig)
'''

nPins = 100
vThreshold = 200
vRanges = [[251,275],[276,300],[301,325],[326,350],[351,375],[376,400],[401,425],[426,450]]
vRanges = vRanges + [[451,475],[476,500],[501,525],[526,550],[551,575],[576,600]]
vRanges = vRanges + [[601,625],[626,650],[651,675],[676,700],[701,725],[726,750],[751,775],[776,800]]
vRanges = vRanges + [[801,825],[826,850],[851,875],[876,900],[901,925],[926,950],[951,975],[976,1000]]

try:
 DesRangesALL, DesKDEALL, KDEfuncALL = getDesKDE(srcData['V'],desData['V'],vRanges,threshold=vThreshold,nPins=nPins)
 ALLFlag = True
except:
 ALLFlag = False
try:
 DesRangesEJT, DesKDEEJT, KDEfuncEJT = getDesKDE(swClassSrc['VEJT'],swClassDes['VEJT'],vRanges,threshold=vThreshold,nPins=nPins)
 EJTFlag = True
except:
 EJTFlag = False
try:
 DesRangesCHO, DesKDECHO, KDEfuncCHO = getDesKDE(swClassSrc['VCHO'],swClassDes['VCHO'],vRanges,threshold=vThreshold,nPins=nPins)
 CHOFlag = True
except:
 CHOFlag = False
try:
 DesRangesSRR, DesKDESRR, KDEfuncSRR = getDesKDE(swClassSrc['VSRR'],swClassDes['VSRR'],vRanges,threshold=vThreshold,nPins=nPins)
 SRRFlag = True
except:
 SRRFlag = False
try:
 DesRangesSBO, DesKDESBO, KDEfuncSBO = getDesKDE(swClassSrc['VSBO'],swClassDes['VSBO'],vRanges,threshold=vThreshold,nPins=nPins)
 SBOFlag = True
except:
 SBOFlag = False

pp = PdfPages('./figures/TMP/SWCatKDE.pdf')
for j in range(len(DesKDEALL)):
 fig = plt.figure(601+j,figsize=(10,10))
 if EJTFlag and DesKDEEJT[j] != []:
  xVals = numpy.linspace(min(DesRangesEJT[j]),max(DesRangesEJT[j]),nPins)
  plt.plot(xVals,DesKDEEJT[j],color='b',label='Ejecta')
 if CHOFlag and DesKDECHO[j] != []:
  xVals = numpy.linspace(min(DesRangesCHO[j]),max(DesRangesCHO[j]),nPins)
  plt.plot(xVals,DesKDECHO[j],color='r',label='Coronal-Hole')
 if SRRFlag and DesKDESRR[j] != []:
  xVals = numpy.linspace(min(DesRangesSRR[j]),max(DesRangesSRR[j]),nPins)
  plt.plot(xVals,DesKDESRR[j],color='m',label='Sector-Reversal')
 if SBOFlag and DesKDESBO[j] != []:
  xVals = numpy.linspace(min(DesRangesSBO[j]),max(DesRangesSBO[j]),nPins)
  plt.plot(xVals,DesKDESBO[j],color='g',label='Streamer-Belt')
 if ALLFlag and DesKDEALL[j] != []:
  xVals = numpy.linspace(min(DesRangesALL[j]),max(DesRangesALL[j]),nPins)
  plt.plot(xVals,DesKDEALL[j],color='k',label='Uncategorized')
  plt.suptitle('Source Speed Ranges = [' + str(vRanges[j][0]) + ',' + str(vRanges[j][1]) + '] (km/s)', fontsize = 20)
  plt.title('KDE of Propagated Solar Wind for 4 hours Slots', fontsize = 20)
  plt.xlabel('Solar Wind Speed at Destination (km/s)', fontsize = 15)
  plt.ylabel('Kernel Density Estimation (KDE)', fontsize = 15)
  plt.legend(loc = 'best')
  plt.axvline(x=vRanges[j][0])
  plt.axvline(x=vRanges[j][1])
  pp.savefig(fig)
 plt.close(fig)
pp.close()

vDestStd = numpy.zeros(len(vRanges))

pp = PdfPages('./figures/TMP/KDE.pdf')
for j in range(len(DesKDEALL)):
 fig = plt.figure(701+j,figsize=(10,10))
 if ALLFlag and DesKDEALL[j] != []:
  xVals = numpy.linspace(min(DesRangesALL[j]),max(DesRangesALL[j]),nPins)
  vDestStd[j] = numpy.std(xVals)
  plt.plot(xVals,DesKDEALL[j],color='k')
  plt.suptitle('Source Speed Ranges = [' + str(vRanges[j][0]) + ',' + str(vRanges[j][1]) + '] (km/s)', fontsize = 20)
  plt.title('KDE of Propagated Solar Wind for 4 hours Slots', fontsize = 20)
  plt.xlabel('Solar Wind Speed at Destination (km/s)', fontsize = 15)
  plt.ylabel('Kernel Density Estimation (KDE)',fontsize = 15)
  plt.axvline(x=vRanges[j][0])
  plt.axvline(x=vRanges[j][1])
  pp.savefig(fig)
 plt.close(fig)
pp.close()

sDate       = datetime.datetime(2003, 1,25, 0, 0, 0)
eDate       = datetime.datetime(2003, 3,10,23,59,59)
locDateList = dateList(sDate, eDate, shift = 'hour')
if ACEDataFlag:
 dataSRC = 'ace'
 swParams = ['Magnitude','Np','Vp','Tpr','BGSEc','V_GSE','SC_pos_GSE']
 print('Reading ACE Dataset')
 aceData  = getACEdata(locDateList,'/home/ehab/SWData',swParams,['16sec','1min'],dataStat='clean')
 print('Adjusting ACE Dataset')
 aceData  = aceDataAdjust(aceData)


uniDateList = dateList(sDate, eDate, shift = 'minute')
sEpochID    = bisect.bisect_left(aceData['plasmaEpoch'], sDate)
eEpochID    = bisect.bisect_left(aceData['plasmaEpoch'], eDate)
aceEpoch    = numpy.array(aceData['plasmaEpoch'][sEpochID:eEpochID])
EE,VV       = removeNaN(aceData['V'][sEpochID:eEpochID],aceEpoch)
aceVmod     = mapDataToEpoch(uniDateList,EE,VV,interpKind='linear')
aceEE       = uniDateList
MedFilterTS = 30
aceVV       = swMedFilter(aceEE,aceVmod,MedFilterTS*60)

vpstd, vbase, vmstd, vepoch = getSurrogate(vRanges,KDEfuncALL,aceVV,aceEE,nSamples=1000)
pp = PdfPages('./figures/TMP/speed.pdf')
fig = plt.figure(801,figsize=(15,15))
plt.plot(vepoch,vbase,'k-')
plt.hold(True)
plt.fill_between(vepoch,vmstd,vpstd,alpha=0.5,facecolor='blue')
plt.title('Solar Wind Ensemble at IMP8 Jan-Mar 2003')
plt.xlabel('Epoch')
plt.ylabel('Solar Wind Speed')
plt.xticks(rotation='45')
pp.savefig(fig)
plt.close(fig)
pp.close()

#srcDataEJT = {'epoch':swClassSrc['EEJT'],'V':swClassSrc['VEJT'],'N':swClassSrc['NEJT'],'T':swClassSrc['TEJT'],'B':swClassSrc['BEJT']}
#srcDataSBO = {'epoch':swClassSrc['ESBO'],'V':swClassSrc['VSBO'],'N':swClassSrc['NSBO'],'T':swClassSrc['TSBO'],'B':swClassSrc['BSBO']}
#srcDataSRR = {'epoch':swClassSrc['ESRR'],'V':swClassSrc['VSRR'],'N':swClassSrc['NSRR'],'T':swClassSrc['TSRR'],'B':swClassSrc['BSRR']}
#srcDataCHO = {'epoch':swClassSrc['ECHO'],'V':swClassSrc['VCHO'],'N':swClassSrc['NCHO'],'T':swClassSrc['TCHO'],'B':swClassSrc['BCHO']}

#desDataEJT = {'epoch':swClassDes['EEJT'],'V':swClassDes['VEJT'],'N':swClassDes['NEJT'],'T':swClassDes['TEJT'],'B':swClassDes['BEJT']}
#desDataSBO = {'epoch':swClassDes['ESBO'],'V':swClassDes['VSBO'],'N':swClassDes['NSBO'],'T':swClassDes['TSBO'],'B':swClassDes['BSBO']}
#desDataSRR = {'epoch':swClassDes['ESRR'],'V':swClassDes['VSRR'],'N':swClassDes['NSRR'],'T':swClassDes['TSRR'],'B':swClassDes['BSRR']}
#desDataCHO = {'epoch':swClassDes['ECHO'],'V':swClassDes['VCHO'],'N':swClassDes['NCHO'],'T':swClassDes['TCHO'],'B':swClassDes['BCHO']}


sys.exit()


