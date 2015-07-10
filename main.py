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
import spacepy.time as spt
import matplotlib.pyplot as plt
from scipy.signal import medfilt
from mpl_toolkits.mplot3d import Axes3D

from scipy.integrate import simps, trapz
from scipy.stats import pearsonr, gaussian_kde

from swdatanal import getDistrib, getIndices, omniDataCorr, ccorr, kdeBW
from swdatanal import getSolarWindType, getTimeLag, swMedFilter
from getswdata import getOMNIfiles, getOMNIdata, getOMNIparams
from getswdata import getACEfiles, getACEdata, getACEparams, aceDataAdjust
from getswdata import getIMP8files,getIMP8data, getIMP8params, imp8DataAdjust
from getswdata import getWINDfiles,getWINDdata, getWINDparams
from getswdata import getGeotailfiles,getGeotaildata, getGeotailparams
from getswdata import dataClean, dateShift, dateList, epochMatch, removeNaN


startDate   = datetime.datetime(1999, 7,21, 0, 0, 0)
endDate     = datetime.datetime(1999, 7,28, 0, 0, 0)
locDateList = dateList(startDate, endDate, shift = 'hour')

iCounter = 0
sTime = []; eTime = []
while 1:
 sTime.extend([dateShift(startDate, hours=iCounter*6)])
 eTime.extend([dateShift(startDate, hours=(iCounter+1)*6)])
 if eTime[iCounter] >= endDate: break
 iCounter = iCounter + 1

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

#plt.figure(1)
#plt.plot(imp8Data['magneticEpoch'],imp8Data['SCxGSE'])
#plt.show()

for i in range(len(sTime)):
 sEpochBIDI  = bisect.bisect_left(imp8Data['magneticEpoch'], sTime[i])
 eEpochBIDI  = bisect.bisect_left(imp8Data['magneticEpoch'], eTime[i])
 sEpochPIDI  = bisect.bisect_left(imp8Data['plasmaEpoch'], sTime[i])
 eEpochPIDI  = bisect.bisect_left(imp8Data['plasmaEpoch'], eTime[i])
 impPEpoch   = imp8Data['plasmaEpoch'][sEpochPIDI:eEpochPIDI]
 impBEpoch   = imp8Data['magneticEpoch'][sEpochBIDI:eEpochBIDI]
 
 sEpochBIDA  = bisect.bisect_left(aceData['magneticEpoch'], sTime[i])
 eEpochBIDA  = bisect.bisect_left(aceData['magneticEpoch'], eTime[i])
 sEpochPIDA  = bisect.bisect_left(aceData['plasmaEpoch'], sTime[i])
 eEpochPIDA  = bisect.bisect_left(aceData['plasmaEpoch'], eTime[i])
 acePEpoch   = aceData['plasmaEpoch'][sEpochPIDA:eEpochPIDA]
 aceBEpoch   = aceData['magneticEpoch'][sEpochBIDA:eEpochBIDA]
 
 TT,VV = removeNaN(acePEpoch,aceData['V'][sEpochPIDA:eEpochPIDA])
 aceVmod=epochMatch(impPEpoch,TT,VV)
 TT,VV = removeNaN(impPEpoch,imp8Data['V'][sEpochPIDI:eEpochPIDI])
 impVmod=epochMatch(impPEpoch,TT,VV)

#TT,BB = removeNaN(aceBEpoch,aceData['Bx'][sEpochPIDA:eEpochPIDA])
#aceBxmod=epochMatch(impPEpoch,TT,BB)
#TT,BB = removeNaN(impPEpoch,imp8Data['Bx'][sEpochPIDI:eEpochPIDI])
#impBxmod=epochMatch(impPEpoch,TT,BB)

#TT,BB = removeNaN(aceBEpoch,aceData['By'][sEpochPIDA:eEpochPIDA])
#aceBymod=epochMatch(impPEpoch,TT,BB)
#TT,BB = removeNaN(impPEpoch,imp8Data['By'][sEpochPIDI:eEpochPIDI])
#impBymod=epochMatch(impPEpoch,TT,BB)

#TT,BB = removeNaN(aceBEpoch,aceData['Bz'][sEpochPIDA:eEpochPIDA])
#aceBzmod=epochMatch(impPEpoch,TT,BB)
#TT,BB = removeNaN(impPEpoch,imp8Data['Bz'][sEpochPIDI:eEpochPIDI])
#impBzmod=epochMatch(impPEpoch,TT,BB)

#TT,VV = removeNaN(acePEpoch,aceData['Vx'])
#aceVxmod=epochMatch(imp8Data['plasmaEpoch'][sEpochPIDI:eEpochPIDI],TT,VV)

#TT,NN  = removeNaN(acePEpoch,aceData['N'][sEpochPIDA:eEpochPIDA])
#aceNmod=epochMatch(impPEpoch,TT,NN)
#TT,NN = removeNaN(impPEpoch,imp8Data['N'][sEpochPIDI:eEpochPIDI])
#impNmod=epochMatch(impPEpoch,TT,NN)

 plt.figure(100)

 plt.subplot(4,1,1)
 plt.plot(impPEpoch,impVmod)
 plt.plot(impPEpoch,aceVmod)

#plt.subplot(4,1,2)
#plt.plot(impPEpoch,impNmod)
#plt.plot(impPEpoch,aceNmod)

#plt.subplot(4,1,3)
#plt.plot(impPEpoch,impBzmod)
#plt.plot(impPEpoch,aceBzmod)

#plt.subplot(4,1,4)
#plt.plot(impPEpoch,aceVxmod)

 plt.show()
sys.exit()

aceXmod=epochMatch(imp8Data['plasmaEpoch'],aceData['magneticEpoch'],aceData['SCxGSE'])
impXmod=epochMatch(imp8Data['plasmaEpoch'],imp8Data['magneticEpoch'],imp8Data['SCxGSE'])

aceYmod=epochMatch(imp8Data['plasmaEpoch'],aceData['magneticEpoch'],aceData['SCyGSE'])
impYmod=epochMatch(imp8Data['plasmaEpoch'],imp8Data['magneticEpoch'],imp8Data['SCyGSE'])

aceZmod=epochMatch(imp8Data['plasmaEpoch'],aceData['magneticEpoch'],aceData['SCzGSE'])
impZmod=epochMatch(imp8Data['plasmaEpoch'],imp8Data['magneticEpoch'],imp8Data['SCzGSE'])

imp8DataMD = {'epoch':imp8Data['plasmaEpoch']}
imp8DataMD['V'] = swMedFilter(imp8DataMD['epoch'],impVmod,45*60)
imp8DataMD['N'] = swMedFilter(imp8DataMD['epoch'],impNmod,45*60)
imp8DataMD['Bz'] = swMedFilter(imp8DataMD['epoch'],impBzmod,45*60)
imp8DataMD['SCxGSE'] = impXmod
imp8DataMD['SCyGSE'] = impYmod
imp8DataMD['SCzGSE'] = impZmod

aceDataMD = {'epoch':imp8Data['plasmaEpoch']}
aceDataMD['V'] = swMedFilter(aceDataMD['epoch'],aceVmod,45*60)
aceDataMD['Vx'] = swMedFilter(aceDataMD['epoch'],aceVxmod,45*60)
aceDataMD['N'] = swMedFilter(aceDataMD['epoch'],aceNmod,45*60)
aceDataMD['Bx'] = swMedFilter(aceDataMD['epoch'],aceBxmod,45*60)
aceDataMD['By'] = swMedFilter(aceDataMD['epoch'],aceBymod,45*60)
aceDataMD['Bz'] = swMedFilter(aceDataMD['epoch'],aceBzmod,45*60)
aceDataMD['SCxGSE'] = aceXmod
aceDataMD['SCyGSE'] = aceYmod
aceDataMD['SCzGSE'] = aceZmod

RE = 6371.0
plt.figure(200)
plt.subplot(3,1,1)
plt.plot(imp8DataMD['epoch'],impXmod/RE,label='IMP8')
plt.plot(aceDataMD['epoch'],aceXmod/RE,label='ACE')
plt.ylabel('X-Position ($R_E$)')
plt.title('ACE and IMP8 Location')
plt.legend(loc='best')
plt.subplot(3,1,2)
plt.plot(imp8DataMD['epoch'],impYmod/RE,label='Imp8')
plt.plot(aceDataMD['epoch'],aceYmod/RE,label='ACE')
plt.ylabel('Y-Position ($R_E$)')
plt.legend(loc='best')
plt.subplot(3,1,3)
plt.plot(imp8DataMD['epoch'],impZmod/RE,label='IMP8')
plt.plot(aceDataMD['epoch'],aceZmod/RE,label='ACE')
plt.ylabel('Z-Position ($R_E$)')
plt.xlabel('Epoch Time')
plt.legend(loc='best')

destPos = {'X':impXmod,'Y':impYmod,'Z':impZmod}
lagging, corEpoch = getTimeLag(aceDataMD,destPos,method='standard')
lagging = dataClean(lagging,[1e-12],['<='])

sEpoch1ID  = bisect.bisect_left(imp8DataMD['epoch'], corEpoch[0])
eEpoch1ID  = bisect.bisect_left(imp8DataMD['epoch'],imp8DataMD['epoch'][-1])
sEpoch2ID  = 0
eEpoch2ID  = eEpoch1ID-sEpoch1ID

plt.figure(300)
plt.subplot(4,1,1)
plt.plot(imp8DataMD['epoch'][sEpoch1ID:eEpoch1ID],imp8DataMD['V'][sEpoch1ID:eEpoch1ID],label='IMP8')
plt.plot(aceDataMD['epoch'][sEpoch1ID:eEpoch1ID],aceDataMD['V'][sEpoch1ID:eEpoch1ID],label='ACE')
plt.plot(corEpoch[sEpoch2ID:eEpoch2ID],aceDataMD['V'][sEpoch2ID:eEpoch2ID],label='propagated ACE')
plt.legend(loc='best')
plt.ylabel('Solar Wind Speed (km/s)')
plt.title('Solar Wind - Time Lag')
plt.subplot(4,1,2)
plt.plot(imp8DataMD['epoch'][sEpoch1ID:eEpoch1ID],imp8DataMD['N'][sEpoch1ID:eEpoch1ID],label='IMP8')
plt.plot( aceDataMD['epoch'][sEpoch1ID:eEpoch1ID], aceDataMD['N'][sEpoch1ID:eEpoch1ID],label='ACE')
plt.plot(corEpoch[sEpoch2ID:eEpoch2ID],aceDataMD['N'][sEpoch2ID:eEpoch2ID],label='propagated ACE')
plt.legend(loc='best')
plt.ylabel('Solar Wind Density (N/cc)')
plt.subplot(4,1,3)
plt.plot(imp8DataMD['epoch'][sEpoch1ID:eEpoch1ID],imp8DataMD['Bz'][sEpoch1ID:eEpoch1ID],label='IMP8')
plt.plot( aceDataMD['epoch'][sEpoch1ID:eEpoch1ID], aceDataMD['Bz'][sEpoch1ID:eEpoch1ID],label='ACE')
plt.plot(corEpoch[sEpoch2ID:eEpoch2ID],aceDataMD['Bz'][sEpoch2ID:eEpoch2ID],label='propagated ACE')
plt.legend(loc='best')
plt.ylabel('$B_z$ (nT)')
plt.subplot(4,1,4)
plt.plot( aceDataMD['epoch'][sEpoch1ID:eEpoch1ID], aceDataMD['Vx'][sEpoch1ID:eEpoch1ID],label='ACE')
plt.plot(corEpoch[sEpoch2ID:eEpoch2ID],aceDataMD['Vx'][sEpoch2ID:eEpoch2ID],label='propagated ACE')
plt.legend(loc='best')
plt.ylabel('$V_x$ (Km/s)')
plt.xlabel('Epoch')
plt.suptitle('1999-01-15')

plt.figure(400)
plt.subplot(3,1,1)
plt.plot(abs(ccorr(imp8DataMD['V'][sEpoch1ID:eEpoch1ID],aceDataMD['V'][sEpoch2ID:eEpoch2ID])))
plt.ylabel('Velocity Correlation')
plt.title('Correlation between Solar Wind Parameters at ACE (after lagging) and IMP8')
plt.subplot(3,1,2)
plt.plot(abs(ccorr(imp8DataMD['N'][sEpoch1ID:eEpoch1ID],aceDataMD['N'][sEpoch2ID:eEpoch2ID])))
plt.ylabel('Density Correlation')
plt.subplot(3,1,3)
plt.plot(abs(ccorr(imp8DataMD['Bz'][sEpoch1ID:eEpoch1ID],aceDataMD['Bz'][sEpoch2ID:eEpoch2ID])))
plt.ylabel('Bz Correlation')
plt.xlabel('Index')

plt.figure(500)
plt.plot(corEpoch,lagging)
plt.title('Time Lag between ACE and IMP8 as a function of Epoch time')
plt.xlabel('Epoch')
plt.ylabel('Time Lag (seconds)')

aceKDE = gaussian_kde(aceDataMD['V'][sEpoch2ID:eEpoch2ID], bw_method=kdeBW)
imp8KDE = gaussian_kde(imp8DataMD['V'][sEpoch1ID:eEpoch1ID], bw_method=kdeBW)
v_eval = numpy.linspace(500, 610, num=110)

plt.figure(600)
plt.plot(v_eval, aceKDE(v_eval), 'b-', label="ACE-KDE")
plt.plot(v_eval, imp8KDE(v_eval), 'r-', label="IMP8-KDE")
plt.title('Kernel Density Estimation for Solar Wind Velocity at ACE and IMP8')
plt.xlabel('Solar Wind Velocity (Km/s)')
plt.ylabel('Kernel Density Estimation')
plt.legend(loc='best')

plt.show()
sys.exit()

RE = 6371.0
plt.figure(100)
plt.subplot(3,1,1)
plt.plot(imp8Data['magneticEpoch'],imp8Data['SCxGSE']/RE)
plt.plot(aceData['magneticEpoch'],aceData['SCxGSE']/RE)
plt.subplot(3,1,2)
plt.plot(imp8Data['magneticEpoch'],imp8Data['SCyGSE']/RE)
plt.plot(aceData['magneticEpoch'],aceData['SCyGSE']/RE)
plt.subplot(3,1,3)
plt.plot(imp8Data['magneticEpoch'],imp8Data['SCzGSE']/RE)
plt.plot(aceData['magneticEpoch'],aceData['SCzGSE']/RE)

imp8DataMD = imp8Data.copy()
imp8DataMD['V'] = swMedFilter(imp8Data['plasmaEpoch'],imp8Data['V'],45*60)
imp8DataMD['N'] = swMedFilter(imp8Data['plasmaEpoch'],imp8Data['N'],45*60)
imp8DataMD['Bz'] = swMedFilter(imp8Data['magneticEpoch'],imp8Data['Bz'],45*60)

aceDataMD = aceData.copy()
aceDataMD['V'] = swMedFilter(aceData['plasmaEpoch'],aceData['V'],45*60)
aceDataMD['Vx'] = swMedFilter(aceData['plasmaEpoch'],aceData['Vx'],45*60)
aceDataMD['Vy'] = swMedFilter(aceData['plasmaEpoch'],aceData['Vy'],45*60)
aceDataMD['Vz'] = swMedFilter(aceData['plasmaEpoch'],aceData['Vz'],45*60)
aceDataMD['N'] = swMedFilter(aceData['plasmaEpoch'],aceData['N'],45*60)
aceDataMD['Bz'] = swMedFilter(aceData['magneticEpoch'],aceData['Bz'],45*60)

plt.figure(200)
plt.subplot(3,1,1)
plt.plot(imp8DataMD['plasmaEpoch'],imp8DataMD['V'],label='IMP8')
plt.plot(aceDataMD['plasmaEpoch'],aceDataMD['V'],label='ACE')
plt.subplot(3,1,2)
plt.plot(imp8DataMD['plasmaEpoch'],imp8DataMD['N'],label='IMP8')
plt.plot( aceDataMD['plasmaEpoch'], aceDataMD['N'],label='ACE')
plt.subplot(3,1,3)
plt.plot(imp8DataMD['magneticEpoch'],imp8DataMD['Bz'],label='IMP8')
plt.plot( aceDataMD['magneticEpoch'], aceDataMD['Bz'],label='ACE')

iCounter=0
RE = 6371.0
plt.figure(200)
for i in range(0,len(imp8Data['SCxGSE']),3000):
 iCounter = iCounter + 1
 if iCounter > 1: break
 if imp8Data['SCxGSE'][i] >= 0:
  destPos = {'X':imp8Data['SCxGSE'][i],'Y':imp8Data['SCyGSE'][i],'Z':imp8Data['SCzGSE'][i]}
 #destPos = {'X':imp8Data['SCxGSE'][1],'Y':imp8Data['SCyGSE'][-1],'Z':imp8Data['SCzGSE'][i]}
  lagging, corEpoch = getTimeLag(aceDataMD,destPos,method='standard')
  lagging = dataClean(lagging,[1e-12],['<='])
  plt.subplot(3,1,1)
  plt.plot(corEpoch,aceDataMD['V'],label='propagated ACE')
  plt.legend(loc='best')
  plt.ylabel('Solar Wind Speed (km/s)')
  plt.title('Solar Wind - Time Lag')
  plt.subplot(3,1,2)
  plt.plot(corEpoch,aceDataMD['N'],label='propagated ACE')
  plt.legend(loc='best')
  plt.ylabel('Solar Wind Density (N/cc)')
  plt.xlabel('Epoch')
# plt.plot(corcEpoch,aceData['V'],label='$IMP8_x$ = ' + str(int(imp8Data['SCxGSE'][i]/RE)) + ' $R_E$')
# plt.plot(aceData['plasmaEpoch'],numpy.array(lagging)/60.0,label='$IMP8_x$ = ' + str(int(imp8Data['SCxGSE'][i]/RE))) + ' $R_E$'
plt.suptitle('1999-01-15')


plt.show()

sys.exit()


