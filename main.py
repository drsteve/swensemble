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
from scipy.stats import pearsonr

from swdatanal import getDistrib, getIndices, omniDataCorr, ccorr
from swdatanal import getSolarWindType, getTimeLag, swMedFilter
from getswdata import getOMNIfiles, getOMNIdata, getOMNIparams
from getswdata import getACEfiles, getACEdata, getACEparams, aceDataAdjust
from getswdata import getIMP8files,getIMP8data, getIMP8params, imp8DataAdjust
from getswdata import getWINDfiles,getWINDdata, getWINDparams
from getswdata import getGeotailfiles,getGeotaildata, getGeotailparams
from getswdata import dataClean, dateShift, dateList, epochMatch, removeNaN


startDate   = datetime.datetime(1999, 1,15,13, 0, 0)
endDate     = datetime.datetime(1999, 1,15,21,59,59)
locDateList = dateList(startDate, endDate, shift = 'hour')

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
  imp8Data['N']   = numpy.array(imp8Data['proton_density_fit'])
  imp8Data['V']   = numpy.array(imp8Data['V_fit'])
  imp8Data['T']   = numpy.array(imp8Data['protonV_thermal_fit'])
  imp8Data['Bx']  = numpy.array(imp8Data['B_Vector_GSE'][0])
  imp8Data['By']  = numpy.array(imp8Data['B_Vector_GSE'][1])
  imp8Data['Bz']  = numpy.array(imp8Data['B_Vector_GSE'][2])
  BB = []
  for i in range(len(imp8Data['Bx'])):
   BB.extend([math.sqrt(imp8Data['Bx'][i]**2+imp8Data['By'][i]**2+imp8Data['Bz'][i]**2)])
  imp8Data['B']   = BB
  imp8Data['SCxGSE'] = numpy.array(imp8Data['SC_Pos_GSE'][0])*RE
  imp8Data['SCyGSE'] = numpy.array(imp8Data['SC_Pos_GSE'][1])*RE
  imp8Data['SCzGSE'] = numpy.array(imp8Data['SC_Pos_GSE'][2])*RE

dataSRC = 'ace'
if dataSRC == 'ace':
 aceParams = getACEparams(locDateList,'/home/ehab/SWData',aceSet='16sec',aceDat='magnetic')
 aceParams = getACEparams(locDateList,'/home/ehab/SWData',aceSet='1min',aceDat='plasma')
 if aceParams != '--empty--':
  swParams = ['Magnitude','Np','Vp','Tpr','BGSEc','V_GSE','SC_pos_GSE']
  aceData  = getACEdata(locDateList,'/home/ehab/SWData',swParams,['16sec','1min'],dataStat='clean')
  aceData['B']   = numpy.array(aceData['Magnitude'])
  aceData['N']   = numpy.array(aceData['Np'])
  aceData['V']   = numpy.array(aceData['Vp'])
  aceData['T']   = numpy.array(aceData['Tpr'])
  aceData['Bx']  = numpy.array(aceData['BGSEc'][0])
  aceData['By']  = numpy.array(aceData['BGSEc'][1])
  aceData['Bz']  = numpy.array(aceData['BGSEc'][2])
  aceData['Vx']  = numpy.array(aceData['V_GSE'][0])
  aceData['Vy']  = numpy.array(aceData['V_GSE'][1])
  aceData['Vz']  = numpy.array(aceData['V_GSE'][2])
  aceData['SCxGSE'] = numpy.array(aceData['SC_pos_GSE'][0])
  aceData['SCyGSE'] = numpy.array(aceData['SC_pos_GSE'][1])
  aceData['SCzGSE'] = numpy.array(aceData['SC_pos_GSE'][2])

  if swCats:
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

aceDataAdjust(aceData,aceData['plasmaEpoch'],epochLen='fixed')
imp8DataAdjust(imp8Data,imp8Data['plasmaEpoch'],epochLen='fixed')
sys.exit()

TT,VV = removeNaN(aceData['plasmaEpoch'],aceData['V'])
aceVmod=epochMatch(imp8Data['plasmaEpoch'],TT,VV)
TT,BB = removeNaN(aceData['magneticEpoch'],aceData['Bz'])
aceBzmod=epochMatch(imp8Data['plasmaEpoch'],TT,BB)
plt.subplot(3,1,1)
plt.plot(imp8Data['plasmaEpoch'],aceBzmod)
plt.subplot(3,1,2)
plt.plot(imp8Data['plasmaEpoch'],aceVmod)
plt.subplot(3,1,3)
plt.plot(aceData['magneticEpoch'],aceData['Bz'])
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


