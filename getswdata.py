'''Module to get data sets from CDAweb and SPDF
'''
import os, glob, httplib, bisect
import datetime as dt
import numpy as np
from spacepy import pycdf

def checkConnection():
    conn = httplib.HTTPConnection("www.google.com")
    try:
        conn.request("HEAD", "/")
        return True
    except:
        return False


def getOMNIfiles(omniDates,dataLoc,omniSet='hourly'):
    omniSet = omniSet.lower()

    if dataLoc[len(dataLoc)-1] == "/":
        dataLoc = dataLoc[:len(dataLoc)]
    else:
        dataLoc = dataLoc[:len(dataLoc)] + "/"

    fnames = []

    for i in range(len(omniDates)):
        if i > 0:
            if omniDates[i] == omniDates[i-1]:
                continue
            elif omniDates[i].year == omniDates[i-1].year:
                if omniDates[i].month == omniDates[i-1].month:
                    continue
                elif omniSet == 'hourly':
                    if omniDates[i].month <= 6 and omniDates[i-1].month <= 6:
                        continue
                    elif omniDates[i].month >= 7 and omniDates[i-1].month >= 7:
                        continue

        omniYear = str(omniDates[i].year)
        if omniSet == 'hourly' and omniDates[i].month <= 6:
            omniDate = str(omniDates[i].year) + "0101"
        elif omniSet == 'hourly' and omniDates[i].month >= 7:
            omniDate = str(omniDates[i].year) + "0701"
        else:
            omniDate = str(omniDates[i].year) + "%02d" % (omniDates[i].month) + "01"

        if omniSet=='hourly' and int(omniYear) >= 1963:
             fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hourly/' + str(omniYear) + '/omni2_h0_mrg1hr_' + str(omniDate) + '*.cdf'
        elif omniSet=='5min' and int(omniYear) >= 1981:
             fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_5min/' + str(omniYear) + '/omni_hro_5min_' + str(omniDate) + '*.cdf'
        elif omniSet=='1min' and int(omniYear) >= 1981:
             fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/' + str(omniYear) + '/omni_hro_1min_' + str(omniDate) + '*.cdf'
        else:
             fullPath = ""
             print(omniSet + '-' + omniYear + ' OMNI data is not available. Only 1hour (1963:), 5min (1981:), 1min (1981:).')
             return ""

        if sorted(glob.glob(dataLoc + fullPath)) == []:
            print("OMNI data could not be found locally, try to download OMNI data from destination ....")
            if checkConnection():
                webPath = "ftp://" + fullPath[:len(fullPath)]
                print("Connect to " + webPath)
                envComm = "wget -P /home/ehab/SWData" + " -r -l1 -A.cdf " + webPath
                print(envComm)
                os.system(envComm)
            else:
                print("No internet connection found ... Check you network settings.")

        localPath = dataLoc + fullPath
        fnames.extend(glob.glob(localPath))

    return sorted(fnames)

def getOMNIparams(omniDates, dataLoc, omniSet='hourly'):
    OMNIfnames  = getOMNIfiles(omniDates, dataLoc, omniSet)
    if len(OMNIfnames) > 0:
        cdfKeys = pycdf.CDF(OMNIfnames[0])
    else:
        cdfKeys = "--empty--"
    return cdfKeys


def getOMNIdata(omniDates,dataLoc,dataList,omniSet='hourly',dataStat='raw'):
    OMNIfnames = getOMNIfiles(omniDates,dataLoc,omniSet)
    OMNIparams = getOMNIparams(omniDates,dataLoc,omniSet)

    dataStorage = []
    for j in range(len(dataList)+1):
     dataStorage.append([])

    aFlage = True
    for j in range(len(dataList)):
        if dataList[j] in OMNIparams:
             for i in range(len(OMNIfnames)):
                 cdf = pycdf.CDF(OMNIfnames[i])
                 dataStorage[j].extend(list(cdf[dataList[j]]))
                 if aFlage: dataStorage[len(dataList)].extend(list(cdf['Epoch']))
             if aFlage: aFlage = False
             cdf.close()
        else:
          print dataList[j], ' is not available in OMNI CDF File!'

    swData = {'Epoch':dataStorage[len(dataList)]}
    for j in range(len(dataList)):
     if dataStat == 'clean':
      if dataList[j] == 'V': dataStorage[j] = dataClean(dataStorage[j], [2500], ['>='])
      if dataList[j] == 'N': dataStorage[j] = dataClean(dataStorage[j], [999.0, 0.0], ['>=', '<='])
      if dataList[j] == 'T': dataStorage[j] = dataClean(dataStorage[j], [1.0e7], ['>='])
      if dataList[j] == 'ABS_B': dataStorage[j] = dataClean(dataStorage[j], [0.0], ['<='])
      if dataList[j] == 'BX_GSE': dataStorage[j] = dataClean(dataStorage[j], [0.0], ['<='])
      if dataList[j] == 'BY_GSE': dataStorage[j] = dataClean(dataStorage[j], [0.0], ['<='])
      if dataList[j] == 'BZ_GSE': dataStorage[j] = dataClean(dataStorage[j], [0.0], ['<='])
     swData[dataList[j]] = dataStorage[j]

    return swData


def omniDataAdjust(Data, epoch=[], epochLen='not fixed'):
    BB = []
    epochDiff = max(np.diff(Data['Epoch']))
    if epochDiff.seconds >= 3600:
        omniSet = 'hourly'
    else:
        omniSet = 'minutes'
    if epochLen == 'fixed' and epoch != []:
        DataAdj={'epoch': epoch}
        if omniSet == 'hourly':
            DataAdj['B']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['ABS_B'], interpKind='linear'))
            DataAdj['Bx']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['BX_GSE'], interpKind='linear'))
            DataAdj['By']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['BY_GSE'], interpKind='linear'))
            DataAdj['Bz']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['BZ_GSE'], interpKind='linear'))
            DataAdj['N']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['N'], interpKind='linear'))
            DataAdj['V']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['V'], interpKind='linear'))
            DataAdj['T']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['T'], interpKind='linear'))
        else:
            DataAdj['Bx']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['BX_GSE'], interpKind='linear'))
            DataAdj['By']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['BY_GSE'], interpKind='linear'))
            DataAdj['Bz']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['BZ_GSE'], interpKind='linear'))
            for i in range(len(DataAdj['Bx'])):
                BB.extend([np.sqrt(DataAdj['Bx'][i]**2 + DataAdj['By'][i]**2 + DataAdj['Bz'][i]**2)])
            DataAdj['B']      = np.array(BB)
            DataAdj['N']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['proton_density'], interpKind='linear'))
            DataAdj['V']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['flow_speed'], interpKind='linear'))
            DataAdj['T']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['T'], interpKind='linear'))
            DataAdj['Vx']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['Vx'], interpKind='linear'))
            DataAdj['Vy']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['Vy'], interpKind='linear'))
            DataAdj['Vz']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['Vz'], interpKind='linear'))
    else:
         DataAdj={'epoch':Data['Epoch']}
         if omniSet == 'hourly':
            DataAdj['B']      = np.array(Data['ABS_B'])
            DataAdj['Bx']     = np.array(Data['BX_GSE'])
            DataAdj['By']     = np.array(Data['BY_GSE'])
            DataAdj['Bz']     = np.array(Data['BZ_GSE'])
            DataAdj['N']      = np.array(Data['N'])
            DataAdj['V']      = np.array(Data['V'])
            DataAdj['T']      = np.array(Data['T'])
         else:
            DataAdj['Bx']     = np.array(Data['BX_GSE'])
            DataAdj['By']     = np.array(Data['BY_GSE'])
            DataAdj['Bz']     = np.array(Data['BZ_GSE'])
            for i in range(len(DataAdj['Bx'])):
                BB.extend([np.sqrt(DataAdj['Bx'][i]**2 + DataAdj['By'][i]**2 + DataAdj['Bz'][i]**2)])
            DataAdj['B']      = np.array(BB)
            DataAdj['N']      = np.array(Data['proton_density'])
            DataAdj['V']      = np.array(Data['flow_speed'])
            DataAdj['T']      = np.array(Data['T'])
            DataAdj['Vx']     = np.array(Data['Vx'])
            DataAdj['Vy']     = np.array(Data['Vy'])
            DataAdj['Vz']     = np.array(Data['Vz'])

    return DataAdj


def getACEfiles(aceDates, dataLoc, aceSet='hourly', aceDat='magnetic'):
    aceSet = aceSet.lower()
    aceDat = aceDat.lower()

    if dataLoc[len(dataLoc)-1] == "/":
        dataLoc = dataLoc[:len(dataLoc)]
    else:
        dataLoc = dataLoc[:len(dataLoc)] + "/"

    fnames = []

    for i in range(len(aceDates)):
        if i > 0:
            if aceDates[i] == aceDates[i-1]:
                continue
            elif aceDates[i].year == aceDates[i-1].year:
                if aceDates[i].month == aceDates[i-1].month:
                    if aceDates[i].day == aceDates[i-1].day:
                        continue

        aceYear = str(aceDates[i].year)
        aceDate = str(aceDates[i].year) + "%02d" % (aceDates[i].month) + "%02d" % (aceDates[i].day)
        if aceDat == 'plasma':
            if aceSet == 'hourly' and int(aceYear) >= 1998:
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/ace/swepam/level_2_cdaweb/swe_h2/' + str(aceYear) + '/ac_h2_swe_' + str(aceDate) + '*.cdf'
            elif aceSet == '1min' and int(aceYear) >= 1998:
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/ace/swepam/level_2_cdaweb/swe_h0/' + str(aceYear) + '/ac_h0_swe_' + str(aceDate) + '*.cdf'
            else:
                fullPath = ""
                print(aceSet + ' ACE Plasma data is not available. Only 1hour (1998:), 23min (1998:), 5min (2011:).')
                return fullPath
        
        if aceDat == 'magnetic':
            if aceSet == 'hourly' and int(aceYear) >= 1998:
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/ace/mag/level_2_cdaweb/mfi_h2/' + str(aceYear) + '/ac_h2_mfi_' + str(aceDate) + '*.cdf'
            elif aceSet == '4min' and int(aceYear) >= 1998:
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/ace/mag/level_2_cdaweb/mfi_h1/' + str(aceYear) + '/ac_h1_mfi_' + str(aceDate) + '*.cdf'
            elif aceSet == '16sec' and int(aceYear) >= 1998:
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/ace/mag/level_2_cdaweb/mfi_h0/' + str(aceYear) + '/ac_h0_mfi_' + str(aceDate) + '*.cdf'
            elif aceSet == '1sec' and int(aceYear) >= 1998:
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/ace/mag/level_2_cdaweb/mfi_h3/' + str(aceYear) + '/ac_h3_mfi_' + str(aceDate) + '*.cdf'
            else:
                fullPath = ""
                print(aceSet + '-' + aceYear + \
                    ' ACE Magnetic data is not available. Only 1hour (1998:), 5min (2011:), 4min (1998:), 16sec (1998:), 1sec (1998:).')
                return fullPath 

        if sorted(glob.glob(dataLoc + fullPath)) == []:
            print("ACE data could not be found locally, try to download ACE data from destination ....")
            if checkConnection():
                webPath = "ftp://" + fullPath[:len(fullPath)]
                print("Connect to " + webPath)
                envComm = "wget -P " + dataLoc + " -r -l1 -A.cdf " + webPath
                print(envComm)
                os.system(envComm)
            else:
                print("No internet connection found ... Check you network settings.")

        localPath = dataLoc + fullPath
        fnames.extend(glob.glob(localPath))

    return sorted(fnames)


def getACEparams(aceDates,dataLoc,aceSet='hourly',aceDat='magnetic'):
    aceSet = aceSet.lower()
    aceDat = aceDat.lower()

    if aceDat == 'magnetic':
        ACEfnames  = getACEfiles(aceDates,dataLoc,aceSet,aceDat='magnetic')
    elif aceDat == 'plasma':
        ACEfnames  = getACEfiles(aceDates,dataLoc,aceSet,aceDat='plasma')
    if len(ACEfnames) > 0:
        cdfKeys = pycdf.CDF(ACEfnames[0])
    else:
        cdfKeys = "--empty--"
    return cdfKeys


def getACEdata(aceDates,dataLoc,dataList,aceSet=['1sec','1min'],dataStat='raw'):
    from numpy import transpose, array

    dataStat = dataStat.lower()

    if len(aceSet) == 2:
        aceSet[0] = aceSet[0].lower()
        aceSet[1] = aceSet[1].lower()
        ACEBparams = getACEparams(aceDates,dataLoc,aceSet[0],aceDat='magnetic')
        ACEPparams = getACEparams(aceDates,dataLoc,aceSet[1],aceDat='plasma')
    elif len(aceSet) == 1:
        aceSet[0] = aceSet[0].lower()
        ACEBparams = getACEparams(aceDates,dataLoc,aceSet[0],aceDat='magnetic')
        ACEPparams = getACEparams(aceDates,dataLoc,aceSet='1min',aceDat='plasma')

    for j in range(len(dataList)):
        if dataList[j] in ACEBparams: bFlag = True
        if dataList[j] in ACEPparams: pFlag = True

    dataStorage = []

    if bFlag and pFlag:
        bpFlag = True; bFlag  = False; pFlag  = False
        ACEBfnames = getACEfiles(aceDates, dataLoc, aceSet[0], aceDat='magnetic')
        if len(aceSet) == 2:
            ACEPfnames = getACEfiles(aceDates, dataLoc, aceSet[1], aceDat='plasma')
        elif len(aceSet) == 1:
            ACEPfnames = getACEfiles(aceDates, dataLoc, aceSet='1min', aceDat='plasma')

        for j in range(len(dataList)+2):
            dataStorage.append([])

        for i in range(len(ACEBfnames)):
            cdf = pycdf.CDF(ACEBfnames[i])
            dataStorage[len(dataList)].extend(list(cdf['Epoch']))
            for j in range(len(dataList)):
                if dataList[j] in ACEBparams:
                    dataStorage[j].extend(list(cdf[dataList[j]]))
            cdf.close()

        for i in range(len(ACEPfnames)):
            cdf = pycdf.CDF(ACEPfnames[i])
            dataStorage[len(dataList)+1].extend(list(cdf['Epoch']))
            for j in range(len(dataList)):
                if dataList[j] in ACEPparams:
                    dataStorage[j].extend(list(cdf[dataList[j]]))
            cdf.close()

        sEpochBID  = bisect.bisect_left(dataStorage[len(dataList)], aceDates[0])
        eEpochBID  = bisect.bisect_left(dataStorage[len(dataList)], aceDates[len(aceDates)-1])
        sEpochPID  = bisect.bisect_left(dataStorage[len(dataList)+1], aceDates[0])
        eEpochPID  = bisect.bisect_left(dataStorage[len(dataList)+1], aceDates[len(aceDates)-1])
        swData = {'magneticEpoch':dataStorage[len(dataList)][sEpochBID:eEpochBID], 'plasmaEpoch':dataStorage[len(dataList)+1][sEpochPID:eEpochPID]}

    elif bFlag:
        ACEBfnames = getACEfiles(aceDates, dataLoc, aceSet[0], aceDat='magnetic')

        for j in range(len(dataList)+1):
            dataStorage.append([])

        for i in range(len(ACEBfnames)):
            cdf = pycdf.CDF(ACEBfnames[i])
            dataStorage[len(dataList)].extend(list(cdf['Epoch']))
            for j in range(len(dataList)):
                if dataList[j] in ACEBparams:
                    dataStorage[j].extend(list(cdf[dataList[j]]))
            cdf.close()
        sEpochBID  = bisect.bisect_left(dataStorage[len(dataList)], aceDates[0])
        eEpochBID  = bisect.bisect_left(dataStorage[len(dataList)], aceDates[len(aceDates)-1])
        swData = {'magneticEpoch':dataStorage[len(dataList)][sEpochBID:eEpochBID]}

    elif pFlag:
        if len(aceSet) == 2:
            ACEPfnames = getACEfiles(aceDates, dataLoc, aceSet[1], aceDat='plasma')
        elif len(aceSet) == 1:
            ACEPfnames = getACEfiles(aceDates, dataLoc, aceSet='1min', aceDat='plasma')

        for j in range(len(dataList)+1):
            dataStorage.append([])

        for i in range(len(ACEPfnames)):
            cdf = pycdf.CDF(ACEPfnames[i])
            dataStorage[len(dataList)].extend(list(cdf['Epoch']))
            for j in range(len(dataList)):
                if dataList[j] in ACEPparams:
                    dataStorage[j].extend(list(cdf[dataList[j]]))
            cdf.close()
        sEpochPID  = bisect.bisect_left(dataStorage[len(dataList)], aceDates[0])
        eEpochPID  = bisect.bisect_left(dataStorage[len(dataList)], aceDates[len(aceDates)-1])
        swData = {'plasmaEpoch':dataStorage[len(dataList)][sEpochPID:eEpochPID]}

    for j in range(len(dataList)):
        if dataList[j] in ['BGSEc','V_GSE','SC_pos_GSE']:
            swDataTMP = transpose(dataStorage[j])
            if dataList[j] in ACEBparams:
                swData[dataList[j]] = [swDataTMP[0][sEpochBID:eEpochBID], swDataTMP[1][sEpochBID:eEpochBID], swDataTMP[2][sEpochBID:eEpochBID]]
                if dataStat == 'clean':
                    if dataList[j] == 'BGSEc':
                        swData[dataList[j]][0] = dataClean(swData[dataList[j]][0], [999, -999], ['>', '<'])
                        swData[dataList[j]][1] = dataClean(swData[dataList[j]][1], [999, -999], ['>', '<'])
                        swData[dataList[j]][2] = dataClean(swData[dataList[j]][2], [999, -999], ['>', '<'])
            elif dataList[j] in ACEPparams:
                swData[dataList[j]] = [swDataTMP[0][sEpochPID:eEpochPID], swDataTMP[1][sEpochPID:eEpochPID], swDataTMP[2][sEpochPID:eEpochPID]]
                if dataStat == 'clean':
                    if dataList[j] == 'V_GSE':
                        swData[dataList[j]][0] = dataClean(abs(swData[dataList[j]][0]), [2500, 100], ['>', '<'])
                        swData[dataList[j]][1] = dataClean(abs(swData[dataList[j]][1]), [2500, 100], ['>', '<'])
                        swData[dataList[j]][2] = dataClean(abs(swData[dataList[j]][2]), [2500, 100], ['>', '<'])
        else:
            if dataList[j] in ACEBparams:
                swData[dataList[j]] = dataStorage[j][sEpochBID:eEpochBID]
            elif dataList[j] in ACEPparams:
                swData[dataList[j]] = dataStorage[j][sEpochPID:eEpochPID]
            if dataStat == 'clean':
                if dataList[j] == 'Np': swData[dataList[j]] = dataClean(swData[dataList[j]], [999, 0], ['>=', '<='])
                if dataList[j] == 'Vp': swData[dataList[j]] = dataClean(swData[dataList[j]], [2500, 100], ['>', '<'])
                if dataList[j] == 'Magnitude': swData[dataList[j]] = dataClean(swData[dataList[j]], [-999, 0, 999], ['<=', '=', '>='])
                if dataList[j] == 'Tpr': swData[dataList[j]] = dataClean(swData[dataList[j]], [1e8, 0], ['>=', '<='])

    return swData


'''
    if len(aceSet) == 2:
     aceSet[0] = aceSet[0].lower()
     aceSet[1] = aceSet[1].lower()
     ACEBfnames = getACEfiles(aceDates,dataLoc,aceSet[0],aceDat='magnetic')
     ACEBparams = getACEparams(aceDates,dataLoc,aceSet[0],aceDat='magnetic')
     ACEPfnames = getACEfiles(aceDates,dataLoc,aceSet[1],aceDat='plasma')
     ACEPparams = getACEparams(aceDates,dataLoc,aceSet[1],aceDat='plasma')
    elif len(aceSet) == 1:
     aceSet[0] = aceSet[0].lower()
     ACEBfnames = getACEfiles(aceDates,dataLoc,aceSet[0],aceDat='magnetic')
     ACEBparams = getACEparams(aceDates,dataLoc,aceSet[0],aceDat='magnetic')
     ACEPfnames = getACEfiles(aceDates,dataLoc,aceSet='1min',aceDat='plasma')
     ACEPparams = getACEparams(aceDates,dataLoc,aceSet='1min',aceDat='plasma')
    dataStorage = []
    for j in range(len(dataList)+2):
     dataStorage.append([])

    for j in range(len(dataList)):
     if dataList[j] in ACEBparams: bpFlag = True
     if dataList[j] in ACEPparams: ppFlag = True

    bFlage = True
    pFlage = True
    for j in range(len(dataList)):
     if (dataList[j] in ACEBparams):
      for i in range(len(ACEBfnames)):
       cdf = pycdf.CDF(ACEBfnames[i])
       dataStorage[j].extend(list(cdf[dataList[j]]))
       if bFlage: dataStorage[len(dataList)].extend(list(cdf['Epoch']))
      if bFlage: bFlage = False
      cdf.close()
     elif (dataList[j] in ACEPparams):
      for i in range(len(ACEPfnames)):
       cdf = pycdf.CDF(ACEPfnames[i])
       dataStorage[j].extend(list(cdf[dataList[j]]))
       if pFlage: dataStorage[len(dataList)+1].extend(list(cdf['Epoch']))
      if pFlage: pFlage = False
      cdf.close()
     else:
       print dataList[j], ' is not available in ACE CDF File!'

    sEpochBID  = bisect_left(dataStorage[len(dataList)], aceDates[0])
    eEpochBID  = bisect_left(dataStorage[len(dataList)], aceDates[len(aceDates)-1])
    sEpochPID  = bisect_left(dataStorage[len(dataList)+1], aceDates[0])
    eEpochPID  = bisect_left(dataStorage[len(dataList)+1], aceDates[len(aceDates)-1])

    swData = {'magneticEpoch':dataStorage[len(dataList)][sEpochBID:eEpochBID], 'plasmaEpoch':dataStorage[len(dataList)+1][sEpochPID:eEpochPID]}
'''

def aceDataAdjust(Data,epoch=[],epochLen='not fixed'):
    if epochLen == 'fixed' and epoch == Data['plasmaEpoch']:
        DataAdj={'epoch':epoch}
        DataAdj['B']      = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['Magnitude'],interpKind='linear'))
        DataAdj['Bx']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['BGSEc'][0],interpKind='linear'))
        DataAdj['By']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['BGSEc'][1],interpKind='linear'))
        DataAdj['Bz']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['BGSEc'][2],interpKind='linear'))
        DataAdj['SCxGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_pos_GSE'][0],interpKind='linear'))
        DataAdj['SCyGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_pos_GSE'][1],interpKind='linear'))
        DataAdj['SCzGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_pos_GSE'][2],interpKind='linear'))
        DataAdj['N']      = np.array(Data['Np'])
        DataAdj['V']      = np.array(Data['Vp'])
        DataAdj['T']      = np.array(Data['Tpr'])
        DataAdj['Vx']     = np.array(Data['V_GSE'][0])
        DataAdj['Vy']     = np.array(Data['V_GSE'][1])
        DataAdj['Vz']     = np.array(Data['V_GSE'][2])
    elif epochLen == 'fixed' and epoch == Data['magneticEpoch']:
        DataAdj={'epoch':epoch}
        DataAdj['B']      = np.array(Data['Magnitude'])
        DataAdj['Bx']     = np.array(Data['BGSEc'][0])
        DataAdj['By']     = np.array(Data['BGSEc'][1])
        DataAdj['Bz']     = np.array(Data['BGSEc'][2])
        DataAdj['SCxGSE'] = np.array(Data['SC_pos_GSE'][0])
        DataAdj['SCyGSE'] = np.array(Data['SC_pos_GSE'][1])
        DataAdj['SCzGSE'] = np.array(Data['SC_pos_GSE'][2])
        DataAdj['N']      = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['Np'],interpKind='linear'))
        DataAdj['V']      = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['Vp'],interpKind='linear'))
        DataAdj['T']      = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['Tpr'],interpKind='linear'))
        DataAdj['Vx']     = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['V_GSE'][0],interpKind='linear'))
        DataAdj['Vy']     = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['V_GSE'][1],interpKind='linear'))
        DataAdj['Vz']     = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['V_GSE'][2],interpKind='linear'))
    elif epochLen == 'fixed' and epoch != []:
        DataAdj={'epoch':epoch}
        DataAdj['B']      = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['Magnitude'],interpKind='linear'))
        DataAdj['Bx']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['BGSEc'][0],interpKind='linear'))
        DataAdj['By']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['BGSEc'][1],interpKind='linear'))
        DataAdj['Bz']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['BGSEc'][2],interpKind='linear'))
        DataAdj['SCxGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_pos_GSE'][0],interpKind='linear'))
        DataAdj['SCyGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_pos_GSE'][1],interpKind='linear'))
        DataAdj['SCzGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_pos_GSE'][2],interpKind='linear'))
        DataAdj['N']      = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['Np'],interpKind='linear'))
        DataAdj['V']      = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['Vp'],interpKind='linear'))
        DataAdj['T']      = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['Tpr'],interpKind='linear'))
        DataAdj['Vx']     = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['V_GSE'][0],interpKind='linear'))
        DataAdj['Vy']     = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['V_GSE'][1],interpKind='linear'))
        DataAdj['Vz']     = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['V_GSE'][2],interpKind='linear'))
    else:
        DataAdj={'plasmaEpoch':Data['plasmaEpoch'],'magneticEpoch':Data['magneticEpoch']}
        DataAdj['B']      = np.array(Data['Magnitude'])
        DataAdj['Bx']     = np.array(Data['BGSEc'][0])
        DataAdj['By']     = np.array(Data['BGSEc'][1])
        DataAdj['Bz']     = np.array(Data['BGSEc'][2])
        DataAdj['SCxGSE'] = np.array(Data['SC_pos_GSE'][0])
        DataAdj['SCyGSE'] = np.array(Data['SC_pos_GSE'][1])
        DataAdj['SCzGSE'] = np.array(Data['SC_pos_GSE'][2])
        DataAdj['N']      = np.array(Data['Np'])
        DataAdj['V']      = np.array(Data['Vp'])
        DataAdj['T']      = np.array(Data['Tpr'])
        DataAdj['Vx']     = np.array(Data['V_GSE'][0])
        DataAdj['Vy']     = np.array(Data['V_GSE'][1])
        DataAdj['Vz']     = np.array(Data['V_GSE'][2])

    return DataAdj


def getIMP8files(imp8Dates,dataLoc,imp8Set='15sec',imp8Dat='magnetic'):
    import os
    import glob

    imp8Set = imp8Set.lower()
    imp8Dat = imp8Dat.lower()

    if dataLoc[len(dataLoc)-1] == "/":
        dataLoc = dataLoc[:len(dataLoc)]
    else:
        dataLoc = dataLoc[:len(dataLoc)] + "/"

    fnames = []

    for i in range(len(imp8Dates)):
        if i > 0:
            if imp8Dates[i] == imp8Dates[i-1]:
                continue
            elif imp8Dates[i].year == imp8Dates[i-1].year:
                if imp8Dates[i].month == imp8Dates[i-1].month:
                    if imp8Dates[i].day == imp8Dates[i-1].day:
                        continue

        imp8Year = str(imp8Dates[i].year)
        imp8Date = str(imp8Dates[i].year) + "%02d" % (imp8Dates[i].month) + "%02d" % (imp8Dates[i].day)

        if imp8Dat == 'magnetic':
            if imp8Set == '15sec' and (int(imp8Year) >= 1973 and int(imp8Year) <= 2000):
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/imp/imp8/mag/mag_15sec_cdaweb/' + str(imp8Year) + '/i8_15sec_mag_' + str(imp8Date) + '*.cdf'
            elif imp8Set == '320msec' and (int(imp8Year) >= 1973 and int(imp8Year) <= 2000):
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/imp/imp8/mag/mag_320msec_cdaweb/' + str(imp8Year) + '/i8_320msec_mag_' + str(imp8Date) + '*.cdf'
            else:
                fullPath = ""
                print(imp8Set + '-' + imp8Year + ' IMP8 magnetic data is not available. Only 320msec and 15sec (1973-2000).')
                return fullPath
        elif imp8Dat == 'plasma':
            if int(imp8Year) >= 1973 and int(imp8Year) <= 2006:
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/imp/imp8/plasma_mit/mitplasma_h0/' + str(imp8Year) + '/i8_h0_mitplasma_' + str(imp8Date) + '*.cdf'
            else:
                fullPath = ""
                print(imp8Set + '-' + imp8Year + ' IMP8 magnetic data is not available. Only 4min (1973-2006).')
                return fullPath
        if sorted(glob.glob(dataLoc + fullPath)) == []:
            print("IMP8 data could not be found locally, try to download IMP8 data from destination ....")
            if checkConnection():
                webPath = "ftp://" + fullPath[:len(fullPath)]
                print("Connect to " + webPath)
                envComm = "wget -P " + dataLoc + " -r -l1 -A.cdf " + webPath
                print(envComm)
                os.system(envComm)
            else:
                print("No internet connection found ... Check you network settings.")

        localPath = dataLoc + fullPath
        fnames.extend(glob.glob(localPath))

    return sorted(fnames)

def getIMP8params(imp8Dates,dataLoc,imp8Set='15sec',imp8Dat='magnetic'):
    imp8Set = imp8Set.lower()
    imp8Dat = imp8Dat.lower()

    if imp8Dat == 'magnetic':
        IMP8fnames  = getIMP8files(imp8Dates,dataLoc,imp8Set,imp8Dat='magnetic')
    elif imp8Dat == 'plasma':
        IMP8fnames  = getIMP8files(imp8Dates,dataLoc,imp8Set,imp8Dat='plasma')
    if len(IMP8fnames) > 0:
        cdfKeys = pycdf.CDF(IMP8fnames[0])
    else:
        cdfKeys = "--empty--"
    return cdfKeys


def getIMP8data(imp8Dates, dataLoc, dataList, imp8Set='15sec', dataStat='raw'):
    imp8Set = imp8Set.lower()
    dataStat = dataStat.lower()

    IMP8Bfnames = getIMP8files(imp8Dates,dataLoc,imp8Set,imp8Dat='magnetic')
    IMP8Bparams = getIMP8params(imp8Dates,dataLoc,imp8Set,imp8Dat='magnetic')
    IMP8Pfnames = getIMP8files(imp8Dates,dataLoc,imp8Set,imp8Dat='plasma')
    IMP8Pparams = getIMP8params(imp8Dates,dataLoc,imp8Set,imp8Dat='plasma')

    dataStorage = []
    for j in range(len(dataList)+2):
        dataStorage.append([])

    bFlage = True
    pFlage = True
    for j in range(len(dataList)):
    	if dataList[j] in IMP8Bparams:
    		for i in range(len(IMP8Bfnames)):
    			cdf = pycdf.CDF(IMP8Bfnames[i])
    			if bFlage: dataStorage[len(dataList)].extend(list(cdf['Epoch']))
    			fData = list(cdf[dataList[j]])
    			dataStorage[j].extend(fData)
    		if bFlage: bFlage = False
    		cdf.close()
    	elif dataList[j] in IMP8Pparams:
    		for i in range(len(IMP8Pfnames)):
    			cdf = pycdf.CDF(IMP8Pfnames[i])
    			if pFlage: dataStorage[len(dataList)+1].extend(list(cdf['Epoch']))
    			fData = list(cdf[dataList[j]])
    			dataStorage[j].extend(fData)
    		if pFlage: pFlage = False
    		cdf.close()
    	else:
    	    print dataList[j], ' is not available in IMP8 CDF File!'

    sEpochBID  = bisect.bisect_left(dataStorage[len(dataList)], imp8Dates[0])
    eEpochBID  = bisect.bisect_left(dataStorage[len(dataList)], imp8Dates[len(imp8Dates)-1])
    sEpochPID  = bisect.bisect_left(dataStorage[len(dataList)+1], imp8Dates[0])
    eEpochPID  = bisect.bisect_left(dataStorage[len(dataList)+1], imp8Dates[len(imp8Dates)-1])

    swData = {'magneticEpoch': dataStorage[len(dataList)][sEpochBID:eEpochBID],
              'plasmaEpoch':dataStorage[len(dataList)+1][sEpochPID:eEpochPID]}

    for j in range(len(dataList)):
        if dataList[j] in ['B_Vector_GSE', 'xyzgse', 'SC_Pos_GSE', 'SC_Pos_GSM']:
            swDataTMP = np.transpose(dataStorage[j])
            if dataList[j] in IMP8Bparams:
                swData[dataList[j]] = [swDataTMP[0][sEpochBID:eEpochBID], swDataTMP[1][sEpochBID:eEpochBID], swDataTMP[2][sEpochBID:eEpochBID]]
            elif dataList[j] in IMP8Pparams:
                swData[dataList[j]] = [swDataTMP[0][sEpochPID:eEpochPID], swDataTMP[1][sEpochPID:eEpochPID], swDataTMP[2][sEpochPID:eEpochPID]]
            if dataList[j] == 'B_Vector_GSE' and dataStat == 'clean':
                swData[dataList[j]][0] = dataClean(swData[dataList[j]][0], [9999, -9999], ['>=', '<='])
                swData[dataList[j]][1] = dataClean(swData[dataList[j]][1], [9999, -9999], ['>=', '<='])
                swData[dataList[j]][2] = dataClean(swData[dataList[j]][2], [9999, -9999], ['>=', '<='])
        else:
            if dataList[j] in IMP8Bparams:
                swData[dataList[j]] = dataStorage[j][sEpochBID:eEpochBID]
            elif dataList[j] in IMP8Pparams:
                swData[dataList[j]] = dataStorage[j][sEpochPID:eEpochPID]
            if dataStat == 'clean':
                if dataList[j] == 'proton_density_fit':
                    swData[dataList[j]] = dataClean(swData[dataList[j]], [99, 0], ['>=', '<='])
                elif dataList[j] == 'protonV_thermal_fit':
                    swData[dataList[j]] = dataClean(swData[dataList[j]], [79, 0], ['>=', '<='])
                elif dataList[j] == 'protonV_thermal_mom':
                    swData[dataList[j]] = dataClean(swData[dataList[j]], [79, 0], ['>=', '<='])
                elif dataList[j] == 'V_fit':
                    swData[dataList[j]] = dataClean(swData[dataList[j]], [2500, 100], ['>=', '<='])

    return swData


def imp8DataAdjust(Data,epoch=[],epochLen='not fixed'):
    from scipy.constants import proton_mass as mi
    from scipy.constants import k as KB
    BB = []
    RE = 6371.0
    if epochLen == 'fixed' and epoch == Data['plasmaEpoch']:
        DataAdj={'epoch':epoch}
        DataAdj['Bx']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['B_Vector_GSE'][0],interpKind='cubic'))
        DataAdj['By']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['B_Vector_GSE'][1],interpKind='cubic'))
        DataAdj['Bz']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['B_Vector_GSE'][2],interpKind='cubic'))
        for i in range(len(DataAdj['Bx'])):
            BB.extend([np.sqrt(DataAdj['Bx'][i]**2+DataAdj['By'][i]**2+DataAdj['Bz'][i]**2)])
        DataAdj['B']      = np.array(BB)
        DataAdj['SCxGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_Pos_GSE'][0],interpKind='cubic'))*RE
        DataAdj['SCyGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_Pos_GSE'][1],interpKind='cubic'))*RE
        DataAdj['SCzGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_Pos_GSE'][2],interpKind='cubic'))*RE
        DataAdj['N']      = np.array(Data['proton_density_fit'])
        DataAdj['V']      = np.array(Data['V_fit'])
        DataAdj['T']      = (mi/KB)*(array(Data['protonV_thermal_fit'])*1000.0)**2
    elif epochLen == 'fixed' and epoch == Data['magneticEpoch']:
        DataAdj={'epoch':epoch}
        DataAdj['Bx']     = np.array(Data['B_Vector_GSE'][0])
        DataAdj['By']     = np.array(Data['B_Vector_GSE'][1])
        DataAdj['Bz']     = np.array(Data['B_Vector_GSE'][2])
        for i in range(len(DataAdj['Bx'])):
            BB.extend([np.sqrt(DataAdj['Bx'][i]**2+DataAdj['By'][i]**2+DataAdj['Bz'][i]**2)])
        DataAdj['B']      = np.array(BB)
        DataAdj['SCxGSE'] = np.array(Data['SC_Pos_GSE'][0])*RE
        DataAdj['SCyGSE'] = np.array(Data['SC_Pos_GSE'][1])*RE
        DataAdj['SCzGSE'] = np.array(Data['SC_Pos_GSE'][2])*RE
        DataAdj['N']      = np.array(epochMatch(aceDataAdj['epoch'],Data['plasmaEpoch'],Data['proton_density_fit'],interpKind='linear'))
        DataAdj['V']      = np.array(epochMatch(aceDataAdj['epoch'],Data['plasmaEpoch'],Data['V_fit'],interpKind='linear'))
        DataAdj['T']      = np.array(epochMatch(aceDataAdj['epoch'],Data['plasmaEpoch'],(mi/KB)*(Data['protonV_thermal_fit']*1000.0)**2,interpKind='linear'))
    elif epochLen == 'fixed' and epoch != []:
        DataAdj={'epoch':epoch}
        DataAdj['Bx']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['B_Vector_GSE'][0],interpKind='linear'))
        DataAdj['By']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['B_Vector_GSE'][1],interpKind='linear'))
        DataAdj['Bz']     = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['B_Vector_GSE'][2],interpKind='linear'))
        for i in range(len(DataAdj['Bx'])):
            BB.extend([np.sqrt(DataAdj['Bx'][i]**2+DataAdj['By'][i]**2+DataAdj['Bz'][i]**2)])
        DataAdj['B']      = np.array(BB)
        DataAdj['SCxGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_Pos_GSE'][0],interpKind='linear'))*RE
        DataAdj['SCyGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_Pos_GSE'][1],interpKind='linear'))*RE
        DataAdj['SCzGSE'] = np.array(epochMatch(DataAdj['epoch'],Data['magneticEpoch'],Data['SC_Pos_GSE'][2],interpKind='linear'))*RE
        DataAdj['N']      = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['proton_density_fit'],interpKind='linear'))
        DataAdj['V']      = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],Data['V_fit'],interpKind='linear'))
        DataAdj['T']      = np.array(epochMatch(DataAdj['epoch'],Data['plasmaEpoch'],(mi/KB)*(Data['protonV_thermal_fit']*1000.0)**2,interpKind='linear'))
    else:
        DataAdj={'plasmaEpoch':Data['plasmaEpoch'],'magneticEpoch':Data['magneticEpoch']}
        DataAdj['Bx']     = np.array(Data['B_Vector_GSE'][0])
        DataAdj['By']     = np.array(Data['B_Vector_GSE'][1])
        DataAdj['Bz']     = np.array(Data['B_Vector_GSE'][2])
        for i in range(len(DataAdj['Bx'])):
            BB.extend([np.sqrt(DataAdj['Bx'][i]**2+DataAdj['By'][i]**2+DataAdj['Bz'][i]**2)])
        DataAdj['B']      = np.array(BB)
        DataAdj['SCxGSE'] = np.array(Data['SC_Pos_GSE'][0])*RE
        DataAdj['SCyGSE'] = np.array(Data['SC_Pos_GSE'][1])*RE
        DataAdj['SCzGSE'] = np.array(Data['SC_Pos_GSE'][2])*RE
        DataAdj['N']      = np.array(Data['proton_density_fit'])
        DataAdj['V']      = np.array(Data['V_fit'])
        DataAdj['T']      = (mi/KB)*(array(Data['protonV_thermal_fit'])*1000.0)**2

    return DataAdj
    

def getWINDfiles(windDates,dataLoc,windSet='1min',windDat='merged'):
    windSet = windSet.lower()
    windDat = windDat.lower()

    if dataLoc[len(dataLoc)-1] == "/":
        dataLoc = dataLoc[:len(dataLoc)]
    else:
        dataLoc = dataLoc[:len(dataLoc)] + "/"

    fnames = []

    for i in range(len(windDates)):
        if i > 0:
            if windDates[i] == windDates[i-1]:
                continue
            elif windDates[i].year == windDates[i-1].year:
                if windDates[i].month == windDates[i-1].month:
                    if windDates[i].day == windDates[i-1].day:
                        continue

        windYear = str(windDates[i].year)
        windDate = str(windDates[i].year) + "%02d" % (windDates[i].month) + "%02d" % (windDates[i].day)

        if windDat == 'merged' and int(windYear) >= 1994:
            fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/wind/swe/swe_h1/' + str(windYear) + '/wi_h1_swe_' + str(windDate) + '*.cdf'

        elif windDat == 'magnetic' and int(windYear) >= 1994:
            if windSet == '1min':
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/wind/mfi/mfi_h0/' + str(windYear) + '/wi_h0_mfi_' + str(windDate) + '*.cdf'
            elif windSet == '1sec':
                fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/wind/mfi/mfi_h2/' + str(windYear) + '/wi_h2_mfi_' + str(windDate) + '*.cdf'
            else:
                fullPath = ""
                print(windSet + '-' + windYear + ' WIND magnetic data is not available. Only 1min (1994:).')
                return fullPath

        elif windDat == 'plasma' and (int(windYear) >= 1994 and int(windYear) <= 2001):
            fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/wind/swe/swe_h0/' + str(windYear) + '/wi_h0_swe_' + str(windDate) + '*.cdf'
        elif windDat == 'plasma' and int(windYear) >= 2002:
            fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/wind/swe/swe_h5/' + str(windYear) + '/wi_h5_swe_' + str(windDate) + '*.cdf'
        else:
            fullPath = ""
            print(windSet + '-' + windYear + ' WIND Electron PLasma data is not available. Only 1 min (2002:).')
            return fullPath

        if sorted(glob.glob(dataLoc + fullPath)) == []:
            print("WIND data could not be found locally, try to download WIND data from destination ....")
            if checkConnection():
                webPath = "ftp://" + fullPath[:len(fullPath)]
                print("Connect to " + webPath)
                envComm = "wget -P " + dataLoc + " -r -l1 -A.cdf " + webPath
                print(envComm)
                os.system(envComm)
            else:
                print("No internet connection found ... Check you network settings.")

        localPath = dataLoc + fullPath
        fnames.extend(glob.glob(localPath))

    return sorted(fnames)


def getWINDparams(windDates,dataLoc,windSet='1min',windDat='merged'):
    windSet = windSet.lower()
    windDat = windDat.lower()

    if windDat == 'merged':
        WINDfnames  = getWINDfiles(windDates,dataLoc,windSet,windDat='merged')
    elif windDat == 'magnetic':
        WINDfnames  = getWINDfiles(windDates,dataLoc,windSet,windDat='magnetic')
    elif windDat == 'plasma':
        WINDfnames  = getWINDfiles(windDates,dataLoc,windSet,windDat='plasma')

    if len(WINDfnames) > 0:
        cdfKeys = pycdf.CDF(WINDfnames[0])
    else:
        cdfKeys = "--empty--"
    return cdfKeys


def getWINDdata(windDates,dataLoc,dataList,windSet='1min',windDat='merged',dataStat='raw'):
    windSet = windSet.lower()
    windDat = windDat.lower()

    WINDparams = []; WINDBparams = []; WINDPparams = []
    if windDat == 'merged':
    	WINDfnames = getWINDfiles(windDates,dataLoc,windSet)
    	WINDparams = getWINDparams(windDates,dataLoc,windSet)
    else:
    	WINDBfnames = getWINDfiles(windDates,dataLoc,windSet,windDat='magnetic')
    	WINDBparams = getWINDparams(windDates,dataLoc,windSet,windDat='magnetic')
    	WINDPfnames = getWINDfiles(windDates,dataLoc,windSet,windDat='Plasma')
    	WINDPparams = getWINDparams(windDates,dataLoc,windSet,windDat='Plasma')

    dataStorage = []
    for j in range(len(dataList)+2):
        dataStorage.append([])

    aFlage = True
    bFlage = True
    pFlage = True
    for j in range(len(dataList)):
    	if dataList[j] in WINDparams:
    		for i in range(len(WINDfnames)):
    			cdf = pycdf.CDF(WINDfnames[i])
    			dataStorage[j].extend(list(cdf[dataList[j]]))
    			if aFlage: dataStorage[len(dataList)].extend(list(cdf['Epoch']))
    			if aFlage: dataStorage[len(dataList)+1].extend(list(cdf['Epoch']))
    		if aFlage: aFlage = False
    		cdf.close()
    	elif dataList[j] in WINDBparams:
    		for i in range(len(WINDBfnames)):
    			cdf = pycdf.CDF(WINDBfnames[i])
    			if bFlage: dataStorage[len(dataList)].extend(list(cdf['Epoch']))
    			fData = list(cdf[dataList[j]])
    			dataStorage[j].extend(fData)
    		if bFlage: bFlage = False
    		cdf.close()
    	elif dataList[j] in WINDPparams:
    		for i in range(len(WINDPfnames)):
    			cdf = pycdf.CDF(WINDPfnames[i])
    			if pFlage: dataStorage[len(dataList)+1].extend(list(cdf['Epoch']))
    			fData = list(cdf[dataList[j]])
    			dataStorage[j].extend(fData)
    		if pFlage: pFlage = False
    		cdf.close()
    	else:
    	    print dataList[j], ' is not available in OMNI CDF File!'

    swData = {'magneticEpoch':dataStorage[len(dataList)], 'plasmaEpoch':dataStorage[len(dataList)+1]}
    for j in range(len(dataList)):
    	if windDat == 'merged' and dataStat == 'clean':
    		if dataList[j] == 'Proton_V_nonlin': dataStorage[j] = dataClean(dataStorage[j], [2500], ['>='])
    		if dataList[j] == 'Proton_Np_nonlin': dataStorage[j] = dataClean(dataStorage[j], [999.0, 0.0], ['>=', '<='])
    		if dataList[j] == 'BZ': dataStorage[j] = dataClean(dataStorage[j], [999.0], ['>='])
    	swData[dataList[j]] = dataStorage[j]

    return swData


def windDataAdjust(Data,epoch=[],epochLen='not fixed'):
    from scipy.constants import proton_mass as mi
    from scipy.constants import k as KB
    BB = []
    RE = 6371.0
    if epochLen == 'fixed' and epoch != []:
        DataAdj={'epoch':epoch}
        DataAdj['Bx']     = np.array(epochMatch(DataAdj['epoch'],Data['Epoch'],Data['BX'],interpKind='linear'))
        DataAdj['By']     = np.array(epochMatch(DataAdj['epoch'],Data['Epoch'],Data['BY'],interpKind='linear'))
        DataAdj['Bz']     = np.array(epochMatch(DataAdj['epoch'],Data['Epoch'],Data['BZ'],interpKind='linear'))
        for i in range(len(DataAdj['Bx'])):
            BB.extend([np.sqrt(DataAdj['Bx'][i]**2+DataAdj['By'][i]**2+DataAdj['Bz'][i]**2)])
        DataAdj['B']      = np.array(BB)
        DataAdj['N']      = np.array(epochMatch(DataAdj['epoch'],Data['Epoch'],Data['Proton_Np_nonlin'],interpKind='linear'))
        DataAdj['V']      = np.array(epochMatch(DataAdj['epoch'],Data['Epoch'],Data['Proton_V_nonlin'],interpKind='linear'))
        DataAdj['T']      = np.array(epochMatch(DataAdj['epoch'],Data['Epoch'],(mi/KB)*(Data['Proton_W_nonlin']/1000.0)**2,interpKind='linear'))
        DataAdj['SCxGSE'] = np.array(Data['xgse'])*RE
        DataAdj['SCyGSE'] = np.array(Data['ygse'])*RE
        DataAdj['SCzGSE'] = np.array(Data['zgse'])*RE
        DataAdj['Vx']     = np.array(epochMatch(DataAdj['epoch'],Data['Epoch'],Data['Proton_VX_nonlin'],interpKind='linear'))
        DataAdj['Vy']     = np.array(epochMatch(DataAdj['epoch'],Data['Epoch'],Data['Proton_VY_nonlin'],interpKind='linear'))
        DataAdj['Vz']     = np.array(epochMatch(DataAdj['epoch'],Data['Epoch'],Data['Proton_VZ_nonlin'],interpKind='linear'))
    else:
        DataAdj={'epoch':Data['Epoch']}
        DataAdj['Bx']     = np.array(Data['BX'])
        DataAdj['By']     = np.array(Data['BY'])
        DataAdj['Bz']     = np.array(Data['BZ'])
        for i in range(len(DataAdj['Bx'])):
            BB.extend([np.sqrt(DataAdj['Bx'][i]**2+DataAdj['By'][i]**2+DataAdj['Bz'][i]**2)])
        DataAdj['B']      = np.array(BB)
        DataAdj['N']      = np.array(Data['Proton_Np_nonlin'])
        DataAdj['V']      = np.array(Data['Proton_V_nonlin'])
        DataAdj['T']      = (mi/KB)*(array(Data['Proton_W_nonlin'])/1000.0)**2
        DataAdj['SCxGSE'] = np.array(Data['xgse'])*RE
        DataAdj['SCyGSE'] = np.array(Data['ygse'])*RE
        DataAdj['SCzGSE'] = np.array(Data['zgse'])*RE
        DataAdj['Vx']     = np.array(Data['Proton_VX_nonlin'])
        DataAdj['Vy']     = np.array(Data['Proton_VY_nonlin'])
        DataAdj['Vz']     = np.array(Data['Proton_VZ_nonlin'])

    return DataAdj


def getGeotailfiles(geotailDates,dataLoc,geotailSet='1min'):
    geotailSet = geotailSet.lower()

    if dataLoc[len(dataLoc)-1] == "/":
        dataLoc = dataLoc[:len(dataLoc)]
    else:
        dataLoc = dataLoc[:len(dataLoc)] + "/"

    fnames = []

    for i in range(len(geotailDates)):
        if i > 0:
            if geotailDates[i] == geotailDates[i-1]:
             continue
            elif geotailDates[i].year == geotailDates[i-1].year:
             if geotailDates[i].month == geotailDates[i-1].month:
              if geotailDates[i].day == geotailDates[i-1].day:
               continue

        geotailYear = str(geotailDates[i].year)
        geotailDate = str(geotailDates[i].year) + "%02d" % (geotailDates[i].month) + "%02d" % (geotailDates[i].day)

        if geotailSet=='1min' and (int(geotailYear) >= 1995 and int(geotailYear) <= 2005):
            fullPath = 'cdaweb.gsfc.nasa.gov/pub/data/geotail/merged/magplasma_sw_1min_cdaweb/' + str(geotailYear) + \
                          '/ge_1min_mag_plasma_sw_only_' + str(geotailDate) + '*.cdf'
        else:
            fullPath = ""
            print(geotailSet + '-' + geotailYear + ' Geotail data is not available. Only 1min (1995:2005).')
            return fullPath

        if sorted(glob.glob(dataLoc + fullPath)) == []:
            print("Geotail data could not be found locally, try to download OMNI data from destination ....")
            if checkConnection():
                webPath = "ftp://" + fullPath[:len(fullPath)]
                print("Connect to " + webPath)
                envComm = "wget -P /home/ehab/SWData" + " -r -l1 -A.cdf " + webPath
                print(envComm)
                os.system(envComm)
            else:
                print("No internet connection found ... Check you network settings.")

        localPath = dataLoc + fullPath
        fnames.extend(glob.glob(localPath))

    return sorted(fnames)


def getGeotailparams(geotailDates,dataLoc,geotailSet='1min'):
    Geotailfnames  = getGeotailfiles(geotailDates,dataLoc,geotailSet)
    if len(Geotailfnames) > 0:
     cdfKeys = pycdf.CDF(Geotailfnames[0])
    else:
     cdfKeys = "--empty--"
    return cdfKeys


def getGeotaildata(geotailDates,dataLoc,dataList,geotailSet='1min',dataStat='raw'):
    Geotailfnames = getGeotailfiles(geotailDates,dataLoc,geotailSet)
    Geotailparams = getGeotailparams(geotailDates,dataLoc,geotailSet)

    dataStorage = []
    for j in range(len(dataList)+1):
        dataStorage.append([])

    aFlage = True
    for j in range(len(dataList)):
    	if dataList[j] in Geotailparams:
    		for i in range(len(Geotailfnames)):
    			cdf = pycdf.CDF(Geotailfnames[i])
    			dataStorage[j].extend(list(cdf[dataList[j]]))
    			if aFlage: dataStorage[len(dataList)].extend(list(cdf['Epoch']))
    		if aFlage: aFlage = False
    		cdf.close()
    	else:
    	    print dataList[j], ' is not available in Geotail CDF File!'

    swData = {'Epoch': dataStorage[len(dataList)]}
    for j in range(len(dataList)):
    	if dataStat == 'clean':
    		if dataList[j] == 'V': dataStorage[j] = dataClean(dataStorage[j], [2500,-2500], ['>=','<='])
    		if dataList[j] == 'N': dataStorage[j] = dataClean(dataStorage[j], [999.0,0.0], ['>=','<='])
    		if dataList[j] == 'T': dataStorage[j] = dataClean(dataStorage[j], [1.0e7], ['>='])
    		if dataList[j] == 'ABS_B': dataStorage[j] = dataClean(dataStorage[j], [999.0,-999.0], ['>=','<='])
    		if dataList[j] == 'BX_GSE': dataStorage[j] = dataClean(dataStorage[j], [999.0,-999.0], ['>=','<='])
    		if dataList[j] == 'BY_GSE': dataStorage[j] = dataClean(dataStorage[j], [999.0,-999.0], ['>=','<='])
    		if dataList[j] == 'BZ_GSE': dataStorage[j] = dataClean(dataStorage[j], [999.0,-999.0], ['>=','<='])
    	swData[dataList[j]] = dataStorage[j]

    return swData


def geotailDataAdjust(Data,epoch=[],epochLen='not fixed'):
    BB = []
    RE = 6371.0
    if epochLen == 'fixed' and epoch != []:
     DataAdj={'epoch':epoch}
     DataAdj['B']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['ABS_B'], interpKind='linear'))
     DataAdj['Bx']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['BX_GSE'], interpKind='linear'))
     DataAdj['By']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['BY_GSE'], interpKind='linear'))
     DataAdj['Bz']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['BZ_GSE'], interpKind='linear'))
     DataAdj['N']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['N'], interpKind='linear'))
     DataAdj['V']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['V'], interpKind='linear'))
     DataAdj['T']      = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['T'], interpKind='linear'))
     DataAdj['SCxGSE'] = np.array(Data['X'])*RE
     DataAdj['SCyGSE'] = np.array(Data['Y'])*RE
     DataAdj['SCzGSE'] = np.array(Data['Z'])*RE
     DataAdj['Vx']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['VX_GSE'], interpKind='linear'))
     DataAdj['Vy']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['VY_GSE'], interpKind='linear'))
     DataAdj['Vz']     = np.array(epochMatch(DataAdj['epoch'], Data['Epoch'], Data['VZ_GSE'], interpKind='linear'))
    else:
     DataAdj={'epoch':Data['Epoch']}
     DataAdj['B']      = np.array(Data['ABS_B'])
     DataAdj['Bx']     = np.array(Data['BX_GSE'])
     DataAdj['By']     = np.array(Data['BY_GSE'])
     DataAdj['Bz']     = np.array(Data['BZ_GSE'])
     DataAdj['N']      = np.array(Data['N'])
     DataAdj['V']      = np.array(Data['V'])
     DataAdj['T']      = np.array(Data['T'])
     DataAdj['SCxGSE'] = np.array(Data['X'])*RE
     DataAdj['SCyGSE'] = np.array(Data['Y'])*RE
     DataAdj['SCzGSE'] = np.array(Data['Z'])*RE
     DataAdj['Vx']     = np.array(Data['VX_GSE'])
     DataAdj['Vy']     = np.array(Data['VY_GSE'])
     DataAdj['Vz']     = np.array(Data['VZ_GSE'])

     return DataAdj


def dataClean(myList,threshold,operand):
    for i in range(len(myList)):
    	for j in range(len(threshold)):
    		if operand[j] == '>=':
    			if myList[i] >= threshold[j]:
    			    myList[i] = float('nan')
    		elif operand[j] == '<=':
    			if myList[i] <= threshold[j]:
    			    myList[i] = float('nan')
    		elif operand[j] == '==':
    			if myList[i] == threshold[j]:
    			    myList[i] = float('nan')
    		elif operand[j] == '<':
    			if myList[i] < threshold[j]:
    			    myList[i] = float('nan')
    		elif operand[j] == '>':
    			if myList[i] > threshold[j]:
    			    myList[i] = float('nan')
    return myList


def dateShift(oldDate, years=0, months=0, days=0, hours=0, minutes=0, seconds=0):
    month = oldDate.month - 1 + months
    year = oldDate.year + month / 12 + years
    month = month % 12 + 1
    day = oldDate.day
    hour = oldDate.hour
    minute = oldDate.minute
    second = oldDate.second
    newDate = dt.datetime(year,month,day,hour,minute,second) + dt.timedelta(days=days,hours=hours,minutes=minutes,seconds=seconds)

    return newDate


def dateList(sDate, eDate, shift = 'day'):
    shift = shift.lower()
    dList = []
    if shift == 'day' or shift == 'days':
    	diff = eDate - sDate
    	for i in range(diff.days + 1):
    	    dList = dList + [dateShift(sDate, days = i)]
    elif shift == 'month' or shift == 'months':
    	if eDate.year - sDate.year == 0:
    	    diff = eDate.month - sDate.month
    	elif eDate.year - sDate.year != 0:
    	    diff = (eDate.year - sDate.year) * 12 + abs(eDate.month - sDate.month)
    	for i in range(diff + 1):
    	    dList = dList + [dateShift(sDate, months = i)]
    elif shift == 'year' or shift == 'years':
    	diff = eDate.year - sDate.year
    	for i in range(diff + 1):
    	    dList = dList + [dateShift(sDate, years = i)]
    elif shift == 'hour' or shift == 'hours':
    	diff = eDate - sDate
    	for i in range(diff.days*24 + diff.seconds / 3600 + 1):
    	    dList = dList + [dateShift(sDate, hours = i)]
    elif shift == 'minute' or shift == 'minutes':
    	diff = eDate - sDate
    	for i in range(diff.days*(24*60) + diff.seconds / 60 + 1):
    	    dList = dList + [dateShift(sDate, minutes = i)]
    elif shift == 'second' or shift == 'seconds':
    	diff = eDate - sDate
    	for i in range(diff.days*(24*60*60) + diff.seconds + 1):
    	    dList = dList + [dateShift(sDate, seconds = i)]

    return dList   


def removeNaN(data,epoch=[]):
    if epoch == [] and data != []:
    	newData=[]
    	for i in range(len(data)):
    		if str(data[i]) != 'nan':
    		    newData.extend([data[i]])
    	return newData, []
    elif epoch != [] and len(data) == len(epoch):
    	newEpoch=[]; newData=[]
    	for i in range(len(epoch)):
    		if str(data[i]) != 'nan':
    		    newEpoch.extend([epoch[i]])
    		    newData.extend([data[i]])
    	return newEpoch,newData
    else:
        return [], []

'''
def removeNaN(data,epoch=[]):
    newEpoch=[]; newData=[]
    for i in range(len(epoch)):
     if str(data[i]) != 'nan':
      newEpoch.extend([epoch[i]])
      newData.extend([data[i]])
    return newEpoch,newData
'''

def mapDataToEpoch(refEpoch, epoch, data, interpKind='linear'):
    from time import mktime
    from scipy.interpolate import interp1d

    interpKind = interpKind.lower()
    refEpochStamp = []
    for i in range(len(refEpoch)):
        refEpochStamp.extend([mktime(refEpoch[i].timetuple())])
    epochStamp = []
    for i in range(len(epoch)):
        epochStamp.extend([mktime(epoch[i].timetuple())])

    f = interp1d(epochStamp, data, bounds_error=False, fill_value=0.0, kind = interpKind)

    return f(refEpochStamp)

#def commonEpoch(epoch1,epoch2):
#   from bisect import bisect_left
#   if epoch1[0] <= epoch2[0]:
#    sEpoch = bisect_left(epoch1, epoch2[0])
#   elif epoch1[0] > epoch2[0]:
#    sEpoch = bisect_left(epoch1, epoch1[0])

#   if epoch1[-1] <= epoch2[-1]:
#    eEpoch = bisect_left(epoch1, epoch1[-1])
#   elif epoch1[-1] > epoch2[-1]:
#    eEpoch = bisect_left(epoch1, epoch2[-1])

#   return epoch1[sEpoch:eEpoch+1], range(sEpoch,eEpoch+1,1)

def commonEpoch(epoch1,epoch2):
    from swdatanal import search
    cmnEpoch = []
    indEpoch = []
    for epoch in epoch2:
    	if epoch <= epoch1[-1] and epoch >= epoch1[0]:
    		ind = bisect.bisect_left(epoch1, epoch)
    		target = epoch1[ind]
    		if search(cmnEpoch,epoch1[ind]) == []:
    			cmnEpoch.extend([epoch1[ind]])
    			indEpoch.extend([ind])
    return cmnEpoch, indEpoch


def epochShift(usEpoch,sLag):
    sEpoch = []
    for i in range(len(usEpoch)):
        sEpoch.extend([usEpoch[i] + dt.timedelta(0, sLag[i])])
    return sEpoch

#def epochShift(usEpoch,usParam,sLag):
#   from bisect import bisect_left
#   from numpy import zeros
#   from datetime import timedelta
#   cEpoch = []
#   for i in range(len(usEpoch)):
#    cEpoch.extend([usEpoch[i] + timedelta(0,sLag[i])])
#   if cEpoch[0] >= usEpoch[0]:
#    sEpochID  = bisect_left(usEpoch, cEpoch[0])
#    eEpochID  = bisect_left(usEpoch,usEpoch[-1])
#   elif cEpoch[0] < usEpoch[0]:
#    sEpochID  = bisect_left(usEpoch,usEpoch[0])
#    eEpochID  = bisect_left(usEpoch, cEpoch[-1])

#   sEpoch = usEpoch[sEpochID:eEpochID]
#   sParam = usParam[sEpochID:eEpochID]

#   return sEpoch,sParam


 
