
def checkConnection():
    import httplib
    conn = httplib.HTTPConnection("www.google.com")
    try:
        conn.request("HEAD", "/")
        return True
    except:
        return False

def getOMNIfiles(omniDates,dataLoc,omniSet='hourly'):
    import os
    import glob

    if dataLoc[len(dataLoc)-1] == "/":
     dataLoc = dataLoc[:len(dataLoc)]
    else:
     dataLoc = dataLoc[:len(dataLoc)] + "/"

    fnames = []

    for i in range(len(omniDates)):
     if len(omniDates[i]) <= 4:
      omniYear = omniDates[i][:4]
      omniDate = omniYear
     elif len(omniDates[i]) > 4: 
      omniYear  = omniDates[i][:4]
      if omniSet == 'hourly' and int(omniDates[i][5:7]) <= 6:
       omniDate = omniYear + "0101"
      elif omniSet == 'hourly' and int(omniDates[i][5:7]) >= 7:
       omniDate = omniYear + "0701"
      else:
       omniDate = omniYear + omniDates[i][5:7] + "01"

     if omniYear == omniDate:
      if omniSet=='hourly':
       fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hourly/' + str(omniYear) + '/*.cdf'
      elif omniSet=='5min':
       fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_5min/' + str(omniYear) + '/*.cdf'
      elif omniSet=='1min':
       fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/' + str(omniYear) + '/*.cdf'
      else:
       fullPath = ""
       print(omniSet + ' OMNI data is not available.')
     else:
      if omniSet=='hourly':
       fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hourly/' + str(omniYear) + '/omni2_h0_mrg1hr_' + str(omniDate) + '_v01.cdf'
      elif omniSet=='5min':
       fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_5min/' + str(omniYear) + '/omni_hro_5min_' + str(omniDate) + '_v01.cdf'
      elif omniSet=='1min':
       fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/' + str(omniYear) + '/omni_hro_1min_' + str(omniDate) + '_v01.cdf'
      else:
       fullPath = ""
       print(omniSet + ' OMNI data is not available.')

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


def dataClean(oldList,threshold):
    newList = []
    for i in oldList:
     if i <= threshold:
      newList.extend([i]) 
     else:
      newList.extend([float('nan')]) 
    return newList

