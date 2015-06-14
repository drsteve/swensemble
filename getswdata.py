
def checkConnection():
    import httplib
    conn = httplib.HTTPConnection("www.google.com")
    try:
        conn.request("HEAD", "/")
        return True
    except:
        return False

def getOMNIfiles(omniYear,dataLoc,omniSet='hourly'):
    import os
    import glob
    if dataLoc[len(dataLoc)-1] == "/":
     dataLoc = dataLoc[:len(dataLoc)]
    else:
     dataLoc = dataLoc[:len(dataLoc)] + "/"

    if omniSet=='hourly':
     fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hourly/' + str(omniYear) + '/*.cdf'
    elif omniSet=='5min':
     fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_5min/' + str(omniYear) + '/*.cdf'
    elif omniSet=='1min':
     fullPath = 'spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hro_1min/' + str(omniYear) + '/*.cdf'
    else:
     fullPath = ""
     print(omniSet + ' OMNI data is not available.')

    localPath = dataLoc + fullPath
    fnames = sorted(glob.glob(localPath))

    if fnames == []:
     print "OMNI data could not be found locally, try to download OMNI data from destination ...."
     if checkConnection():
      webPath = "ftp://" + fullPath[:len(fullPath)-5]
      print "Connect to " + webPath
      envComm = "wget -P /home/ehab/SWData" + " -r -l1 -A.cdf " + webPath
      print envComm
      os.system(envComm)

    localPath = dataLoc + fullPath
    fnames = sorted(glob.glob(localPath))
    return fnames

