
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
    from datetime import datetime, timedelta

    month = oldDate.month - 1 + months
    year = oldDate.year + month / 12 + years
    month = month % 12 + 1
    day = oldDate.day
    hour = oldDate.hour
    minute = oldDate.minute
    second = oldDate.second
    newDate = datetime(year,month,day,hour,minute,second) + timedelta(days=days,hours=hours,minutes=minutes,seconds=seconds)

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
     for i in range(diff.seconds / 3600 + 1):
      dList = dList + [dateShift(sDate, hours = i)]
    elif shift == 'minute' or shift == 'minutes':
     diff = eDate - sDate
     for i in range(diff.seconds / 60 + 1):
      dList = dList + [dateShift(sDate, minutes = i)]
    elif shift == 'second' or shift == 'seconds':
     diff = eDate - sDate
     for i in range(diff.seconds + 1):
      dList = dList + [dateShift(sDate, seconds = i)]

    return dList 
     
