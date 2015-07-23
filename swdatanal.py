def findContiguousData(x, delta, minLength=None):
    import numpy as np
    import datetime as dt
    """Find all intervals of contiguous data in x exceeding min_length

    Contiguous data are defined as neighbouring points separated by delta
    i.e., x[i+1] - x[i] = delta

    If min_length is undefined then min_length = Delta
    
    Inputs:
    =======
    x, input series of times. Can be datetime objects or serial time
    delta, expected resolution of x.
    minLength [defaults to delta], minimum length of contiguous interval required.
    
    Returns:
    ========
    istart, iend - indices of starts and ends of contiguous blocks
    
    Author:
    =======
    Original -- Mervyn Freeman, British Antarctic Survey
    Port to Python by Steve Morley, Los Alamos National Lab.
    smorley@lanl.gov/morley_steve@hotmail.com
    """

    if not minLength:
        minLength = delta

    if type(x)==list:
        x = np.array(x)
    #now ensure type consistency of array contents for datetime input
    if isinstance(x[0], dt.datetime):
        try:
            assert type(delta) == type(minLength)
            assert isinstance(delta, dt.timedelta)
        except:
            return 'findContiguousData: inconsistent data types for time objects'

    #Calculate distance between neighbouring data points in array x
    dx = x[1:]-x[0:-1]

    #Find positions i where X neighbours are non-contiguous
    #i.e., X(i+1) - X(i) > Delta
    #Store in array igaps
    igaps, = (dx > delta).nonzero()
    igaps1, nigaps = igaps+1, len(igaps)

    #Now find intervals of contiguous data exceeding min_length
    #Contiguous data interval starts at the end of a non-contiguous interval
    #and ends at the start of the next non-contiguous interval.

    #Start of series X is start of first potentially contiguous interval
    istart_c = [0]
    #Find starts of other potentially contiguous intervals
    #from ends of non-contiguous intervals
    end_nc = x[igaps+1].tolist()
    start_c = end_nc
    istart_c = igaps1.tolist()
    start_c.insert(0, x[0])
    istart_c.insert(0, 0)
    #Find ends of potentially contiguous intervals
    #from starts of non-contiguous intervals
    end_c = x[igaps].tolist() #start_nc
    iend_c = igaps.tolist()
    #Add end of series X as end of last potentially contiguous interval
    end_c.append(x[-1])
    iend_c.append(len(x)-1)

    #Find lengths of all potentially contiguous intervals
    length = [cEnd - start_c[i] for i, cEnd in enumerate(end_c)]
    #Find those whose length exceeds min_length
    ilong, = (np.array(length) > minLength).nonzero()

    #Return start and end indices of these intervals
    istart = [istart_c[ind] for ind in ilong]
    iend = [iend_c[ind] for ind in ilong]

    return istart, iend

def epochBlock(epoch, data, blockLen, gapLen, fixStep = True):
    from numpy import diff
    from datetime import timedelta
    from bisect import bisect_left
    from getswdata import removeNaN

    if len(epoch) != len(data):
     print '(epoch) length MUST equal (data) length.'
     return ''

    epochLength = epoch[-1] - epoch[0]
    nHours = epochLength.days*24 + epochLength.seconds/3600
    if nHours < blockLen:
     print '(blockLen) is shorter than the available (data) block of time'
     return ''

    blockstart = []
    cEpoch = epoch[0]
    while cEpoch < epoch[-(blockLen-1)]:
     try:
      sEpochID  = bisect_left(epoch, cEpoch + timedelta(0,0))
      eEpochID  = bisect_left(epoch, cEpoch + timedelta(0,blockLen*3600))
     except:
      if fixStep:
       cEpoch = cEpoch + timedelta(0,blockLen*3600)
      else:
       cEpoch = cEpoch + timedelta(0,3600)
      continue
     TT,DD = removeNaN(data[sEpochID:eEpochID],epoch = epoch[sEpochID:eEpochID])
     if diff(TT) != []:
      try:
       dt = max(diff(TT))
       if dt.seconds <= gapLen:
        blockstart.append(cEpoch)
        cEpoch = cEpoch + timedelta(0,blockLen*3600)
       else:
        if fixStep:
         cEpoch = cEpoch + timedelta(0,blockLen*3600)
        else:
         cEpoch = cEpoch + timedelta(0,3600)
      except:
       if fixStep:
        cEpoch = cEpoch + timedelta(0,blockLen*3600)
       else:
        cEpoch = cEpoch + timedelta(0,3600)
     else:
       if fixStep:
        cEpoch = cEpoch + timedelta(0,blockLen*3600)
       else:
        cEpoch = cEpoch + timedelta(0,3600)

    return blockstart

def dataFilter(data,filterVal,filterCondition):
    condStatus = []
    for val in data:
     if filterCondition == '==':
      if val == filterVal: condStatus.extend([val])
     elif filterCondition == '<=':
      if val <= filterVal: condStatus.extend([val])
     elif filterCondition == '>=':
      if val >= filterVal: condStatus.extend([val])
     elif filterCondition == '<':
      if val <  filterVal: condStatus.extend([val])
     elif filterCondition == '>':
      if val >  filterVal: condStatus.extend([val])
    return condStatus


def findCorrEpoch(epoch1,epoch2):
    from datetime import datetime, timedelta
    from numpy import diff,array
    cEpoch = []
    for i in range(len(epoch1)):
     epochDiff = abs(epoch1[i]-array(epoch2))
     minDiff   = min(epochDiff)
     if minDiff.seconds <= 900:
      j = search(epochDiff, minDiff)[0]
      cEpoch.extend([min(epoch1[i],epoch2[j])])
    return cEpoch


def getDistrib(data, nbins=0, stride=0, bins=[], norm = False):
    from scipy.stats import histogram, histogram2
    from numpy import arange

    if nbins>0:
     stride = (max(data)-min(data))/nbins
     bins = arange(min(data)-stride,max(data)+stride,stride)
     dist = histogram2(data,bins)
     if norm:
      dist = map(float, dist)
      dist = [dist[i]/sum(dist) for i in range(len(dist))]
     return dist, bins, stride
    elif stride>0:
     bins = arange(min(data)-stride,max(data)+stride,stride)
     dist = histogram2(data,bins)
     if norm:
      dist = map(float, dist)
      dist = [dist[i]/sum(dist) for i in range(len(dist))]
     return dist, bins
    elif len(bins)>0:
     dist = histogram2(data,bins)
     if norm:
      dist = map(float, dist)
      dist = [dist[i]/sum(dist) for i in range(len(dist))]
     return dist
    else:
     nbins = 10
     stride = (max(data)-min(data))/nbins
     bins = arange(min(data)-stride,max(data)+stride,stride)
     dist = histogram2(data,bins)
     if norm:
      dist = map(float, dist)
      dist = [dist[i]/sum(dist) for i in range(len(dist))]
     return dist, bins


def omniDataCorr(srefDate, erefDate, startDate, endDate, epochs, SWP, binStride, CorrTime = 'Day', CorrType = 'kstest'):
    import numpy
    import bisect
    import datetime
    from scipy.stats import ks_2samp, pearsonr
    from getswdata import getOMNIfiles, dataClean, dateShift, dateList

    CorrTime = CorrTime.lower()
    CorrType = CorrType.lower()

    if endDate < startDate:
     print('(swdatanal.omniDataCorr).Error: Dates are not applicable')
     SWPDatRng=0; cepochs=0; KSVals=0; KSDist=0; aepochs=0
     return SWPDatRng, cepochs, KSVals, KSDist, aepochs

    sEpochID  = bisect.bisect_left(epochs, startDate)
    eEpochID  = bisect.bisect_left(epochs, endDate)
    cepochs   = epochs[sEpochID:eEpochID]
    SWPDatRng = SWP[sEpochID:eEpochID]
    if SWP[sEpochID:eEpochID] == []:
     print('(swdatanal.omniDataCorr).Error: No data avaliable for the designated date(s) and/or time(s).')
     SWPDatRng=0; cepochs=0; KSVals=0; KSDist=0; aepochs=0
     return SWPDatRng, cepochs, KSVals, KSDist, aepochs
    _, bins   = getDistrib(filter(lambda v: v==v, SWPDatRng), stride = binStride, norm = False)

    sEpochID = bisect.bisect_left(epochs, srefDate)
    eEpochID = bisect.bisect_left(epochs, erefDate)
    SWPV01   = SWP[sEpochID:eEpochID]
    SWPD01   = getDistrib(filter(lambda v: v==v, SWPV01), bins=bins, norm=True)

    if CorrTime == 'day':
     aepochs = []; KSVals = []; KSDist = []
     sEpoch = datetime.datetime(startDate.year,startDate.month,startDate.day, 0, 0, 0)
     eEpoch = dateShift(sEpoch, hours = 23, minutes = 59, seconds = 59)
     for i in range((endDate-startDate).days+1):
      aepochs  = aepochs + [dateShift(sEpoch,0,0,i,0,0,0)]
      sEpochID = bisect.bisect_left(epochs, dateShift(sEpoch,0,0,i,0,0,0))
      eEpochID = bisect.bisect_left(epochs, dateShift(eEpoch,0,0,i,0,0,0))

      SWPV02 = SWP[sEpochID:eEpochID]
      SWPD02 = getDistrib(filter(lambda v: v==v, SWPV02), bins=bins, norm=True)

      if CorrType == 'kstest':
       KSVals = KSVals + [ks_2samp(SWPV01, SWPV02)]
       KSDist = KSDist + [ks_2samp(SWPD01, SWPD02)]
      elif CorrType == 'pearson':
       KSVals = KSVals + [pearsonr(SWPV01, SWPV02)]
       KSDist = KSDist + [pearsonr(SWPD01, SWPD02)]

     KSVals = numpy.array(KSVals)
     KSDist = numpy.array(KSDist)

    return SWPDatRng, cepochs, KSVals, KSDist, aepochs
   

def getSolarWindType(swData,nCats = 4, gplot = True):
    from  numpy import array, log10, dot, sign
    import matplotlib.pyplot as plt
    import scipy.constants

    Tp = (scipy.constants.k*array(swData['T']))/scipy.constants.physical_constants['electron volt'][0]
    Sp = Tp/((array(swData['N']))**0.667)
    Tr = ((array(swData['V'])/258.0)**3.113)/Tp
    VA = 21.8*array(swData['B'])/((array(swData['N']))**0.5)

    dx = log10(Sp)
    dy = log10(VA)
    dz = log10(Tr)

    VEJT=[];VCHO=[];VSRR=[];VSBO=[]
    NEJT=[];NCHO=[];NSRR=[];NSBO=[]
    TEJT=[];TCHO=[];TSRR=[];TSBO=[]
    BEJT=[];BCHO=[];BSRR=[];BSBO=[]
    EEJT=[];ECHO=[];ESRR=[];ESBO=[]
    if nCats == 4:
     for i in range(len(dx)):
      if dy[i] > 0.277 * dx[i] + 0.055 * dz[i] + 1.83:           # Ejecta
       VEJT.extend([swData['V'][i]])
       NEJT.extend([swData['N'][i]])
       TEJT.extend([swData['T'][i]])
       BEJT.extend([swData['B'][i]])
       EEJT.extend([swData['epoch'][i]])
      elif dx[i] > -0.525 * dz[i] - 0.676 * dy[i] + 1.74:       # Coronal-Hole_Origin
       VCHO.extend([swData['V'][i]])
       NCHO.extend([swData['N'][i]])
       TCHO.extend([swData['T'][i]])
       BCHO.extend([swData['B'][i]])
       ECHO.extend([swData['epoch'][i]])
      elif dx[i] < -0.658 * dy[i] - 0.125 * dz[i] + 1.04:       # Sector-Reversal-Region
       VSRR.extend([swData['V'][i]])
       NSRR.extend([swData['N'][i]])
       TSRR.extend([swData['T'][i]])
       BSRR.extend([swData['B'][i]])
       ESRR.extend([swData['epoch'][i]])
      else:                                                     # Streamer-Belt-Origin
       VSBO.extend([swData['V'][i]])
       NSBO.extend([swData['N'][i]])
       TSBO.extend([swData['T'][i]])
       BSBO.extend([swData['B'][i]])
       ESBO.extend([swData['epoch'][i]])
    elif nCats == 3:
     for i in range(len(dx)):
      if dy[i] > 0.277 * dx[i] + 0.055 * dz[i] + 1.83:           # Ejecta
       VEJT.extend([swData['V'][i]])
       NEJT.extend([swData['N'][i]])
       TEJT.extend([swData['T'][i]])
       BEJT.extend([swData['B'][i]])
       EEJT.extend([swData['epoch'][i]])
      elif dx[i] > -0.525 * dz[i] - 0.676 * dy[i] + 1.74:       # Coronal-Hole_Origin
       VCHO.extend([swData['V'][i]])
       NCHO.extend([swData['N'][i]])
       TCHO.extend([swData['T'][i]])
       BCHO.extend([swData['B'][i]])
       ECHO.extend([swData['epoch'][i]])
      else:                                                     # Streamer-Belt-Origin
       VSBO.extend([swData['V'][i]])
       NSBO.extend([swData['N'][i]])
       TSBO.extend([swData['T'][i]])
       BSBO.extend([swData['B'][i]])
       ESBO.extend([swData['epoch'][i]])

    swCats = {'Sp':array(Sp), 'Tr':array(Tr), 'VA':array(VA)}

    swCats['VEJT'] = array(VEJT)
    swCats['NEJT'] = array(NEJT)
    swCats['TEJT'] = array(TEJT)
    swCats['BEJT'] = array(BEJT)
    swCats['EEJT'] = array(EEJT)

    swCats['VCHO'] = array(VCHO)
    swCats['NCHO'] = array(NCHO)
    swCats['TCHO'] = array(TCHO)
    swCats['BCHO'] = array(BCHO)
    swCats['ECHO'] = array(ECHO)

    swCats['VSRR'] = array(VSRR)
    swCats['NSRR'] = array(NSRR)
    swCats['TSRR'] = array(TSRR)
    swCats['BSRR'] = array(BSRR)
    swCats['ESRR'] = array(ESRR)

    swCats['VSBO'] = array(VSBO)
    swCats['NSBO'] = array(NSBO)
    swCats['TSBO'] = array(TSBO)
    swCats['BSBO'] = array(BSBO)
    swCats['ESBO'] = array(ESBO)

    if gplot:
     plt.figure(9001)
     plt.subplot(3,1,1)
     plt.plot(swCats['EEJT'],swCats['VEJT'], 'bo', label = 'Ejecta')
     plt.plot(swCats['ECHO'],swCats['VCHO'], 'ro', label = 'Conronal-Hole')
     if nCats == 4: plt.plot(swCats['ESRR'],swCats['VSRR'], 'mo', label = 'Sector-Reversal')
     plt.plot(swCats['ESBO'],swCats['VSBO'], 'go', label = 'Streamer-Belt')
     plt.ylabel('Flow Velocity (km/s)')
     plt.title('Solar Wind Categorization')
     plt.legend(loc='best')
     plt.subplot(3,1,2)
     plt.plot(swCats['EEJT'],swCats['NEJT'], 'bo', label = 'Ejecta')
     plt.plot(swCats['ECHO'],swCats['NCHO'], 'ro', label = 'Conronal-Hole')
     if nCats == 4: plt.plot(swCats['ESRR'],swCats['NSRR'], 'mo', label = 'Sector-Reversal')
     plt.plot(swCats['ESBO'],swCats['NSBO'], 'go', label = 'Streamer-Belt')
     plt.ylabel('Flow Density (N/cc)')
     plt.legend(loc='best')
     plt.subplot(3,1,3)
     plt.plot(swCats['EEJT'],swCats['BEJT'], 'bo', label = 'Ejecta')
     plt.plot(swCats['ECHO'],swCats['BCHO'], 'ro', label = 'Conronal-Hole')
     if nCats == 4: plt.plot(swCats['ESRR'],swCats['BSRR'], 'mo', label = 'Sector-Reversal')
     plt.plot(swCats['ESBO'],swCats['BSBO'], 'go', label = 'Streamer-Belt')
     plt.ylabel('$B_z$ (nT)')
     plt.xlabel('Date')
     plt.legend(loc='best')
     plt.title('Solar Wind Categorization')

   #return list(Sp), list(Tr), list(VA), SWPClass
    return swCats


def getTimeLag(epoch,srcData,destPos,method='flat'):
    from math import atan2, tan, degrees, radians
    from datetime import timedelta
    from numpy import arange, isfinite

    method = method.lower()
    propgLag=[]
    epochLag=[]
    if method == 'flat':
     for i in range(len(srcData['epoch'])):
      if srcData['Vx'][i] != 0:
       timeLag = srcData['SCxGSE'][i]-destPos['X'][i]
       timeLag = timeLag/abs(srcData['Vx'][i])
       propgLag.extend([timeLag])
       epochLag.extend([epoch[i] + timedelta(0,propgLag[i])])
      elif srcData['Vx'][i] == 0 and i > 0:
       propgLag.extend([propgLag[i-1]])
       epochLag.extend([epoch[i] + timedelta(0,propgLag[i])])
      elif srcData['Vx'][i] == 0:
       print 'I found Vx = 0 at index = ', i
    return propgLag, epochLag


def kdeBW(obj, fac=1./5):
    from numpy import power
    """
       We use Scott's Rule, multiplied by a constant factor 
       to calculate the KDE Bandwidth.
    """
    return power(obj.n, -1./(obj.d+4)) * fac


def getIndices(inList, item):
    from numpy import isnan
    indices = [i for i in range(len(inList)) if isnan(inList[i])]
    return indices


def getDesKDE(srcSC,desSC,srcRanges,nPins=10):
    from scipy.stats import gaussian_kde
    from numpy import linspace,array
    desRanges = []
    desKDE = []
    nRanges = len(srcRanges)

    for i in range(nRanges):
     desRanges.append([])
     desKDE.append([])
     for j in range(len(srcSC)):
      if srcSC[j] >= srcRanges[i][0] and srcSC[j] <= srcRanges[i][1]:
       desRanges[i].extend([desSC[j]])
     desRanges[i]
    #desRanges[i] = dataClean(desRanges[i],[250],['<'])
    #desRanges[i] = removeNaN(desRanges[i])
     desRanges[i] = rejectOutliers(array(desRanges[i]))
     try:
      if len(desRanges[i]) > 1:
       jKDE = gaussian_kde(desRanges[i], bw_method=kdeBW)
       jVAL = linspace(min(desRanges[i]),max(desRanges[i]),nPins)
       desKDE[i].extend(jKDE(jVAL))
      else:
       desKDE[i] = []
     except:
      desKDE[i] = []
    return desRanges,desKDE
    '''
    for j in range(nRanges):
     for i in range(len(srcSC)):
      if srcSC[i] >= srcRanges[j][0] and srcSC[i] <= srcRanges[j][1]:
       desRanges[j].extend([desSC[i]])
     desRanges[j] = rejectOutliers(array(desRanges[j]))
     try:
      if len(desRanges[j]) > 1:
       jKDE = gaussian_kde(desRanges[j], bw_method=kdeBW)
       jVAL = linspace(min(desRanges[j]),max(desRanges[j]),nPins)
       desKDE[j].extend(jKDE(jVAL))
      else:
       desKDE[j] = []
     except:
      desKDE[j] = []
    return desRanges,desKDE
    '''

def swMedFilter(swEpoch,swParam,nSeconds):
    from scipy.signal import medfilt
    from getswdata import dateShift

    epochDiff = swEpoch[-1] - swEpoch[0]
    epochSize = epochDiff.days*(24*60*60) + epochDiff.seconds
    if epochSize < nSeconds:
     print 'Epoch size is too small'
     return ""

    eFilter = epochSize/nSeconds
    if eFilter%2 == 0: eFilter = eFilter + 1
    swParamMF = medfilt(swParam,eFilter)
 
    return swParamMF

def rejectOutliers(data,m=20.):
    from numpy import median
    d = abs(data - median(data))
    mdev = median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

def ccorr(x, y):
    from numpy.fft import fft, ifft
    from numpy     import argmax
    """Periodic correlation, implemented using the FFT.
       x and y must be real sequences with the same length.
    """
    xyccorr = ifft(fft(x) * fft(y).conj())
    xyccorr = xyccorr/max(abs(xyccorr))
    return xyccorr, argmax(abs(xyccorr))


def xcorr(x, y, method = 'pearsonr', shift = 5):
    from numpy import correlate, array, argmax, arange, linspace
    from scipy.stats import spearmanr, pearsonr
    method = method.lower()
    vCorr = []
    tCorr = []
    if shift == 0:
     if len(x) >= len(y):
      if method == 'pearsonr':
       v,t =  pearsonr(x[0:len(y)],y)
      elif method == 'spearmanr':
       v,t = spearmanr(x[0:len(y)],y)
     elif len(x) < len(y):
      if method == 'pearsonr':
       v,t =  pearsonr(x,y[0:len(x)])
      elif method == 'spearmanr':
       v,t = spearmanr(x,y[0:len(x)])
      vCorr.append(v)
    elif shift > 0 and shift < len(y):
     iCounter = 0
     for j in range(0,len(y)-shift,shift):
      vCorr.append([])
      tCorr.append([])
      for i in range(0,len(x)):
       if len(x[i:i+len(y[j:j+shift])]) == len(y[j:j+shift]):
        if method == 'pearsonr':
         v,t =  pearsonr(x[i:i+len(y[j:j+shift])],y[j:j+shift])
        elif method == 'spearmanr':
         v,t = spearmanr(x[i:i+len(y[j:j+shift])],y[j:j+shift])
        vCorr[iCounter].extend([v])
        tCorr[iCounter].extend([i])
      iCounter = iCounter + 1
        

    return array(vCorr), array(tCorr)
    

def search(a,val):
    ind = []
    for i in range(len(a)):
     if a[i] == val: ind = ind + [i]
    return ind

    
def normalize(inList):
    s = sum(inList)
    return map(lambda x: float(x)/s, inList)


#def xcorr(x, y, k, normalize=True):
#   import numpy as np
#   n = x.shape[0]

#   # initialize the output array
#   out = np.empty((2 * k) + 1, dtype=np.double)
#   lags = np.arange(-k, k + 1)

#   # pre-compute E(x), E(y)
#   mu_x = x.mean()
#   mu_y = y.mean()

#   # loop over lags
#   for ii, lag in enumerate(lags):

#       # use slice indexing to get 'shifted' views of the two input signals
#       if lag < 0:
#           xi = x[:lag]
#           yi = y[-lag:]
#       elif lag > 0:
#           xi = x[:-lag]
#           yi = y[lag:]
#       else:
#           xi = x
#           yi = y

#       # x - mu_x; y - mu_y
#       xdiff = xi - mu_x
#       ydiff = yi - mu_y

#       # E[(x - mu_x) * (y - mu_y)]
#       out[ii] = xdiff.dot(ydiff) / n

#       # NB: xdiff.dot(ydiff) == (xdiff * ydiff).sum()

#   if normalize:
#       # E[(x - mu_x) * (y - mu_y)] / (sigma_x * sigma_y)
#       out /=  np.std(x) * np.std(y)

#   return out, lags


