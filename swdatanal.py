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
    while cEpoch < epoch[-5]:
     sEpochID  = bisect_left(epoch, cEpoch + timedelta(0,0))
     eEpochID  = bisect_left(epoch, cEpoch + timedelta(0,(blockLen-1)*3600))
     TT,DD = removeNaN(epoch[sEpochID:eEpochID],data[sEpochID:eEpochID])
     if diff(TT) != []:
      dt = max(diff(TT))
      if dt.seconds <= gapLen:
       blockstart.append(cEpoch)
       cEpoch = cEpoch + timedelta(0,blockLen*3600)
      else:
       if fixStep:
        cEpoch = cEpoch + timedelta(0,blockLen*3600)
       else:
        cEpoch = cEpoch + timedelta(0,3600)

    return blockstart

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


def getIndices(inList, item):
    from numpy import isnan
    indices = [i for i in range(len(inList)) if isnan(inList[i])]
    return indices


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
   

def getSolarWindType(SWPList):
    from  numpy import array, log10, dot, sign
    import scipy.constants

    Tp = (scipy.constants.k*array(SWPList['T']))/scipy.constants.physical_constants['electron volt'][0]
    Sp = Tp/((array(SWPList['N']))**0.667)
    Tr = ((array(SWPList['V'])/258.0)**3.113)/Tp
    VA = 21.8*array(SWPList['B'])/((array(SWPList['N']))**0.5)

    dx = log10(Sp)
    dy = log10(VA)
    dz = log10(Tr)

    SWPClass = []
    for i in range(len(dx)):
     if dy[i] > 0.277 * dx[i] +0.055 * dz[i] + 1.83:           # Ejecta
      SWPClass = SWPClass + ['EJT']
     elif dx[i] > -0.525 * dz[i] - 0.676 * dy[i] + 1.74:       # Coronal-Hole_Origin
      SWPClass = SWPClass + ['CHO']
     elif dx[i] < -0.658 * dy[i] - 0.125 * dz[i] + 1.04:       # Sector-Reversal-Region
      SWPClass = SWPClass + ['SRR']
     else:                                                     # Streamer-Belt-Origin
      SWPClass = SWPClass + ['SBO']

    return list(Sp), list(Tr), list(VA), SWPClass


def getTimeLag(srcData,destPos,method='standard'):
    from math import atan2, tan, degrees, radians
    from datetime import timedelta
    from numpy import arange, isfinite

    method = method.lower()
    propgLag=[]
    epochLag=[]
    if method == 'standard':
     for i in range(len(srcData['epoch'])):
     #alpha   = atan2(srcData['Bx'][i],srcData['By'][i])
     #if (degrees(alpha) <= 89.9 and degrees(alpha) >= -89.9) or (degrees(alpha) > 91.1 and degrees(alpha) < 269.9):
     # timeLag = timeLag + (srcData['SCyGSE'][i] - destPos['Y'][i]) * tan(alpha)
     # print 'timeLag = ',timeLag, ', ACE_y = ', srcData['SCyGSE'][i], ', IMP_y = ', destPos['Y'][i], ', tan(B_x/B_y) = ', degrees(alpha)
     #print srcData['SCxGSE'][i], destPos['X'][i], 'timeLag = ', timeLag, 'Vx = ', srcData['Vx'][i]
      if srcData['Vx'][i] != 0:
       timeLag = srcData['SCxGSE'][i]-destPos['X'][i]
       timeLag = timeLag/abs(srcData['Vx'][i])
       propgLag.extend([timeLag])
       epochLag.extend([srcData['epoch'][i] + timedelta(0,propgLag[i])])
      elif srcData['Vx'][i] == 0 and i > 0:
       propgLag.extend([propgLag[i-1]])
       epochLag.extend([srcData['epoch'][i] + timedelta(0,propgLag[i])])
      elif srcData['Vx'][i] == 0:
       print i
    return propgLag, epochLag

def epochShift(usEpoch,usParam,sLag):
    from bisect import bisect_left
    from numpy import zeros
    from datetime import timedelta
    cEpoch = []
    for i in range(len(usEpoch)):
     cEpoch.extend([usEpoch[i] + timedelta(0,sLag[i])])
    if cEpoch[0] >= usEpoch[0]:
     sEpochID  = bisect_left(usEpoch, cEpoch[0])
     eEpochID  = bisect_left(usEpoch,usEpoch[-1])
    elif cEpoch[0] < usEpoch[0]:
     sEpochID  = bisect_left(usEpoch,usEpoch[0])
     eEpochID  = bisect_left(usEpoch, cEpoch[-1])

    sEpoch = usEpoch[sEpochID:eEpochID]
    sParam = usParam[sEpochID:eEpochID]

    return sEpoch,sParam


def kdeBW(obj, fac=1./5):
    from numpy import power
    """
       We use Scott's Rule, multiplied by a constant factor 
       to calculate the KDE Bandwidth.
    """
    return power(obj.n, -1./(obj.d+4)) * fac

def getDesKDE(srcSC,desSC,srcRanges,nPins=10):
    from scipy.stats import gaussian_kde
    from numpy import linspace
    desRanges = []
    desKDE = []
    nRanges = len(srcRanges)
    for i in range(nRanges):
     desRanges.append([])
     desKDE.append([])

    for j in range(nRanges):
     for i in range(len(srcSC)):
      if srcSC[i] >= srcRanges[j][0] and srcSC[i] <= srcRanges[j][1]:
       desRanges[j].extend([desSC[i]])
     if len(desRanges[j]) > 1:
      jKDE = gaussian_kde(desRanges[j], bw_method=kdeBW)
      jVAL = linspace(min(desRanges[j]),max(desRanges[j]),nPins)
      desKDE[j].extend(jKDE(jVAL))
     else:
      desKDE[j] = []
    return desRanges,desKDE

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


def ccorr(x, y):
    from numpy.fft import fft, ifft
    from numpy     import argmax
    """Periodic correlation, implemented using the FFT.
       x and y must be real sequences with the same length.
    """
    xyccorr = ifft(fft(x) * fft(y).conj())
    xyccorr = xyccorr/max(abs(xyccorr))
    return xyccorr, argmax(abs(xyccorr))


def xcorr(x, y, method = 'pearsonr'):
    from numpy import correlate, array, argmax, arange, linspace
    from scipy.stats import spearmanr, pearsonr
    method = method.lower()
    vCorr = []
    tCorr = []
    if len(x) > len(y):
     for i in range(len(x)-len(y)+1):
      if method == 'pearsonr':
       v,t = pearsonr(x[i:i+len(y)],y)
      elif method == 'spearmanr':
       v,t = spearmanr(x[i:i+len(y)],y)
      vCorr.append(v)
      tCorr.append(i)
    elif len(x) == len(y):
     if method == 'pearsonr':
      v,t = pearsonr(x,y)
     elif method == 'spearmanr':
      v,t = spearmanr(x,y)
     vCorr.append(v)
     tCorr.append(i)
    elif len(x) < len(y):
     print 'length of x is smaller than length of y is not allowed'

    return vCorr, tCorr
    

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


