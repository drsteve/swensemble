
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


def omniDataCorr(srefDate, erefDate, startDate, endDate, epochs, SWP, binStride, CorrTime = 'Day'):
    import numpy
    import bisect
    import datetime
    from scipy.stats import ks_2samp
    from getswdata import getOMNIfiles, dataClean, dateShift, dateList

    CorrTime = CorrTime.lower()

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
    SWPD01   = getDistrib(filter(lambda v: v==v, SWPV01), bins=bins, norm=False)

    if CorrTime == 'day':
     aepochs = []; KSVals = []; KSDist = []
     sEpoch = datetime.datetime(startDate.year,startDate.month,startDate.day, 0, 0, 0)
     eEpoch = dateShift(sEpoch, hours = 23, minutes = 59, seconds = 59)
     for i in range((endDate-startDate).days+1):
      aepochs  = aepochs + [dateShift(sEpoch,0,0,i,0,0,0)]
      sEpochID = bisect.bisect_left(epochs, dateShift(sEpoch,0,0,i,0,0,0))
      eEpochID = bisect.bisect_left(epochs, dateShift(eEpoch,0,0,i,0,0,0))

      SWPV02 = SWP[sEpochID:eEpochID]
      SWPD02 = getDistrib(filter(lambda v: v==v, SWPV02), bins=bins, norm=False)

      KSVals = KSVals + [ks_2samp(SWPV01, SWPV02)]
      KSDist = KSDist + [ks_2samp(SWPD01, SWPD02)]

     KSVals = numpy.array(KSVals)
     KSDist = numpy.array(KSDist)

    return SWPDatRng, cepochs, KSVals, KSDist, aepochs
   



 
