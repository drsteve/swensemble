import sys
import math
import numpy
import datetime
from scipy import signal
import matplotlib.pyplot as plt
from swdatanal import xcorr, ccorr, getDesKDE, search
from getswdata import commonEpoch, dateList

startDate = datetime.datetime(1998, 1, 1, 0, 0, 0)
endDate   = datetime.datetime(1998, 1,31,23,59,59)
epoch1    = dateList(startDate, endDate, shift = 'hour')

startDate = datetime.datetime(1998, 1,21, 0, 0, 0)
endDate   = datetime.datetime(1998, 1,21,23,59,59)
epoch2    = dateList(startDate, endDate, shift = 'hour')

print commonEpoch(epoch1,epoch2)
