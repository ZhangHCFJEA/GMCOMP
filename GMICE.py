#!/usr/bin/python
import math
import numpy
#############################################
#Conversions from PGA or PGV to MMI using
#the GMICE equations from Wald et al. [1999] and
#Worden et al. [2012] - WGRW12
#Relationship between peak ground acceleration,
#peak ground velocity, and Modified Mercalli
#Intensity in California
#Brendan Crowell, University of Washington
#Last Updated October 16, 2018
#The PGA and PGV values should be in a numpy array
#vector
#############################################
def MMI_Y(y, mode, props):
    gmice = props.getgmice()

    if (gmice == 'wald99'):
        MMI = Wald99(y,mode)

    if (gmice == 'wgrw12'):
        MMI = WGRW12(y,mode)

    return(MMI)



def Wald99(y, mode):
    if (mode == 0):
        pgalen = len(y)
        MMI = numpy.zeros((pgalen,1))
        pgalarge = numpy.where(numpy.log10(y) >= 1.819)
        pgasmall = numpy.where((numpy.log10(y) < 1.819) & (y > 0))
        pgazero = numpy.where(y == 0)
        MMI[pgalarge] = 3.66*numpy.log10(y[pgalarge])-1.66
        MMI[pgasmall] = 2.20*numpy.log10(y[pgasmall])+1.00
        MMI[pgazero] = -10
    else:
        pgvlen = len(y)
        MMI = numpy.zeros((pgvlen,1))
        pgvlarge = numpy.where(numpy.log10(y) >= 0.763)
        pgvsmall = numpy.where(numpy.log10(y) < 0.763)
        MMI[pgvlarge] = 3.47*numpy.log10(y[pgvlarge])+2.35
        MMI[pgvsmall] = 2.10*numpy.log10(y[pgvsmall])+3.4
        
    return(MMI)


def WGRW12(y, mode):
    if (mode == 0):
        pgalen = len(y)
        MMI = numpy.zeros((pgalen,1))
        pgalarge = numpy.where(numpy.log10(y) >= 1.57)
        pgasmall = numpy.where((numpy.log10(y) < 1.57) & (y > 0))
        pgazero = numpy.where(y == 0)
        
        MMI[pgalarge] = 3.70*numpy.log10(y[pgalarge])-1.60
        MMI[pgasmall] = 1.55*numpy.log10(y[pgasmall])+1.78
        MMI[pgazero] = -10
    else:
        pgvlen = len(y)
        MMI = numpy.zeros((pgvlen,1))
        pgvlarge = numpy.where(numpy.log10(y) >= 0.53)
        pgvsmall = numpy.where(numpy.log10(y) < 0.53)
        MMI[pgvlarge] = 3.16*numpy.log10(y[pgvlarge])+2.89
        MMI[pgvsmall] = 1.47*numpy.log10(y[pgvsmall])+3.78
        
    return(MMI)

