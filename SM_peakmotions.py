#!/usr/bin/python
import math
import numpy
import obspy
from GMICE import MMI_Y

def pga_measure(file_east,file_north,file_vert,gain_east,gain_north,gain_vert,origintime,props):
    mode=0
    ttot = int(props.getttot())
    east = obspy.read(file_east)
    north = obspy.read(file_north)
    vert = obspy.read(file_vert)

    te_start = east[0].stats.starttime
    te_end = east[0].stats.endtime

    tn_start = north[0].stats.starttime
    tn_end = north[0].stats.endtime

    tv_start = vert[0].stats.starttime
    tv_end = vert[0].stats.endtime

    st = numpy.array([[te_start],[tn_start],[tv_start]])
    stmin = numpy.amin(st)
    et = numpy.array([[te_end],[tn_end],[tv_end]])
    etmax = numpy.amax(et)

    if (etmax < origintime+ttot):
        etime = etmax
    else:
        etime = origintime+ttot


    etrim = east[0].trim(starttime=origintime,endtime=etime,pad=True, fill_value=0)
    ntrim = north[0].trim(starttime=origintime,endtime=etime,pad=True, fill_value=0)
    vtrim = vert[0].trim(starttime=origintime,endtime=etime,pad=True, fill_value=0)

    time = etrim.times()

    data_east = etrim.data/float(gain_east)
    data_north = ntrim.data/float(gain_north)
    data_vert = vtrim.data/float(gain_vert)

    data_east = data_east - data_east[0]
    data_north = data_north - data_north[0]
    data_vert = data_vert - data_vert[0]
    
    eacc = numpy.sqrt(numpy.power(data_east,2))
    nacc = numpy.sqrt(numpy.power(data_north,2))
    vacc = numpy.sqrt(numpy.power(data_vert,2))
    hacc = numpy.sqrt(numpy.power(data_east,2)+numpy.power(data_north,2))
    tacc = numpy.sqrt(numpy.power(data_east,2)+numpy.power(data_north,2)+numpy.power(data_vert,2))

    component = props.getmmicomp()

    if (component == 'e'):
        PGA = numpy.maximum.accumulate(eacc)
        MMI = MMI_Y(numpy.reshape(PGA,(len(PGA),1)), mode, props)
        PGAm = numpy.amax(PGA)
        MMIm = numpy.amax(MMI)
    if (component == 'n'):
        PGA = numpy.maximum.accumulate(nacc)
        MMI = MMI_Y(numpy.reshape(PGA,(len(PGA),1)), mode, props)
        PGAm = numpy.amax(PGA)
        MMIm = numpy.amax(MMI)
    if (component == 'v'):
        PGA = numpy.maximum.accumulate(vacc)
        MMI = MMI_Y(numpy.reshape(PGA,(len(PGA),1)), mode, props)
        PGAm = numpy.amax(PGA)
        MMIm = numpy.amax(MMI)
    if (component == 'h'):
        PGA = numpy.maximum.accumulate(hacc)
        MMI = MMI_Y(numpy.reshape(PGA,(len(PGA),1)), mode, props)
        PGAm = numpy.amax(PGA)
        MMIm = numpy.amax(MMI)
    if (component == 't'):
        PGA = numpy.maximum.accumulate(tacc)
        MMI = MMI_Y(numpy.reshape(PGA,(len(PGA),1)), mode, props)
        PGAm = numpy.amax(PGA)
        MMIm = numpy.amax(MMI)

    if (component == 'm'):
        pgae = numpy.maximum.accumulate(eacc)
        pgae = numpy.reshape(pgae,(len(pgae),1))
        
        pgan = numpy.maximum.accumulate(nacc)
        pgan = numpy.reshape(pgan,(len(pgan),1))
        
        pgav = numpy.maximum.accumulate(vacc)
        pgav = numpy.reshape(pgav,(len(pgav),1))
        
        pgaconcat = numpy.concatenate((pgae,pgan,pgav), axis=1)
        PGA = numpy.amax(pgaconcat,axis=1)
        MMI = MMI_Y(numpy.reshape(PGA,(len(PGA),1)), mode, props)
        PGAm = numpy.amax(PGA)
        MMIm = numpy.amax(MMI)
        


    return(time,PGA,MMI,PGAm,MMIm)
        

    
