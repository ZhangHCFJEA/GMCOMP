#!/usr/bin/python
import math
import GMCOMP_coord_tools
import numpy
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from GMPE_faultdist import strikedip, perp_dist

#This code takes in a srcmod formatted fsp file and computes the
#mean rupture distance to your stations
#Mean rupture distance is defined in Thompson and Baltay [2018]
#Written by Brendan Crowell, University of Washington, 2018
#######################################################
#Variable Definitions
#######################################################
#latc is the latitude of the center of the fault
#lonc is the longitude of the center of the fault
#depc is the depth of the hypocenter
#LEN is the length of the fault plane in km
#WID is the width of the fault plane in km
#strike is the strike of the fault plane
#dip is the dip of the fault plane
#latsta are the latitudes of the stations (format as numpy array)
#lonsta are the longitudes of the stations (format as numpy array)
#######################################################

def srcmod_fspreader(fspfile,sta_lat,sta_lon):
    fault_lat = []
    fault_lon = []
    fault_dep = []
    fault_slip = []
    fault_rake = []

    
    for line in open(fspfile):
        li=line.strip()
        if not li.startswith("%"):
            fs = line.split()
            if not fs:
                pass
            else:
                fault_lat = numpy.append(fault_lat,[float(fs[0])],axis=0)
                fault_lon = numpy.append(fault_lon,[float(fs[1])],axis=0)
                fault_dep = numpy.append(fault_dep,[float(fs[4])],axis=0)
                fault_slip = numpy.append(fault_slip,[float(fs[5])],axis=0)
                fault_rake = numpy.append(fault_rake,[float(fs[7])],axis=0)

    lat0 = numpy.mean(fault_lat)
    lon0 = numpy.mean(fault_lon)
    totslip = numpy.sum(fault_slip)
    rake = numpy.mean(fault_rake)

    [x0,y0]=GMCOMP_coord_tools.ll2utm(fault_lon, fault_lat, lon0, lat0) #convert fault locations to UTM
    [xs,ys]=GMCOMP_coord_tools.ll2utm(sta_lon, sta_lat, lon0, lat0) #convert station locations to UTM
    
    Ztor = numpy.amin(fault_dep)

    atop = numpy.where(Ztor == fault_dep)[0]
    abot = numpy.where(numpy.amax(fault_dep) == fault_dep)[0]
    yravelt = numpy.ravel(y0[atop])
    xravelt = numpy.ravel(x0[atop])
    zravelt = numpy.ravel(fault_dep[atop])
    atops = numpy.where((numpy.amin(yravelt)  == yravelt))[0]
    atopn = numpy.where((numpy.amax(yravelt) == yravelt))[0]
    yravelb = numpy.ravel(y0[abot])
    xravelb = numpy.ravel(x0[abot])
    zravelb = numpy.ravel(fault_dep[abot])
    abots = numpy.where((numpy.amin(yravelb) == yravelb))[0]
    abotn = numpy.where((numpy.amax(yravelb) == yravelb))[0]



    ytopn = yravelt[atopn]
    ytops = yravelt[atops]
    ybotn = yravelb[abotn]
    ybots = yravelb[abots]

    xtopn = xravelt[atopn]
    xtops = xravelt[atops]
    xbotn = xravelb[abotn]
    xbots = xravelb[abots]

    ztopn = zravelt[atopn]
    ztops = zravelt[atops]
    zbotn = zravelb[abotn]
    zbots = zravelb[abots]

    FLen = math.sqrt(math.pow(xtops[0]-xtopn[0],2)+math.pow(ytops[0]-ytopn[0],2))/1000
    FWid = math.sqrt(math.pow(xtops[0]-xbots[0],2)+math.pow(ytops[0]-ybots[0],2)+math.pow((ztops[0]-zbots[0])*1000,2))/1000
    
    XS = [float(xtops[0]), float(xtopn[0]), float(xbotn[0])]
    YS = [float(ytops[0]), float(ytopn[0]), float(ybotn[0])]
    ZS = [1000*float(ztops[0]), 1000*float(ztopn[0]), 1000*float(zbotn[0])]
    [strike, dip] = strikedip(XS, YS, ZS)
    fault_strike = strike*numpy.ones([len(fault_lat),1])
    fault_dip = dip*numpy.ones([len(fault_lat),1])
    
    numsta = len(xs) #number of stations
    Rp = numpy.zeros([numsta,1])
    Rrup = numpy.zeros([numsta,1])
    Rjb = numpy.zeros([numsta,1])
    Rx = numpy.zeros([numsta,1])
    Hw = numpy.zeros([numsta,1])
    
    for i in range (0, numsta):
        dist = numpy.sqrt(numpy.power(x0-xs[i],2)+numpy.power(y0-ys[i],2)+numpy.power(fault_dep*1000,2))/1000
        Rrup[i,0] = numpy.amin(dist)
        
        Rp[i,0] = numpy.power(numpy.sum(fault_slip*numpy.power(dist,-1.4)),-1/1.4)/numpy.power(numpy.sum(fault_slip),-1/1.4)

        xl1 = xtops[0]
        xl2 = xtopn[0]
        yl1 = ytops[0]
        yl2 = ytopn[0]
           
        Rx[i,0] = perp_dist(xl1, yl1, xl2, yl2, xs[i], ys[i])/1000

        if (dip < 90):                          
            point = Point(xs[i],ys[i])               
            polygon = Polygon([(xl1,yl1), (xl1,yl2), (xl2,yl2), (xl2,yl1)])
            Rjb[i,0] = point.distance(polygon)/1000.0
        else:
            Rjb[i,0] = Rx[i,0]
                
        Hw[i,0] = numpy.sign((yl2-yl1)*xs[i]-(xl2-xl1)*ys[i]+xl2*yl1-yl2*xl1)


    return(Rrup,Rjb,Rx,Rp,Ztor,Hw,rake,strike,dip,FLen,FWid,fault_lat,fault_lon,fault_dep,fault_slip,fault_rake,fault_strike,fault_dip)
        

    
