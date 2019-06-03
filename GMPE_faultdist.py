#!/usr/bin/python
import math
import GMCOMP_coord_tools
import numpy
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

#This code computes the approximate rupture distance(Rrup),
#the Joyner Boore distance(Rjb), Rx, and hanging wall flag
#Written by Brendan Crowell, University of Washington, 2017
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
def sourcedistance(lats,lons,deps,ss,ds,latsta,lonsta,site_length,nvert,props):
    p = float(props.getpvalue())
    [x,y]=GMCOMP_coord_tools.ll2utm(lons, lats, numpy.median(lons), numpy.median(lats))
    z = deps

    rake = 180/math.pi*numpy.arctan2(ds,ss)

    fault_slip = numpy.sqrt(numpy.power(ss,2)+numpy.power(ds,2))

    totslip = numpy.sum(fault_slip)

    xravel = numpy.ravel(x)
    yravel = numpy.ravel(y)
    zravel = numpy.ravel(z)
    

    if (nvert == 4):
        Z1 = z[:,0]
        Z2 = z[:,1]
        Z3 = z[:,2]
        Z4 = z[:,3]

        X1 = x[:,0]
        X2 = x[:,1]
        X3 = x[:,2]
        X4 = x[:,3]
        
        Y1 = y[:,0]
        Y2 = y[:,1] 
        Y3 = y[:,2]   
        Y4 = y[:,3]

        fseg = len(x)

##        for i in range (0, fseg):
##            XS = [X1[i], X2[i], X3[i]]
##            YS = [Y1[i], Y2[i], Y3[i]]
##            ZS = [Z1[i], Z2[i], Z3[i]]
##
##            [strike, dip] = strikedip(XS, YS, ZS)

            
        Ztor = numpy.amin(numpy.vstack((Z1,Z2,Z3,Z4)))/1000



        atop = numpy.where(numpy.amin(zravel)  == zravel)[0]
        abot = numpy.where(numpy.amax(zravel) == zravel)[0]
        yravelt = numpy.ravel(yravel[atop])
        xravelt = numpy.ravel(xravel[atop])
        zravelt = numpy.ravel(zravel[atop])
        atops = numpy.where((numpy.amin(yravelt)  == yravelt))[0]
        atopn = numpy.where((numpy.amax(yravelt) == yravelt))[0]
        yravelb = numpy.ravel(yravel[abot])
        xravelb = numpy.ravel(xravel[abot])
        zravelb = numpy.ravel(zravel[abot])
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
        
        XS = [xtops[0], xtopn[0], xbotn[0]]
        YS = [ytops[0], ytopn[0], ybotn[0]]
        ZS = [ztops[0], ztopn[0], zbotn[0]]
        [strike, dip] = strikedip(XS, YS, ZS)
        

        Rp = numpy.zeros([site_length,1])
        Rrup = numpy.zeros([site_length,1])
        Rjb = numpy.zeros([site_length,1])
        Rx = numpy.zeros([site_length,1])
        Hw = numpy.zeros([site_length,1])
        
        for i in range (0, site_length):
            [xs,ys] = GMCOMP_coord_tools.ll2utm(lonsta[i,0],latsta[i,0],numpy.median(lons),numpy.median(lats))
            R1 = numpy.sqrt(numpy.power(xs-X1,2)+numpy.power(ys-Y1,2)+numpy.power(Z1,2))
            R2 = numpy.sqrt(numpy.power(xs-X2,2)+numpy.power(ys-Y2,2)+numpy.power(Z2,2))
            R3 = numpy.sqrt(numpy.power(xs-X3,2)+numpy.power(ys-Y3,2)+numpy.power(Z3,2))
            R4 = numpy.sqrt(numpy.power(xs-X4,2)+numpy.power(ys-Y4,2)+numpy.power(Z4,2))
            R = numpy.vstack((R1,R2,R3,R4))
            dist = numpy.mean(R, axis=0)/1000

            Rrup[i,0] = numpy.amin(dist)
            
            Rp[i,0] = numpy.power(numpy.sum(fault_slip*numpy.power(dist,-1.4)),-1/1.4)/numpy.power(numpy.sum(fault_slip),-1/1.4)

            xl1 = xtops[0]
            xl2 = xtopn[0]
            yl1 = ytops[0]
            yl2 = ytopn[0]
               
            Rx[i,0] = perp_dist(xl1, yl1, xl2, yl2, xs, ys)/1000

            if (dip < 90):                          
                point = Point(xs,ys)               
                polygon = Polygon([(xl1,yl1), (xl1,yl2), (xl2,yl2), (xl2,yl1)])
                Rjb[i,0] = point.distance(polygon)/1000.0
            else:
                Rjb[i,0] = Rx[i,0]
                    
            Hw[i,0] = numpy.sign((yl2-yl1)*xs-(xl2-xl1)*ys+xl2*yl1-yl2*xl1)



##        a1 = numpy.where(numpy.amin(z) == z)[0]
##        #Southern and northern extent of top of fault
##        atops = numpy.where((numpy.amin(yravel) == yravel) & (numpy.amin(zravel) == zravel))[0]
##        atopn = numpy.where((numpy.amax(yravel) == yravel) & (numpy.amin(zravel) == zravel))[0]
##        abots = numpy.where((numpy.amin(yravel) == yravel) & (numpy.amax(zravel) == zravel))[0]
##        abotn = numpy.where((numpy.amax(yravel) == yravel) & (numpy.amax(zravel) == zravel))[0]
##        print (atops,atopn,abots,abotn)
##
##        #If fault is true E-W, ends defined by east west extent
##        if (strike == 90 or strike == 270):
##            atops = numpy.where((numpy.amax(xravel) == xravel) & (numpy.amin(zravel) == zravel))[0]
##            atopn = numpy.where((numpy.amin(xravel) == xravel) & (numpy.amin(zravel) == zravel))[0]
##            abots = numpy.where((numpy.amin(xravel) == xravel) & (numpy.amax(zravel) == zravel))[0]
##            abotn = numpy.where((numpy.amax(xravel) == xravel) & (numpy.amax(zravel) == zravel))[0]
##
##        FLen = math.sqrt(math.pow(xravel[atops[0]]-xravel[atopn[0]],2)+math.pow(yravel[atops[0]]-yravel[atopn[0]],2))/1000
##        FWid = math.sqrt(math.pow(xravel[atops[0]]-xravel[abots[0]],2)+math.pow(yravel[atops[0]]-yravel[abots[0]],2)+math.pow(zravel[atops[0]]-zravel[abots[0]],2))/1000
##
##        Rrup = numpy.zeros([site_length,1])
##        Rjb = numpy.zeros([site_length,1])
##        Rx = numpy.zeros([site_length,1])
##        Rp = numpy.zeros([site_length,1])
##        Hw = numpy.zeros([site_length,1])
##        for i in range (0, site_length):
##            [xs,ys] = GMCOMP_coord_tools.ll2utm(lonsta[i,0],latsta[i,0],numpy.median(lons),numpy.median(lats))
##
##            R1 = numpy.sqrt(numpy.power(xs-X1,2)+numpy.power(ys-Y1,2)+numpy.power(Z1,2))
##            R2 = numpy.sqrt(numpy.power(xs-X2,2)+numpy.power(ys-Y2,2)+numpy.power(Z2,2))
##            R3 = numpy.sqrt(numpy.power(xs-X3,2)+numpy.power(ys-Y3,2)+numpy.power(Z3,2))
##            R4 = numpy.sqrt(numpy.power(xs-X4,2)+numpy.power(ys-Y4,2)+numpy.power(Z4,2))
##            R = numpy.vstack((R1,R2,R3,R4))
##            dist = numpy.mean(R, axis=0)/1000
##
##            Rrup[i,0] = numpy.amin(dist)
##
##            Rp[i,0] = numpy.power(numpy.sum(fault_slip/totslip*numpy.power(dist,p)),1/p)
##            xmin = numpy.min(x)
##            xmax = numpy.max(x)
##            ymin = numpy.min(y)
##            ymax = numpy.max(y)
##
##
##            xsoff = 1000*math.sin(strike*math.pi/180)
##            ysoff = 1000*math.cos(strike*math.pi/180)
##
##            xl1 = xravel[atops[0]]
##            xl2 = xravel[atopn[0]]
##            yl1 = yravel[atops[0]]
##            yl2 = yravel[atopn[0]]
##           
##            Rx[i,0] = perp_dist(xl1, yl1, xl2, yl2, xs, ys)/1000
##            
##            if (dip < 90):                          
##                point = Point(xs,ys)               
##                polygon = Polygon([(xmin,ymin), (xmin,ymax), (xmax,ymax), (xmax,ymin)])
##                Rjb[i,0] = point.distance(polygon)/1000.0
##            else:
##                Rjb[i,0] = Rx[i,0]
##
##        
##            Hw[i,0] = numpy.sign((yl2-yl1)*xs-(xl2-xl1)*ys+xl2*yl1-yl2*xl1)
        

    return(Rrup,Rjb,Rx,Rp,Ztor,Hw,rake,strike,dip,FLen,FWid)


def strikedip(x,y,z):
    x1 = x[0]
    x2 = x[1]
    x3 = x[2]

    y1 = y[0]
    y2 = y[1]
    y3 = y[2]

    z1 = -z[0]
    z2 = -z[1]
    z3 = -z[2]

    u1 = (y1-y2)*(z3-z2) - (y3-y2)*(z1-z2)
    u2 = -((x1-x2)*(z3-z2) - (x3-x2)*(z1-z2))
    u3 = (x1-x2)*(y3-y2) - (x3-x2)*(y1-y2)



    if (u3 >= 0):
        northing = u1
        easting = -u2
    else:
        northing = -u1
        easting = u2

    if (easting > 0):
        strike = math.degrees(math.acos(northing/math.sqrt(math.pow(easting,2)+math.pow(northing,2))))
    else:
        strike = math.degrees(2*math.pi - math.acos(northing/math.sqrt(math.pow(easting,2)+math.pow(northing,2))))

    dip = math.degrees(math.asin(math.sqrt(math.pow(u1,2)+math.pow(u2,2))/math.sqrt(math.pow(u1,2)+math.pow(u2,2)+math.pow(u3,2))))
        

    return (strike, dip)


def perp_dist(x1, y1, x2, y2, x3, y3): # x3,y3 is the point
    px = x2-x1
    py = y2-y1

    norm = px*px + py*py

    u =  ((x3 - x1) * px + (y3 - y1) * py) / float(norm)

    if u > 1:
        u = 1
    elif u < 0:
        u = 0

    x = x1 + u * px
    y = y1 + u * py

    dx = x - x3
    dy = y - y3

    dist = (dx*dx + dy*dy)**.5

    return (dist)


def sourcedistance_point(lateq,loneq,depeq,latsta,lonsta,site_length,props):
    [x,y]=GMCOMP_coord_tools.ll2utm(loneq, lateq, loneq, lateq)
    z = depeq*1000
 
    Fwid = 0
    FLen = 0

    Ztor = depeq


    Rrup = numpy.zeros([site_length,1])
    Rjb = numpy.zeros([site_length,1])
    Rx = numpy.zeros([site_length,1])
    Rp = numpy.zeros([site_length,1])
    Hw = numpy.zeros([site_length,1])
    for i in range (0, site_length):
        [xs,ys] = GMCOMP_coord_tools.ll2utm(lonsta[i,0],latsta[i,0],loneq,lateq)

        Rhyp = numpy.sqrt(numpy.power(xs-x,2)+numpy.power(ys-y,2)+numpy.power(z,2))/1000
        Repi = numpy.sqrt(numpy.power(xs-x,2)+numpy.power(ys-y,2))/1000


        Rrup[i,0] = Rhyp
        Rp[i,0] = Rhyp      
        Rx[i,0] = Repi
        Rjb[i,0] = Repi
        Hw[i,0] = 0
        

    return(Rrup,Rjb,Rx,Rp,Ztor,Hw,Fwid)


def sourcedistance_line(latff,lonff,depff,latsta,lonsta,site_length,props):
    loneq = numpy.mean(lonff)
    lateq = numpy.mean(latff)
    Ztor = numpy.mean(depff)
    
    [x,y]=GMCOMP_coord_tools.ll2utm(lonff, latff, loneq, lateq)
 
    Fwid = 0
    FLen = 0


    Rrup = numpy.zeros([site_length,1])
    Rjb = numpy.zeros([site_length,1])
    Rx = numpy.zeros([site_length,1])
    Rp = numpy.zeros([site_length,1])
    Hw = numpy.zeros([site_length,1])
    
    
    for i in range (0, site_length):
        rx = list()
        for j in range (1, len(latff)):
            [xs,ys] = GMCOMP_coord_tools.ll2utm(lonsta[i,0],latsta[i,0],loneq,lateq)

            rxdum = perp_dist(x[j], y[j], x[j-1], y[j-1], xs, ys)/1000
            rx.append(rxdum)

        Rx[i,0] = numpy.amin(numpy.asarray(rx))

        Rrup[i,0] = Rx[i,0]
        Rp[i,0] = Rx[i,0]     
        Rjb[i,0] = Rx[i,0]
        

    return(Rrup,Rjb,Rx,Rp,Ztor,Hw,Fwid)
