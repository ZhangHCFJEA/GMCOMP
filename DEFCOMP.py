#!/usr/bin/python
import math
import numpy
import SLIP_okadagreen
import GMCOMP_coord_tools
import GMPE_faultdist

def def_grid(lonc, latc, dlon, dlat, ddeg, flon, flat, fdep, fss, fds, strike, dip, W, L, props):
        lons = numpy.linspace(lonc-dlon/2,lonc+dlon/2,num=ddeg)
        lats = numpy.linspace(latc-dlat/2,latc+dlat/2,num=ddeg)
        

        ngp = len(lons)*len(lats) #total number of grid points

        LONS = numpy.zeros([ngp,1])
        LATS = numpy.zeros([ngp,1])

        n = 0
        for i in range (0, len(lons)):
                for j in range (0, len(lats)):
                        LONS[n,0] = lons[i]
                        LATS[n,0] = lats[j]
                        n=n+1
        
        nfs = len(flon)
        
        xrs = numpy.zeros([nfs,ngp])
        yrs = numpy.zeros([nfs,ngp])
        zrs = numpy.zeros([nfs,ngp])



        for i in range (0, nfs):
                for j in range (0, ngp):
                        (x1,y1) = GMCOMP_coord_tools.ll2utm(LONS[j], LATS[j], numpy.median(flon), numpy.median(flat))
                        (x2,y2) = GMCOMP_coord_tools.ll2utm(flon[i], flat[i], numpy.median(flon), numpy.median(flat))
                        xrs[i,j] = (x1-x2)
                        yrs[i,j] = (y1-y2)
                        zrs[i,j] = fdep[i]
        print (fss, fds)

        G = SLIP_okadagreen.greenF(xrs, yrs, zrs, strike, dip, W, L) #Compute Green's functions

        S = numpy.zeros([2*nfs,1])

        for i in range (0, nfs):
                S[2*i,0] = fss[i]
                S[2*i+1,0] = fds[i]

        U = numpy.dot(G,S)

        EN = numpy.zeros([ngp,1])
        NN = numpy.zeros([ngp,1])
        UN = numpy.zeros([ngp,1])
        for i in range (0, ngp):
                EN[i,0]=U[3*i,0]
                NN[i,0]=U[3*i+1,0]
                UN[i,0]=U[3*i+2,0]

        return(LONS, LATS, EN, NN, UN)

def fp_organize(LONFF, LATFF, DEPFF):
        [x,y]=GMCOMP_coord_tools.ll2utm(LONFF, LATFF, numpy.median(LONFF), numpy.median(LATFF))
        z = DEPFF

        xravel = numpy.ravel(x)
        yravel = numpy.ravel(y)
        zravel = numpy.ravel(z)

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

        STR = numpy.zeros([len(X1),1])
        DIP = numpy.zeros([len(X1),1])
        L = numpy.zeros([len(X1),1])
        W = numpy.zeros([len(X1),1])

        XM = numpy.zeros([len(X1),1])
        YM = numpy.zeros([len(X1),1])
        ZM = numpy.zeros([len(X1),1])
        

        for i in range (0, len(X1)):
                XS = [X1[i], X2[i], X3[i]]
                YS = [Y1[i], Y2[i], Y3[i]]
                ZS = [Z1[i], Z2[i], Z3[i]]
                [st, di] = GMPE_faultdist.strikedip(XS, YS, ZS)
                STR[i,0] = st
                DIP[i,0] = di
                
                XS = [X1[i], X2[i], X3[i], X4[i]]
                YS = [Y1[i], Y2[i], Y3[i], Y4[i]]
                ZS = [Z1[i], Z2[i], Z3[i], Z4[i]]
                atop = numpy.where(numpy.amin(ZS)  == ZS)[0]
                abot = numpy.where(numpy.amax(ZS) == ZS)[0]

                L[i,0] = math.sqrt(math.pow(XS[atop[0]]-XS[atop[1]],2)+math.pow(YS[atop[0]]-YS[atop[1]],2))

                dH = numpy.amax(ZS)-numpy.amin(ZS)

                W[i,0] = dH/math.sin(DIP[i,0]*math.pi/180)

                XM[i,0] = numpy.mean(XS)
                YM[i,0] = numpy.mean(YS)
                ZM[i,0] = numpy.mean(ZS)

        [LONM, LATM]= GMCOMP_coord_tools.utm2ll(XM, YM, numpy.median(LONFF), numpy.median(LATFF))

                            

                

        return(LONM, LATM, ZM, STR, DIP, L, W)
        

        

	

