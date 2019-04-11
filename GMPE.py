#!/usr/bin/python
import math
import numpy
##############################################################
#GMPE.py - A concatenated version of individual GMPE's
#Written by Brendan W. Crowell, University of Washington
#Last Edited March 15, 2019
##############################################################
#Variables
#M = magnitude, single value
#Rrup = rupture distance, km, array
#Rjb = Joyner-Boore distance, km, array
#Rx = distance to surface projection of updip edge of slip, km, array
#Ztor = depth to top of rupture plane, km, array
#Zhyp = hypocentral depth, km, single value
#VS30 = shear wave velocity of upper 30 meters at a site, km/s, array
#Z25 = depth to 2.5 km/s shear wave velocity, in km, array
#Hw = is the hanging wall flag, 1 for true, array
#dip = dip of fault plane, in degrees, single value
#W = width of fault plane, in km, single value
#rake = rake angle of earthquake
##############################################################
def Y(M,Rrup,Rjb,Rx,Ztor,Zhyp,VS30,Hw,W,dip,rake,mode,props):
    gmpe = props.getgmpe()
    ######################################################
    #Chiou and Youngs 2008 - CY08
    ######################################################
    if (gmpe == 'cy08'):
        y = gmpe_cy08(M,Rrup,Rjb,Rx,Ztor,VS30,Hw,dip,rake,mode)
        if (mode == 0):
            y = y*980.665
    ######################################################
    ######################################################
    #Campbell and Bozorgnia 2014, NGA-West2
    ######################################################
    if (gmpe == 'cb14'):
        y = gmpe_cb14(M,Rjb,Rrup,Rx,Ztor,Zhyp,VS30,Hw,W,dip,rake,mode)
        if (mode == 0):
            y = y*980.665
    ######################################################
    ######################################################
    #Boore Atkinson 2008
    ######################################################
    if (gmpe == 'ba08'):
        y = gmpe_ba08(M,Rjb,VS30,rake,mode)
        if (mode == 0):
            y = y*980.665
    ######################################################
    ######################################################
    #Bradley 2013 - For New Zealand
    ######################################################
    if (gmpe == 'bradley13'):
        y = gmpe_bradley13(M,Rrup,Rjb,Rx,Rtvz,Ztor,VS30,Z10,Hw,dip)



    return(y)



def gmpe_cy08(M,Rrup,Rjb,Rx,Ztor,VS30,Hw,dip,rake,mode):
    ######################################################
    #Coefficients
    ######################################################
    #period independent coefficients
    c2=1.06
    c3=3.45
    c4=-2.1
    c4a=-0.5
    crb=50.0
    chm=3.0
    cg3=4.0
    if mode == 0:
        #PGA coefficients
        c1=-1.2687
        c1a=0.1
        c1b=-0.2550
        cn=2.996
        cm=4.1840
        c5=6.1600
        c6=0.4893
        c7=0.0512
        c7a=0.0860
        c9=0.7900
        c9a=1.5005
        c10=-0.3218
        cg1=-0.00804
        cg2=-0.00785
        phi1=-0.4417
        phi2=-0.1417
        phi3=-0.007010
        phi4=0.102151
        phi5=0.2289
        phi6=0.014996
        phi7=580.0
        phi8=0.0700
    if mode == 1:
        #pgv coefficients
        c1=2.2884
        c1a=0.1094
        c1b=-0.0626
        cn=1.648
        cm=4.2979
        c5=5.1700
        c6=0.4407
        c7=0.0207
        c7a=0.0437
        c9=0.3079
        c9a=2.6690
        c10=-0.1166
        cg1=-0.00275
        cg2=-0.00625
        phi1=-0.7861
        phi2=-0.0699
        phi3=-0.008444
        phi4=5.41000
        phi5=0.2899
        phi6=0.006718
        phi7=459.0
        phi8=0.1138
    ######################################################
    #Compute Preliminary Values
    ######################################################
    #No Z1.0 database exists, so we use the emperical function
    #provided in CY08 (eq 1) that relates to VS30
    Z10 = numpy.exp(28.5-3.82/8*numpy.log(numpy.power(VS30,8)+math.pow(378.7,8)))
    V1 = 1800.0
    b = numpy.where(VS30 <= 1130.0,phi2*(numpy.exp(phi3*(VS30-360.0)) - math.exp(phi3*(1130.0-360.0))), phi2*(numpy.exp(phi3*(1130.0-360.0)) - math.exp(phi3*(1130.0-360.0))))
    ##################################################
    #Style of faulting term
    ##################################################
    #Reverse faulting flag is 30 to 150
    #Normal faulting flag is -120 to -60 or 240 to 300
    NS=0
    RS=0
    if rake < 0:
        rake = 360+rake
    if  rake >= 30 and rake <= 150:
        NS=0
        RS=1
    if rake >= 240 and rake <= 300:
        NS=1
        RS=0
    ######################################################
    #Reference values
    ######################################################
    t1 = c1
    ######################################################
    t2 = c1a*RS + c1b*NS
    ######################################################
    #if an aftershock, this is irrelevant. Not coding it.
    t3 = c7*(Ztor - 4)
    ######################################################
    t4 = c2*(M-6)
    ######################################################
    t5 = (c2-c3)/cn*math.log(1+math.exp(cn*(cm-M)))
    ######################################################
    if (M-chm) > 0:
        t6 = c4*numpy.log(Rrup+c5*math.cosh(c6*(M-chm)))
    else:
        t6 = c4*numpy.log(Rrup+c5*math.cosh(0.0))
    ######################################################
    t7 = (c4a-c4)*numpy.log(numpy.sqrt(Rrup**2+crb**2))
    ######################################################
    if (M-cg3) > 0:
        t8 = (cg1 + cg2/math.cosh(M-cg3))*Rrup
    else:
        t8 = (cg1 + cg2/math.cosh(0.0))*Rrup
    ######################################################
    t9 = numpy.where(Hw == 1, c9*numpy.tanh(Rx*math.cos(dip)*math.cos(dip)/c9a)*(1-numpy.sqrt(numpy.power(Rjb,2)+math.pow(Ztor,2))/(Rrup+0.001)), 0.0)
    ######################################################
    #Reference Value is the sum of first 9 components
    yref = numpy.exp(t1+t2+t3+t4+t5+t6+t7+t8+t9)
    ######################################################
    #Rest of PGA terms
    ######################################################
    t11 = numpy.log(yref)
    ######################################################
    #t21 + t31 is the site response, eq 10
    t21 = numpy.where(numpy.log(VS30/1130) < 0, phi1*numpy.log(VS30/1130),0.0)
    ###############################################
    t31 = b*numpy.log((yref+phi4)/phi4)
    ######################################################
    #t41 is the deep sediment term, eq 11
    t41 = numpy.where(Z10-phi7 >= 0,phi5*(1-1/numpy.cosh(phi6*(Z10-phi7))),phi5*(1-1/math.cosh(phi6*0.0)))
    ######################################################
    #t51 is the shallow sediment term, eq 12
    t51 = numpy.where(Z10-15.0 >= 0, phi8/numpy.cosh(0.15*(Z10-15.0)), phi8/math.cosh(0.15*(0.0)))
    ######################################################         
    y = numpy.exp(t11+t21+t31+t41+t51)

    return(y)


def gmpe_ba08(M,Rjb,VS30,rake,mode):
    a1=0.03
    pga_low=0.06
    a2=0.09
    v1=180.0
    v2=300.0
    vref=760.0
    Mref=4.5
    Rref=1.0
    if mode == 0:
        #PGA coefficients
        blin=-0.360
        b1=-0.640
        b2=-0.14
        c1=-0.66050
        c2=0.11970
        c3=-0.01151
        h=1.35
        e1=-0.53804
        e2=-0.50350
        e3=-0.75472
        e4=-0.50970
        e5=0.28805
        e6=-0.10164
        e7=0.0
        Mh=6.75
    if mode == 1:
        #PGV coefficients
        blin=-0.600
        b1=-0.500
        b2=-0.06
        c1=-0.87370
        c2=0.10060
        c3=-0.00334
        h=2.54
        e1=5.00121
        e2=5.04727
        e3=4.63188
        e4=5.08210
        e5=0.18322
        e6=-0.12736
        e7=0.0
        Mh=8.5
    ##################################################
    #Eq 1 gives general form. There is a magnitude scaling,
    #a distance scaling, and a site amplification term
    #Prelim equations
    ##################################################
    R = numpy.sqrt(numpy.power(Rjb,2)+h*h)
    ##################################################
    #Style of faulting term
    ##################################################
    if rake < 0:
        rake = 360+rake
    if abs(rake >= 0 and rake <= 30):
        SS=1
        U=0
        NS=0
        RS=0
    if abs(rake >= 150 and rake <= 210):
        SS=1
        U=0
        NS=0
        RS=0
    if abs(rake >= 330 and rake <= 360):
        SS=1
        U=0
        NS=0
        RS=0
    if rake > 30 and rake < 150:
        SS=0
        U=0
        NS=0
        RS=1
    if rake > 210 and rake < 330:
        SS=0
        U=0
        NS=1
        RS=0
    ##################################################
    #Magnitude Term (Eq. 5)
    ##################################################
    if M <= Mh:
        fmag = e1*U + e2*SS + e3*NS + e4*RS + e5*(M-Mh) + e6*(math.pow(M-Mh,2))
    if M > Mh:
        fmag = e1*U + e2*SS + e3*NS + e4*RS + e7*(M-Mh)
    ##################################################
    #Distance term (Eq. 3)
    ##################################################
    fdis = (c1 + c2*(M-Mref))*numpy.log(R/Rref) + c3*(R-Rref)
    ##################################################
    #Linear Site response term
    ##################################################
    flin = blin*numpy.log(VS30/vref)
    ##################################################
    #pga4nl term - need to recompute magnitude and distance term
    ##################################################
    pga4nl = pga4nl_calc(M,Rjb,U,SS,NS,RS)
    ##################################################
    #Nonlinear slope term, bnl
    ##################################################
    bnl2 = (b1-b2)*numpy.log(VS30/v2)/math.log(v1/v2) + b2
    bnl3 = b2*numpy.log(VS30/vref)/numpy.log(v2/vref)
    
    bnl = numpy.where(v1 >= VS30,b1,0)
    bnl = numpy.where(VS30 > v1,bnl2,bnl)
    bnl = numpy.where(VS30 > v2,bnl3,bnl)
    bnl = numpy.where(VS30 > vref,0,bnl)
    #Other necessary nonlinear terms
    dx = math.log(a2/a1)
    dy = bnl*math.log(a2/pga_low)
    c = (3*dy - bnl*dx)/math.pow(dx,2)
    d = -(2*dy-bnl*dx)/math.pow(dx,3)
    ##################################################
    #Nonlinear site response term
    ##################################################
    fnl = 0.0*VS30
    fnl1 = bnl*math.log(pga_low/0.1)
    fnl2 = bnl*math.log(pga_low/0.1) + c*numpy.power(numpy.log(pga4nl/a1),2) + d*numpy.power(numpy.log(pga4nl/a1),3)
    fnl3 = bnl*numpy.log(pga4nl/0.1)
    
    fnl = numpy.where(pga4nl <= a1,fnl1,0.0)
    fnl = numpy.where(pga4nl > a1,fnl2,fnl)
    fnl = numpy.where(pga4nl > a2,fnl3,fnl)
    ##################################################
    #All terms
    ##################################################
    y = numpy.exp(fmag+fdis+flin+fnl)
    
    return(y)


def gmpe_cb14(M,Rjb,Rrup,Rx,Ztor,Zhyp,VS30,Hw,W,dip,rake,mode):
    #PGA coefficients
    if mode == 0:
        c0 = -4.416
        c1 = 0.984
        c2 = 0.537
        c3 = -1.499
        c4 = -0.496
        c5 = -2.773
        c6 = 0.248
        c7 = 6.768
        c8 = 0.0
        c9 = -0.212
        c10 = 0.720
        c11 = 1.090
        c12 = 2.186
        c13 = 1.420
        c14 = -0.0064
        c15 = -0.202
        c16 = 0.393
        c17 = 0.0977
        c18 = 0.0333
        c19 = 0.00757
        c20 = -0.0055
        dc20ji = -0.0035
        dc20ch = 0.0036
        dc20 = 0.0
        k1 = 865.0
        k2 = -1.186
        k3 = 1.839
        a2 = 0.167
        h1 = 0.241
        h2 = 1.474
        h3 = -0.715
        h5 = -0.337
        h6 = -0.270

        c = 1.88
        n = 1.18
        h4 = 1.0

        phi1 = -0.4417
        phi2 = -0.1417
        phi3 = -0.00701
        phi4 = 0.102151
        phi5 = 0.2289
        phi6 = 0.014996
        phi7 = 580.0
        phi8 = 0.07
        tau1 = 0.3437
        tau2 = 0.2637
        sigma1 = 0.4458
        sigma2 = 0.3459
        sigma3 = 0.8

        Sj = 0.0
        A1100 = 0.1
    #PGV coefficients
    if mode == 1:
        c0 = -2.895
        c1 = 1.510
        c2 = 0.270
        c3 = -1.299
        c4 = -0.453
        c5 = -2.466
        c6 = 0.204
        c7 = 5.837
        c8 = 0.0
        c9 = -0.168
        c10 = 0.305
        c11 = 1.713
        c12 = 2.602
        c13 = 2.457
        c14 = 0.1060
        c15 = 0.332
        c16 = 0.585
        c17 = 0.0517
        c18 = 0.0327
        c19 = 0.00613
        c20 = -0.0017
        dc20ji = -0.0006
        dc20ch = 0.0017
        dc20 = 0.0
        k1 = 400.0
        k2 = -1.955
        k3 = 1.929
        a2 = 0.596
        h1 = 0.117
        h2 = 1.616
        h3 = -0.733
        h5 = -0.128
        h6 = -0.756

        c = 1.88
        n = 1.18
        h4 = 1.0

        phi1 = -0.4417
        phi2 = -0.1417
        phi3 = -0.00701
        phi4 = 0.102151
        phi5 = 0.2289
        phi6 = 0.014996
        phi7 = 580.0
        phi8 = 0.07
        tau1 = 0.3437
        tau2 = 0.2637
        sigma1 = 0.4458
        sigma2 = 0.3459
        sigma3 = 0.8

        Sj = 0.0
        A1100 = 0.1
    ##################################################
    #Z2.5 Relationship (eq 33)
    #There is no west coast database of Z2.5, so we use
    #the relationship from CB14 to derive it based on Vs30
    ##################################################
    Z25 = numpy.exp(7.089-1.144*numpy.log(VS30))
    ##################################################
    #Magnitude Term (Eq. 2)
    ##################################################
    if M <= 4.5:
        fmag = c0 + c1*M
    if M > 4.5 and M <= 5.5:
        fmag = c0 + c1*M + c2*(M-4.5)
    if M > 5.5 and M <= 6.5:
        fmag = c0 + c1*M + c2*(M-4.5) + c3*(M-5.5)
    if M > 6.5:
        fmag = c0 + c1*M + c2*(M-4.5) + c3*(M-5.5) + c4*(M-6.5)
    ##################################################
    #Geometric attenuation term (Eq. 3)
    ##################################################
    fdis = (c5+c6*M)*numpy.log(numpy.sqrt(numpy.power(Rrup,2)+c7*c7))
    ##################################################
    #Style of faulting term (Eqs 4-6)
    ##################################################
    if rake > 30 and rake < 150:
        FRV = 1
    else:
        FRV = 0
        
    if rake > -150 and rake < -30:
        FNM = 1
    else:
        FNM = 0

    if M <= 4.5:
        ffltm = 0
    if M > 4.5 and M <= 5.5:
        ffltm = M-4.5
    if M > 5.5:
        ffltm = 1

    ffltf = c8*FRV + c9*FNM

    fflt = ffltf*ffltm
    ##################################################
    #Hanging wall term (Eqs. 7-16)
    ##################################################
    fhw = numpy.zeros([len(Hw),1])
    for i in range (0, len(Hw)):
        if Hw[i] == 1:
            #####################
            #Prelim equations
            #####################
            fhngdip = (90-dip)/45
            R1 = W*math.cos(math.radians(dip))
            R2 = 62*M - 350
            f1Rx = h1 + h2*Rx[i]/R1 + h3*math.pow(Rx[i]/R1,2)
            f2Rx = h4 + h5*((Rx[i]-R1)/(R2-R1)) + h6*math.pow((Rx[i]-R1)/(R2-R1),2)
            #####################
            #Equation 8
            #####################
            if (Rx[i] < R1):
                fhngrx = f1Rx
            else:
                if f2Rx > 0:
                    fhngrx = f2Rx
                else:
                    fhngrx = 0
            #####################
            #Equation 13
            #####################
            if Rrup[i] == 0:
                fhngrrup = 1
            else:
                fhngrrup = (Rrup[i]-Rjb[i])/Rrup[i]
            #####################
            #Equation 14
            #####################
            if M <= 5.5:
                fhngm = 0
            if M > 5.5 and M <= 6.5:
                fhngm = (M-5.5)*(1+a2*(M-6.5))
            if M > 6.5:
                fhngm = 1+a2*(M-6.5)
            #####################
            #Equation 15
            #####################
            if Ztor <= 16.66:
                fhngz = 1-0.06*Ztor
            else:
                fhngz = 0.0
                
            #####################
            #Equation 16
            #####################
            fhngdip = (90-dip)/45
            
            #####################
            #Equation 7
            #####################
            fhw[i] = c10*fhngrx*fhngrrup*fhngm*fhngz*fhngdip
    ##################################################
    #Shallow Site Response Term (Eqs. 17-19)
    ##################################################
    fsite = numpy.zeros([len(VS30),1])
    for i in range (0, len(VS30)):
        if (VS30[i] <= k1):
            fsiteG = c11*numpy.log(VS30[i]/k1) + k2* ( numpy.log(A1100+c*numpy.power(VS30[i]/k1,n))-numpy.log(A1100+c))
        if (VS30[i] > k1):
            fsiteG = (c11+k2*n)*numpy.log(VS30[i]/k1)
        #For J site classification
        if (VS30[i] <= 200.0):
            fsiteJ = (c12+k2*n)*(numpy.log(VS30[i]/k1)-numpy.log(200.0/k1))
        else:
            fsiteJ = (c13 + k2*n)*numpy.log(VS30[i]/k1)

        fsite[i] = fsiteG+Sj*fsiteJ
    ##################################################
    #Basin response term (Eq. 20)
    ##################################################
    fsed = numpy.zeros([len(Z25),1])
    for i in range (0, len(Z25)):
        if Z25[i] <= 1:
            fsed[i] = (c14+c15*Sj)*(Z25[i]-1)
        if Z25[i] > 1 and Z25[i] <= 3:
            fsed[i] = 0
        if Z25[i] > 3:
            fsed[i] = c16*k3*math.exp(-0.75)*(1-math.exp(-0.25*(Z25[i]-3)))
    ##################################################
    #Hypocentral Depth term (Eqs. 21-23)
    ##################################################           
    if Zhyp <= 7:
        fhyph = 0
    if (Zhyp > 7 and Zhyp <= 20):
        fhyph = Zhyp-7
    if (Zhyp > 20):
        fhyph = 13.0
        
    if M <= 5.5:
        fhypm = c17
    if M > 5.5 and M <= 6.5:
        fhypm = c17 + (c18-c17)*(M-5.5)
    if M > 6.5:
        fhypm = c18

    fhyp = fhyph*fhypm
    ##################################################
    #Fault dip term (Eq 24)
    ################################################## 
    if M <= 4.5:
        fdip = c19*dip
    if M > 4.5 and M <= 5.5:
        fdip = c19*(5.5-M)*dip
    if M > 5.5:
        fdip = 0
    ##################################################
    #Anelastic attenuation term (Eq 25)
    ################################################## 
    fatn = numpy.zeros([len(Rrup),1])
    for i in range (0, len(Rrup)):
        if Rrup[i] > 80:
            fatn[i] = (c20+dc20)*(Rrup[i] - 80.0)
    ######################################################

    y = numpy.exp(fmag+fdis+fflt+fhw+fsite+fsed+fhyp+fdip+fatn)
    
    return(y)


def gmpe_bradley13(M,Rrup,Rjb,Rx,Rtvz,Ztor,VS30,Z10,Hw,dip):
    c1 = -1.1985
    c1a = 0.1
    c1b = -0.455
    c2 = 1.06
    c3 = 1.500
    cm = 5.85
    cn = 2.996
    c4 = -2.1
    c4a = -0.5
    crb = 50.0
    c5 = 6.16
    c6 = 0.4893
    chm = 3.0
    c7 = 0.0512
    c8 = 10.0
    c9 = 0.79
    c9a = 1.5005
    cg1 = -0.0096
    cg2 = -0.0048
    cg3 = 4.0
    ctvz = 2.0

    phi1 = -0.4417
    phi2 = -0.1417
    phi3 = -0.00701
    phi4 = 0.102151
    phi5 = 0.2289
    phi6 = 0.014996
    phi7 = 580.0
    phi8 = 0.07
    tau1 = 0.3437
    tau2 = 0.2637
    sigma1 = 0.4458
    sigma2 = 0.3459
    sigma3 = 0.8
    ######################################################
    #Compute Preliminary Values
    ######################################################
    V1 = 1800.0
    if VS30 <= 1130.0:
        b = phi2*(math.exp(phi3*(VS30-360.0)) - math.exp(phi3*(1130.0-360.0)))
    else:
        b = 0.0
    ######################################################
    #New Zealand reference PGA (each t value is different term)
    ######################################################
    t1 = c1
    ######################################################
    t2 = c1a
    ######################################################
    if float(Ztor) <= c8:
        t3 = c7*(Ztor - 4)
    else:
        t3 = c7*(c8-4)
    ######################################################
    t4 = c2*(M-6)
    ######################################################
    t5 = (c2-c3)/cn*math.log(1+math.exp(cn*(cm-M)))
    ######################################################
    if (M-chm) > 0:
        t6 = c4*math.log(Rrup+c5*math.cosh(c6*(M-chm)))
    else:
        t6 = c4*math.log(Rrup+c5*math.cosh(0.0))
    ######################################################
    t7 = (c4a-c4)*math.log(math.sqrt(Rrup**2+crb**2))
    ######################################################
    if (M-cg3) > 0:
        t8 = (cg1 + cg2/math.cosh(M-cg3))*(1+ctvz*Rtvz/Rrup)*Rrup
    else:
        t8 = (cg1 + cg2/math.cosh(0.0))*(1+ctvz*Rtvz/Rrup)*Rrup
    ######################################################
    if Hw == 1:
        t9 = c9*math.tanh(Rx*math.cos(dip)*math.cos(dip)/c9a)*(1-math.sqrt(math.pow(Rjb,2)+math.pow(Ztor,2))/(Rrup+0.001))
    else:
        t9 = 0
    ######################################################

    yref = math.exp(t1+t2+t3+t4+t5+t6+t7+t8+t9)

    ######################################################
    #Rest of PGA terms
    ######################################################
    t11 = math.log(yref)
    ######################################################
    if VS30 < V1:
        t21 = phi1*math.log(VS30/1130)
    else:
        t21 = phi1*math.log(V1/1130)
    ######################################################
    t31 = b*math.log((yref+phi4)/phi4)
    ######################################################
    if Z10-phi7 >= 0:
        t41 = phi5*(1-1/math.cosh(phi6*(Z10-phi7)))
    else:
        t41 = phi5*(1-1/math.cosh(phi6*0.0))
    ######################################################
    if Z10-15.0 >= 0:
        t51 = phi8/math.cosh(0.15*(Z10-15.0))
    else:
        t51 = phi8/math.cosh(0.15*(0.0))
    ######################################################

    y = math.exp(t11+t21+t31+t41+t51)
        
    return(y)

def pga4nl_calc(M,Rjb,U,SS,NS,RS):
    a1=0.03
    pga_low=0.06
    a2=0.09
    v1=180.0
    v2=300.0
    vref=760.0
    Mref=4.5
    Rref=5.0
    
    blin=-0.360
    b1=-0.640
    b2=-0.14
    c1=-0.66050
    c2=0.11970
    c3=-0.01151
    h=1.35
    e1=-0.53804
    e2=-0.50350
    e3=-0.75472
    e4=-0.50970
    e5=0.28805
    e6=-0.10164
    e7=0.0
    Mh=6.75

    R= numpy.sqrt(numpy.power(Rjb,2)+h*h)
    
    if M <= Mh:
        fmag = e1*U + e2*SS + e3*NS + e4*RS + e5*(M-Mh) + e6*(math.pow(M-Mh,2))
    if M > Mh:
        fmag = e1*U + e2*SS + e3*NS + e4*RS + e7*(M-Mh)

    fdis = c1 + c2*(M-Mref)*numpy.log(R/Rref) + c3*(R-Rref)

    ftot = numpy.exp(fmag+fdis)
    
    return(ftot)

