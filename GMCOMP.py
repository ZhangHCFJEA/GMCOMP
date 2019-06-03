#!/usr/bin/python
###########################################################
#GMCOMP - A ground motion comparison package
#to compare results from different source models
#derived for ShakeAlert. Uses point source coreinfo xml
#messages, finite fault xml messages, and raw 
#accelerometer data in miniseed format. Can be used for
#generic ground motion and intensity prediction as long as
#event information is in standard coreinfo format or
#finitefault message format.
########
#A secondary mode compares the surface deformation from
#finite fault messages and srcmod formatted slip models
########
#Brendan Crowell, University of Washington
#Last Updated June 3, 2019
###########################################################
import math
import numpy
import GMPE
import csv
from GMCOMP_paraminit import Properties
import MSG_xmlreader
from MSG_srcmod_reader import srcmod_fspreader
import SM_peakmotions
import obspy
from GMPE_faultdist import sourcedistance, sourcedistance_point, sourcedistance_line
from GMICE import MMI_Y
import GMCOMP_coord_tools
import DEFCOMP
import GMCOMP_traveltime
import GMCOMP_gmtplots
import os
###########################################################
#Read Properties file, create dict 'props'
###########################################################
props = Properties('gmcomp.props')
chanfile = props.getaccchanfile() #where ground motions are computed
accdir = props.getaccdatadir()
outdir = props.getoutdir()
if not os.path.exists(outdir): #if output folder doesn't exist, make it
    os.makedirs(outdir)
xmldir = props.getxmldir() #directory of xml files

#coreinfo files and algorithms
numcorefiles = props.getnumcoremessages()
corefiles = props.getcoremessages()
corealgs = props.getcorealgs()

#eqinfo files and algorithms
numeqinfofiles = props.getnumeqinfo2gmmessages()
eqinfofiles = props.geteqinfo2gmmessages()
eqinfoalgs = props.geteqinfo2gmalgs()

#srcmod files and algorithms and srcmod directory
numsrcmodfiles = props.getnumsrcmodfiles()
srcmodfiles = props.getsrcmodfiles()
srcmodalgs = props.getsrcmodalgs()
srcmoddir = props.getsrcmoddir()

#shakemap files and algorithms
numshakemapfiles = props.getnumshakemapfiles()
shakemapfiles = props.getshakemapfiles()
shakemapalgs = props.getshakemapalgs()

#finite fault files and algorithms
numfffiles = props.getnumffmessages()
fffiles = props.getffmessages()
ffalgs = props.getffalgs()

#line source files and algorithms
numlsfiles = props.getnumlsmessages()
lsfiles = props.getlsmessages()
lsalgs = props.getlsalgs()

#vs30 file
vs30file  = props.getvs30file()

#Origin Time Information, convert to proper obspy format
yr = int(props.getotyr())
mo = int(props.getotmo())
dy = int(props.getotdy())
hr = int(props.getothr())
mn = int(props.getotmn())
sc = float(props.getotsc())
origintime = obspy.core.utcdatetime.UTCDateTime(yr,mo,dy,hr,mn,sc)
#earthquake name to use on files
earthquake=str(props.geteq())
#get mmithreshold and mmiwarnthreshold. The mmiwarnthreshold is
#to define the warning level (i.e. MMI 4 in most cases) to
#determine false/true positives/negatives
#mmithreshold is the instrumental mmi to use for ground motion
#comparisons. Rather than include all available stations, we can
#cull the list to only include stations that exceed this threshold
#but are under the warning threshold. Usually set to MMI=1
mmithreshold=float(props.getmmithreshold())
mmiwarnthreshold=float(props.getmmiwarnthreshold())
#this keeps a list of mmi evolution files and algorithm names
evolfiles=list() 
algorithm=list()
#ttot is the total number of seconds after origin time to
#consider
ttot=int(props.getttot())
###########################################################
###########################################################
#Compute true peak ground motions and evolution
###########################################################
#Open evolution and max recorded motion files and write
#headers for them.
fnamesmevol = outdir + earthquake + '_sm_mmi_evol.txt'
evolfiles.append(fnamesmevol)
algorithm.append('sm')
fnamesmmax = outdir + earthquake + '_sm_mmi_max.txt'

runminiseed = props.getrunmseed()
if (runminiseed == 'yes'):
    fsm = open(fnamesmevol,'w')
    fsmmax = open(fnamesmmax,'w')
    fsm.write('##########################################################################'+'\n')
    fsm.write('#'+'Strong Motion Evolution Results - Raw Recordings'+'\n')
    fsm.write('##########################################################################'+'\n')
    fsm.write('#'+'Earthquake'+' '+earthquake+'\n')
    fsm.write('#'+'Chanfile'+' '+chanfile+'\n')
    fsm.write('#'+'Origin Time'+' '+str(origintime)+'\n')
    fsm.write('#'+'Only MMI > '+ props.getmmithreshold() +' reported'+'\n')
    fsm.write('#'+'MMI and PGA component used ='+' '+props.getmmicomp()+'\n')
    fsm.write('#'+'Values not reported unless MMI is greater than previous reporting'+'\n')
    fsm.write('#'+'site'+' '+'lon'+' '+'lat'+' '+'time after OT(s)'+' '+'pga(cm/s2)'+' '+'mmi'+'\n')
    fsm.write('##########################################################################'+'\n')


    fsmmax.write('##########################################################################'+'\n')
    fsmmax.write('#'+'Strong Motion Max Results - Raw Recordings'+'\n')
    fsmmax.write('##########################################################################'+'\n')
    fsmmax.write('#'+'Earthquake'+' '+earthquake+'\n')
    fsmmax.write('#'+'Chanfile'+' '+chanfile+'\n')
    fsmmax.write('#'+'Origin Time'+' '+str(origintime)+'\n')
    fsmmax.write('#'+'Only MMI > '+ props.getmmithreshold() +' reported'+'\n')
    fsmmax.write('#'+'MMI and PGA component used ='+' '+props.getmmicomp()+'\n')
    fsmmax.write('#'+'site'+' '+'lon'+' '+'lat'+' '+' '+'pga(cm/s2)'+' '+'mmi'+'\n')
    fsmmax.write('##########################################################################'+'\n')
    #Open up chanfile, compute PGA values once 3 channels are available for a station. The channel
    #codes must be the same for the first two letters, i.e. HL or HN or EN
    k = 0
    with open(chanfile, 'rt') as chan:
        reader = csv.reader(chan, delimiter="\t")
        for row in reader:
            if "#" not in row[0]:
                smnetcode = row[0]
                smsite = row[1]
                smloccode = row[2]
                smchannel = row[3]
                smlat = row[4]
                smlon = row[5]
                smelev = row[6]
                smsamplerate = row[7]
                smgain = row[8]
                smunits = row[9]
                print('Loading strong motion station ', str(smsite))
                if (smunits == 'counts/(cm/sec2)'):
                    if (smchannel == 'HNE' or smchannel == 'HLE' or smchannel == 'ENE'):
                        fname_east = accdir + smsite.strip() + '.' + smchannel.strip() + '.mseed'
                        site_east = smsite
                        gain_east = smgain
                        chan_east = smchannel
                        k=k+1
                    if (smchannel == 'HNN' or smchannel == 'HLN' or smchannel == 'ENN'):
                        fname_north = accdir + smsite.strip() + '.' + smchannel.strip() + '.mseed'
                        site_north = smsite
                        gain_north = smgain
                        chan_north = smchannel
                        k=k+1
                    if (smchannel == 'HNZ' or smchannel == 'HLZ' or smchannel == 'ENZ'):
                        fname_vert = accdir + smsite.strip() + '.' + smchannel.strip() + '.mseed'
                        site_vert = smsite
                        gain_vert = smgain
                        chan_vert = smchannel
                        k=k+1
                    if (k == 3 and site_east == site_north and site_east == site_vert):
                        if (chan_east[0:2] == chan_north[0:2] and chan_east[0:2] == chan_vert[0:2]):
                            pga_vals = SM_peakmotions.pga_measure(fname_east,fname_north,fname_vert,gain_east,gain_north,gain_vert,origintime,props)
                            t = pga_vals[0]
                            PGA = pga_vals[1]
                            MMI = pga_vals[2]
                            PGAmax = pga_vals[3]
                            MMImax = pga_vals[4]
                            
                            lon = "{0:.4f}".format(float(smlon))
                            lat = "{0:.4f}".format(float(smlat))
                            pgamax = "{0:.4f}".format(float(PGAmax))
                            mmimax = "{0:.2f}".format(float(MMImax))
                            if (MMImax > mmithreshold):
                                fsmmax.write(site_east+' '+lon+' '+lat+' '+pgamax+' '+mmimax+'\n')
                            
                            mmiind = 0
                            for i in range(1,len(t)):
                                if (MMI[i] > mmithreshold):
                                    if (MMI[i] > mmiind):
                                        tout = "{0:.2f}".format(float(t[i]))
                                        pgaout = "{0:.4f}".format(float(PGA[i]))
                                        mmiout = "{0:.2f}".format(float(MMI[i]))

                                        fsm.write(site_east+' '+lon+' '+lat+' '+tout+' '+pgaout+' '+mmiout+'\n')
                                        mmiind = MMI[i]
                        k=0
                        site_east=[]
                        site_north=[]
                        site_vert=[]
                    elif (k > 3):
                        k=0
    fsm.close()
    fsmmax.close()
###########################################################
#Station Locations to compute ground motions at
###########################################################
#get initial file size. Note, this is off the pga/mmi max file
#not the chan file since the chan file contains sites under
#the mmi threshold
k = 0
with open(fnamesmmax) as f:
    for line in f:
        if line.startswith("#"):
            pass
        else:
            k=k+1

#Create a list of site names and arrays of lon,lat,pga,mmi            
SITE = list()
LONS = numpy.nan*numpy.ones([k,1])
LATS = numpy.nan*numpy.ones([k,1])
PGAS = numpy.nan*numpy.ones([k,1])
MMIS = numpy.nan*numpy.ones([k,1])
k=0
with open(fnamesmmax) as f:
    for line in f:
        if line.startswith("#"):
            pass
        else:
            vals = line.split()
            SITE.append(vals[0])
            LONS[k,0] = float(vals[1])
            LATS[k,0] = float(vals[2])
            PGAS[k,0] = float(vals[3])
            MMIS[k,0] = float(vals[4])
            k=k+1
print ('Number of stations considered ', len(LONS))
###########################################################
#Compute theoretical travel time to each station from
#predefined lat/lon/dep
###########################################################
eqlon = float(props.geteqlon())
eqlat = float(props.geteqlat())
eqdep = float(props.geteqdep())

[tp,ts]=GMCOMP_traveltime.tt(eqlon,eqlat,eqdep,LONS,LATS)
###########################################################
#Read VS30 information at sites
###########################################################
print('Loading VS30 information from file ', str(vs30file))
with open(vs30file) as f:
    for i, l in enumerate(f):
        pass
    vs30_length = i+1

vslon = numpy.nan*numpy.ones([vs30_length,1])
vslat = numpy.nan*numpy.ones([vs30_length,1])
vs30 = numpy.nan*numpy.ones([vs30_length,1])

ind = 0
for line in open(vs30file, 'r'):
    data = line.split()
    vslon[ind,0] = data[1]
    vslat[ind,0] = data[0]
    vs30[ind,0] = data[2]
    ind=ind+1

[x_vs30,y_vs30,z_vs30]=GMCOMP_coord_tools.lla2ecef(vslat,vslon,0)

VS30_val = numpy.nan*numpy.ones([len(LONS),1])
for i in range(0, len(LONS)):
    [x,y,z]=GMCOMP_coord_tools.lla2ecef(float(LATS[i]),float(LONS[i]),0.0)

    dist = numpy.sqrt(numpy.power(x_vs30-x,2)+numpy.power(y_vs30-y,2)+numpy.power(z_vs30-z,2))
    ind = numpy.argmin(dist)
    VS30_val[i,0]=vs30[ind,0]
###########################################################
#Read Core info Messages and compute ground motions
###########################################################
for i in range(0, 1):
    cf = xmldir+str(corefiles[i])
    print('Loading coreinfo file ', str(cf))
    cfoutfile = outdir + str(earthquake) + '_' + str(corealgs[i]) +  '_' + str(props.getgmpe()) + '_' + str(props.getgmice()) +  '_' + str(props.getmmicomp()) + '_out.txt'
    evolfiles.append(cfoutfile)
    algorithm.append(str(corealgs[i]))
    fcf = open(cfoutfile,'w')
    fcf.write('##########################################################################'+'\n')
    fcf.write('#'+'Strong Motion Evolution Results - Coreinfo File'+'\n')
    fcf.write('##########################################################################'+'\n')
    fcf.write('#'+'Earthquake'+' '+earthquake+'\n')
    fcf.write('#'+'Corefile'+' '+str(corefiles[i])+'\n')
    fcf.write('#'+'Algorithm'+' '+str(corealgs[i]) +'\n')
    fcf.write('#'+'Origin Time'+' '+str(origintime)+'\n')
    fcf.write('#'+'Only MMI > '+ props.getmmithreshold() +' used for comparisons'+'\n')
    fcf.write('#'+'MMI and PGA component used ='+' '+props.getmmicomp()+'\n')
    fcf.write('#'+'site'+','+'lon'+','+'lat'+','+'OT+(s)'+','+'pga(cm/s2)'+','+'mmi'+','+'mmibias'+','+'MV'+','+'Mw'+'\n')
    fcf.write('##########################################################################'+'\n')
    k = 0
    for line in open(cf):
        [mess_overview, mess_details] = MSG_xmlreader.cireader(xmldir+str(line.strip()))
        ot_off = obspy.core.utcdatetime.UTCDateTime(mess_details[9])-origintime #CI origin time offset
        rt_off = obspy.core.utcdatetime.UTCDateTime(mess_overview[5])-origintime #CI runtime offset
        CImag = float(mess_details[1])
        CIlon = float(mess_details[5])
        CIlat = float(mess_details[3])
        CIdep = float(mess_details[7])
        CIos = mess_overview[4]
        CIvers = mess_overview[6]
        [rrup, rjb, rx, rp, ztor, hw, W] = sourcedistance_point(CIlat, CIlon, CIdep, LATS, LONS, len(LATS), props)
        Y = GMPE.Y(CImag, rrup, rjb, rx, ztor, CIdep, VS30_val, hw, W, 90, 0, 0, props)
        MMI = MMI_Y(numpy.reshape(Y,(len(Y),1)), 0, props)

        for j in range(0, len(LATS)):
            lon = "{0:.4f}".format(float(LONS[j]))
            lat = "{0:.4f}".format(float(LATS[j]))
            pgaout = "{0:.4f}".format(float(Y[j]))
            mmiout = "{0:.2f}".format(float(MMI[j]))
            mmibias = "{0:.2f}".format(float(MMI[j]-MMIS[j]))
            mw = "{0:.2f}".format(float(CImag))
            tout = "{0:.2f}".format(float(rt_off-ot_off))
            fcf.write(str(SITE[j])+' '+lon+' '+lat+' '+tout+' '+pgaout+' '+mmiout+' '+mmibias+' '+CIvers+' '+mw+'\n')
        k = k+1
    fcf.close()
            
###########################################################
#Read FF Messages and compute ground motion/deformation field
###########################################################
for i in range(0, int(numfffiles)):
    fffile = xmldir+str(fffiles[i])
    print('Loading finite fault file ', str(fffile))
    ffoutfile = outdir + str(earthquake) + '_' + str(ffalgs[i]) +  '_' + str(props.getgmpe()) + '_' + str(props.getgmice()) +  '_' + str(props.getmmicomp()) + '_out.txt'
    evolfiles.append(ffoutfile)
    algorithm.append(ffalgs[i])
    
    fff = open(ffoutfile,'w')
    fff.write('##########################################################################'+'\n')
    fff.write('#'+'Strong Motion Evolution Results - Finite Fault File'+'\n')
    fff.write('##########################################################################'+'\n')
    fff.write('#'+'Earthquake'+' '+earthquake+'\n')
    fff.write('#'+'Corefile'+' '+str(fffiles[i])+'\n')
    fff.write('#'+'Algorithm'+' '+str(ffalgs[i]) +'\n')
    fff.write('#'+'Origin Time'+' '+str(origintime)+'\n')
    fff.write('#'+'Only MMI > '+ props.getmmithreshold() +' reported'+'\n')
    fff.write('#'+'MMI and PGA component used ='+' '+props.getmmicomp()+'\n')
    fff.write('#'+'site'+','+'lon'+','+'lat'+','+'OT+(s)'+','+'pga(cm/s2)'+','+'mmi'+','+'mmibias'+','+'MV'+','+'Mw'+'\n')
    fff.write('##########################################################################'+'\n')
    
    if (props.getdefmodon() == 'yes'):
        defoutfile = outdir + str(earthquake) + '_' + str(ffalgs[i]) +  '_' + 'deformation_out.txt'
        deff = open(defoutfile,'w')
        deff.write('##########################################################################'+'\n')
        deff.write('#'+'Forward Model Deformation - Finite Fault File'+'\n')
        deff.write('##########################################################################'+'\n')
        deff.write('#'+'Earthquake'+' '+earthquake+'\n')
        deff.write('#'+'Algorithm'+' '+str(ffalgs[i]) +'\n')
        deff.write('#'+'Origin Time'+' '+str(origintime)+'\n')
        deff.write('#'+'Grid lon from'+' '+str(float(CIlon)-0.25)+' '+str(float(CIlon)+0.25)+'\n')
        deff.write('#'+'Grid lat from'+' '+str(float(CIlat)-0.25)+' '+str(float(CIlat)+0.25)+'\n')
        deff.write('#'+'lon'+','+'lat'+','+'OT+(s)'+','+'E(m)'+','+'N(m)'+','+'U(m)'+'\n')
        deff.write('##########################################################################'+'\n')

    k=0
    for line in open(fffile):
        k=k+1

    k2=0
    for line in open(fffile):
        k2=k2+1
        [mess_overview, mess_details, gunc_details, LONFF, LATFF, DEPFF, SSFF, DSFF] = MSG_xmlreader.ffreader(xmldir+str(line.strip()))
        ot_off = obspy.core.utcdatetime.UTCDateTime(mess_details[9])-origintime #FF origin time offset
        rt_off = obspy.core.utcdatetime.UTCDateTime(mess_overview[5])-origintime #FF runtime offset
        FFmag = mess_details[1]
        nvert = mess_details[11]
        zhyp = float(mess_details[7])
        FFvers = mess_overview[6]
        
        [rrup, rjb, rx, rp, ztor, hw, rake, strike, dip, L, W] = sourcedistance(LATFF, LONFF, DEPFF, SSFF, DSFF, LATS, LONS, len(LATS), nvert, props)

        rruporrp=props.getrruporrp()
        if (rruporrp == 'rrup'):
            Y = GMPE.Y(FFmag, rrup, rjb, rx, ztor, zhyp, VS30_val, hw, W, dip, numpy.mean(rake), 0, props)
        else:
            Y = GMPE.Y(FFmag, rp, rjb, rx, ztor, zhyp, VS30_val, hw, W, dip, numpy.mean(rake), 0, props)
        MMI = MMI_Y(numpy.reshape(Y,(len(Y),1)), 0, props)

        if (props.getdefmodon() == 'yes'):
            if (k == k2):
                [lonm, latm, depm, strm, dipm, lm, wm]=DEFCOMP.fp_organize(LONFF, LATFF, DEPFF)
                [lon_grid, lat_grid, en, nn, un] = DEFCOMP.def_grid(float(props.getgridloncen()), float(props.getgridlatcen()), float(props.getgridlon()), float(props.getgridlat()), int(props.getgridnum()), lonm, latm, depm, SSFF, DSFF, strm, dipm, wm, lm, props)
                for j in range(0, len(lon_grid)):
                    lon = "{0:.4f}".format(float(lon_grid[j]))
                    lat = "{0:.4f}".format(float(lat_grid[j]))
                    enout = "{0:.4f}".format(float(en[j]))
                    nnout = "{0:.4f}".format(float(nn[j]))
                    unout = "{0:.4f}".format(float(un[j]))
                    tout = "{0:.2f}".format(float(rt_off-ot_off))
                    deff.write(lon+' '+lat+' '+tout+' '+enout+' '+nnout+' '+unout+'\n')
                

        for j in range(0, len(LATS)):
            lon = "{0:.4f}".format(float(LONS[j]))
            lat = "{0:.4f}".format(float(LATS[j]))
            pgaout = "{0:.4f}".format(float(Y[j]))
            mmiout = "{0:.2f}".format(float(MMI[j]))
            mmibias = "{0:.2f}".format(float(MMI[j]-MMIS[j]))
            mw = "{0:.2f}".format(float(FFmag))
            tout = "{0:.2f}".format(float(rt_off-ot_off))
            fff.write(str(SITE[j])+' '+lon+' '+lat+' '+tout+' '+pgaout+' '+mmiout+' '+mmibias+' '+FFvers+' '+mw+'\n')
    fff.close()
    if (props.getdefmodon() == 'yes'):
        deff.close()

###########################################################
#Read FF Line Source Messages and compute ground motion field
###########################################################
for i in range(0, int(numlsfiles)):
    lsfile = xmldir+str(lsfiles[i])
    print('Loading finite fault line source file ', str(lsfile))
    lsoutfile = outdir + str(earthquake) + '_' + str(lsalgs[i]) +  '_' + str(props.getgmpe()) + '_' + str(props.getgmice()) +  '_' + str(props.getmmicomp()) + '_out.txt'
    evolfiles.append(lsoutfile)
    algorithm.append(lsalgs[i])
    
    lsf = open(lsoutfile,'w')
    lsf.write('##########################################################################'+'\n')
    lsf.write('#'+'Strong Motion Evolution Results - Finite Fault Line Source File'+'\n')
    lsf.write('##########################################################################'+'\n')
    lsf.write('#'+'Earthquake'+' '+earthquake+'\n')
    lsf.write('#'+'Corefile'+' '+str(lsfiles[i])+'\n')
    lsf.write('#'+'Algorithm'+' '+str(lsalgs[i]) +'\n')
    lsf.write('#'+'Origin Time'+' '+str(origintime)+'\n')
    lsf.write('#'+'Only MMI > '+ props.getmmithreshold() +' reported'+'\n')
    lsf.write('#'+'MMI and PGA component used ='+' '+props.getmmicomp()+'\n')
    lsf.write('#'+'site'+','+'lon'+','+'lat'+','+'OT+(s)'+','+'pga(cm/s2)'+','+'mmi'+','+'mmibias'+','+'MV'+','+'Mw'+'\n')
    lsf.write('##########################################################################'+'\n')

    for line in open(lsfile):
        [mess_overview, mess_details, LONFF, LATFF, DEPFF] = MSG_xmlreader.ffreader_line(xmldir+str(line.strip()))
        ot_off = obspy.core.utcdatetime.UTCDateTime(mess_details[9])-origintime #FF origin time offset
        rt_off = obspy.core.utcdatetime.UTCDateTime(mess_overview[5])-origintime #FF runtime offset
        FFmag = mess_details[1]
       
        zhyp = float(mess_details[7])
        FFvers = mess_overview[6]
        
        [rrup, rjb, rx, rp, ztor, hw, W] = sourcedistance_line(LATFF, LONFF, DEPFF, LATS, LONS, len(LATS), props)

        Y = GMPE.Y(FFmag, rrup, rjb, rx, ztor, zhyp, VS30_val, hw, W, 90, 0, 0, props)

        MMI = MMI_Y(numpy.reshape(Y,(len(Y),1)), 0, props)

        for j in range(0, len(LATS)):
            lon = "{0:.4f}".format(float(LONS[j]))
            lat = "{0:.4f}".format(float(LATS[j]))
            pgaout = "{0:.4f}".format(float(Y[j]))
            mmiout = "{0:.2f}".format(float(MMI[j]))
            mmibias = "{0:.2f}".format(float(MMI[j]-MMIS[j]))
            mw = "{0:.2f}".format(float(FFmag))
            tout = "{0:.2f}".format(float(rt_off-ot_off))
            lsf.write(str(SITE[j])+' '+lon+' '+lat+' '+tout+' '+pgaout+' '+mmiout+' '+mmibias+' '+FFvers+' '+mw+'\n')
    lsf.close()
###########################################################
#Read eqinfo2gm Messages
###########################################################
for i in range(0, int(numeqinfofiles)):
    eqifile = xmldir+str(eqinfofiles[i])
    print('Loading eqinfo2gm file ', str(eqifile))
    eqioutfile = outdir + str(earthquake) + '_' + str(eqinfoalgs[i]) +  '_' + str(props.getgmpe()) + '_' + str(props.getgmice()) +  '_' + str(props.getmmicomp()) + '_out.txt'
    evolfiles.append(eqioutfile)
    algorithm.append(eqinfoalgs[i])
    feqi = open(eqioutfile,'w')

    feqi.write('##########################################################################'+'\n')
    feqi.write('#'+'Strong Motion Evolution Results - eqinfo2gm File'+'\n')
    feqi.write('##########################################################################'+'\n')
    feqi.write('#'+'Earthquake'+' '+earthquake+'\n')
    feqi.write('#'+'eqinfo2gm File'+' '+str(eqinfofiles[i])+'\n')
    feqi.write('#'+'Algorithm'+' '+str(eqinfoalgs[i]) +'\n')
    feqi.write('#'+'Origin Time'+' '+str(origintime)+'\n')
    feqi.write('#'+'Only MMI > '+ props.getmmithreshold() +' used for comparisons'+'\n')
    feqi.write('#'+'MMI and PGA component used ='+' '+props.getmmicomp()+'\n')
    feqi.write('#'+'site'+','+'lon'+','+'lat'+','+'OT+(s)'+','+'pga(cm/s2)'+','+'mmi'+','+'mmibias'+'\n')
    feqi.write('##########################################################################'+'\n')
    for line in open(eqifile):
        [mess_overview, mess_details, eqinfo_gm, eqinfotype] = MSG_xmlreader.eqinfo2gmreader(xmldir+str(line.strip()))
        
        ot_off = obspy.core.utcdatetime.UTCDateTime(mess_details[9])-origintime #origin time offset
        rt_off = obspy.core.utcdatetime.UTCDateTime(mess_overview[5])-origintime #runtime offset
        
        if (eqinfotype == 'map'):
            SMlon = numpy.asarray(eqinfo_gm[0])
            SMlat = numpy.asarray(eqinfo_gm[1])
            SMpga = numpy.asarray(eqinfo_gm[2])
            SMmmi = numpy.asarray(eqinfo_gm[3])

            for j in range(0, len(LATS)):
                [xsm,ysm]=GMCOMP_coord_tools.ll2utm(SMlon, SMlat, LONS[j], LATS[j])
                [xj,yj]=GMCOMP_coord_tools.ll2utm(LONS[j], LATS[j], LONS[j], LATS[j])
                smdist = numpy.sqrt(numpy.power(xsm-xj,2)+numpy.power(ysm-yj,2))
                
                a1 = numpy.where(numpy.amin(smdist) == smdist)[0]

                lon = "{0:.4f}".format(float(LONS[j]))
                lat = "{0:.4f}".format(float(LATS[j]))
                pgaout = "{0:.4f}".format(float(SMpga[a1]))
                mmiout = "{0:.2f}".format(float(SMmmi[a1]))
                mmibias = "{0:.2f}".format(float(SMmmi[a1]-MMIS[j]))
                tout = "{0:.2f}".format(float(rt_off-ot_off))
                feqi.write(str(SITE[j])+' '+lon+' '+lat+' '+tout+' '+pgaout+' '+mmiout+' '+mmibias+' '+'\n')
            k = k+1
    feqi.close()
###########################################################
#Read srcmod files
###########################################################
for i in range(0, int(numsrcmodfiles)):
    srcmodfile = srcmoddir+str(srcmodfiles[i])
    print('Loading srcmod file ', str(srcmodfile))
    srcmodoutfile = outdir + str(earthquake) + '_' + str(srcmodalgs[i]) +  '_' + str(props.getgmpe()) + '_' + str(props.getgmice()) +  '_' + str(props.getmmicomp()) + '_out.txt'
    evolfiles.append(srcmodoutfile)
    algorithm.append(srcmodalgs[i])
    if (props.getdefmodon() == 'yes'):
        defoutfile = outdir + str(earthquake) + '_' + str(srcmodalgs[i]) +  '_' + 'deformation_out.txt'
        deff = open(defoutfile,'w')
        deff.write('##########################################################################'+'\n')
        deff.write('#'+'Forward Model Deformation - Finite Fault File'+'\n')
        deff.write('##########################################################################'+'\n')
        deff.write('#'+'Earthquake'+' '+earthquake+'\n')
        deff.write('#'+'Algorithm'+' '+str(srcmodalgs[i]) +'\n')
        deff.write('#'+'Origin Time'+' '+str(origintime)+'\n')
        deff.write('#'+'Grid lon from'+' '+str(float(CIlon)-0.25)+' '+str(float(CIlon)+0.25)+'\n')
        deff.write('#'+'Grid lat from'+' '+str(float(CIlat)-0.25)+' '+str(float(CIlat)+0.25)+'\n')
        deff.write('#'+'lon'+','+'lat'+','+'OT+(s)'+','+'E(m)'+','+'N(m)'+','+'U(m)'+'\n')
        deff.write('##########################################################################'+'\n')
    
    smf = open(srcmodoutfile,'w')
    smf.write('##########################################################################'+'\n')
    smf.write('#'+'Strong Motion Evolution Results - Srcmod File'+'\n')
    smf.write('##########################################################################'+'\n')
    smf.write('#'+'Earthquake'+' '+earthquake+'\n')
    smf.write('#'+'SRCMOD file'+' '+str(srcmodfiles[i])+'\n')
    smf.write('#'+'Algorithm'+' '+str(srcmodalgs[i]) +'\n')
    smf.write('#'+'Origin Time'+' '+str(origintime)+'\n')
    smf.write('#'+'Only MMI > '+ props.getmmithreshold() +' reported'+'\n')
    smf.write('#'+'MMI and PGA component used ='+' '+props.getmmicomp()+'\n')
    smf.write('#'+'site'+','+'lon'+','+'lat'+','+'time after OT(s)'+','+'pga(cm/s2)'+','+'mmi'+','+'mmibias'+','+'Message Version'+'\n')
    smf.write('##########################################################################'+'\n')

    [rrup, rjb, rx, rp, ztor, hw, rake, strike, dip, L, W, fault_lat, fault_lon, fault_dep, fault_slip, fault_rake, fault_strike, fault_dip] = srcmod_fspreader(srcmodfile, LATS, LONS)
    rruporrp=props.getrruporrp()

    if (rruporrp == 'rrup'):
        Y = GMPE.Y(FFmag, rrup, rjb, rx, ztor, zhyp, VS30_val, hw, W, dip, numpy.mean(rake), 0, props)
    else:
        Y = GMPE.Y(FFmag, rp, rjb, rx, ztor, zhyp, VS30_val, hw, W, dip, numpy.mean(rake), 0, props)
    #Y = GMPE.Y(FFmag, rrup, rjb, rx, ztor, zhyp, VS30_val, hw, W, dip, rake, 0, props)
    MMI = MMI_Y(numpy.reshape(Y,(len(Y),1)), 0, props)
    SSFF = fault_slip*(-0.9903)
    DSFF = fault_slip*(-0.1392)

    FL = 0.75*numpy.ones([len(fault_slip),1])*1000
    FW = 0.625*numpy.ones([len(fault_slip),1])*1000
    if (props.getdefmodon() == 'yes'):
        [lon_grid, lat_grid, en, nn, un] = DEFCOMP.def_grid(float(props.getgridloncen()), float(props.getgridlatcen()), float(props.getgridlon()), float(props.getgridlat()), int(props.getgridnum()), numpy.asarray(fault_lon), numpy.asarray(fault_lat), 1000*numpy.asarray(fault_dep)+307, numpy.asarray(SSFF), numpy.asarray(DSFF), fault_strike, fault_dip, FW, FL, props)
        for j in range(0, len(lon_grid)):
            lon = "{0:.4f}".format(float(lon_grid[j]))
            lat = "{0:.4f}".format(float(lat_grid[j]))
            enout = "{0:.4f}".format(float(en[j]))
            nnout = "{0:.4f}".format(float(nn[j]))
            unout = "{0:.4f}".format(float(un[j]))
            tout = "{0:.2f}".format(float(rt_off-ot_off))
            deff.write(lon+' '+lat+' '+tout+' '+enout+' '+nnout+' '+unout+'\n')
 

    for j in range(0, len(LATS)):
        lon = "{0:.4f}".format(float(LONS[j]))
        lat = "{0:.4f}".format(float(LATS[j]))
        pgaout = "{0:.4f}".format(float(Y[j]))
        mmiout = "{0:.2f}".format(float(MMI[j]))
        mmibias = "{0:.2f}".format(float(MMI[j]-MMIS[j]))
        tout = str('0')
        smf.write(str(SITE[j])+' '+lon+' '+lat+' '+tout+' '+pgaout+' '+mmiout+' '+mmibias+' '+'\n')
    smf.close()
    
    if (props.getdefmodon() == 'yes'):
        deff.close()

###########################################################
#Read ShakeMap Files
###########################################################
for i in range(0, int(numshakemapfiles)):
    shm = xmldir+str(shakemapfiles[i])
    print('Loading ShakeMap file ', str(shm))
    shmoutfile = outdir + str(earthquake) + '_' + str(shakemapalgs[i]) +  '_' + str(props.getgmpe()) + '_' + str(props.getgmice()) +  '_' + str(props.getmmicomp()) + '_out.txt'
    evolfiles.append(shmoutfile)
    algorithm.append(str(shakemapalgs[i]))
    shmf = open(shmoutfile,'w')
    shmf.write('##########################################################################'+'\n')
    shmf.write('#'+'Strong Motion Evolution Results - ShakeMap File'+'\n')
    shmf.write('##########################################################################'+'\n')
    shmf.write('#'+'Earthquake'+' '+earthquake+'\n')
    shmf.write('#'+'ShakeMap File'+' '+str(shakemapfiles[i])+'\n')
    shmf.write('#'+'Algorithm'+' '+str(shakemapalgs[i]) +'\n')
    shmf.write('#'+'Origin Time'+' '+str(origintime)+'\n')
    shmf.write('#'+'Only MMI > '+ props.getmmithreshold() +' used for comparisons'+'\n')
    shmf.write('#'+'MMI and PGA component used ='+' '+props.getmmicomp()+'\n')
    shmf.write('#'+'site'+','+'lon'+','+'lat'+','+'OT+(s)'+','+'pga(cm/s2)'+','+'mmi'+','+'mmibias'+'\n')
    shmf.write('##########################################################################'+'\n')

    for line in open(shm):
        [SMlon, SMlat, SMpga, SMmmi] = MSG_xmlreader.shakemapreader(xmldir+str(line.strip()))

        SMlon = numpy.asarray(SMlon)
        SMlat = numpy.asarray(SMlat)
        SMpga = numpy.asarray(SMpga)
        SMmmi = numpy.asarray(SMmmi)

        for j in range(0, len(LATS)):
            [xsm,ysm]=GMCOMP_coord_tools.ll2utm(SMlon, SMlat, LONS[j], LATS[j])
            [xj,yj]=GMCOMP_coord_tools.ll2utm(LONS[j], LATS[j], LONS[j], LATS[j])
            smdist = numpy.sqrt(numpy.power(xsm-xj,2)+numpy.power(ysm-yj,2))

            a1 = numpy.where(numpy.amin(smdist) == smdist)[0]

            lon = "{0:.4f}".format(float(LONS[j]))
            lat = "{0:.4f}".format(float(LATS[j]))
            pgaout = "{0:.4f}".format(float(SMpga[a1])*980.665/100)
            mmiout = "{0:.2f}".format(float(SMmmi[a1]))
            mmibias = "{0:.2f}".format(float(SMmmi[a1]-MMIS[j]))
            tout = str('0')
            shmf.write(str(SITE[j])+' '+lon+' '+lat+' '+tout+' '+pgaout+' '+mmiout+' '+mmibias+' '+'\n')

    shmf.close()

###########################################################
#Compare Ground Motion Predictions, output to file
###########################################################
MMImat = numpy.nan*numpy.ones([len(LONS),ttot,len(evolfiles)]) #MMI matrix of values
PGAmat = numpy.nan*numpy.ones([len(LONS),ttot,len(evolfiles)]) #PGA matrix of values
for i in range (0, len(evolfiles)):
    with open(str(evolfiles[i]), 'rt') as g:
        rows = (line.split() for line in g)
        for grow in rows:
            if "#" in str(grow[0]):
                pass
            else:
                asite = numpy.where(grow[0] == numpy.asarray(SITE))[0]
                tout = math.floor(float(grow[3]))
                if tout < ttot:
                    MMImat[asite,int(tout),i] = float(grow[5])
                    PGAmat[asite,int(tout),i] = float(grow[4])



fevolout = outdir + str(earthquake) + '_'  + str(props.getgmpe()) + '_' + str(props.getgmice()) +  '_' + str(props.getmmicomp()) + '_' + str(props.getrruporrp()) + '_evolution.txt'
feo = open(fevolout,'w')
feo.write('##########################################################################'+'\n')
feo.write('#'+'Strong Motion Evolution Results - All Algorithms'+'\n')
feo.write('##########################################################################'+'\n')
feo.write('#'+'Earthquake'+' '+earthquake+'\n')
feo.write('#'+'Origin Time'+' '+str(origintime)+'\n')
feo.write('#'+'Only MMI > '+ props.getmmithreshold() +' reported'+'\n')
feo.write('#'+'MMI Warning Threshold is '+ props.getmmiwarnthreshold()+'\n')
feo.write('#'+'MMI and PGA component used ='+' '+props.getmmicomp()+'\n')
for i in range (0, len(evolfiles)):
    feo.write('#'+'Algorithm'+' '+str(i+1)+' '+str(algorithm[i])+'\n')
feo.write('#'+'num algs'+','+'site'+','+'lon'+','+'lat'+','+'OT+(s)'+','+'\a')
for i in range (0, len(evolfiles)):
    feo.write('algorithm'+str(i+1)+','+'mmi'+str(i+1)+','+'pga'+str(i+1)+','+'\a')
feo.write('MMI_ex_time(s)'+','+'maxMMI'+'\n')
feo.write('##########################################################################'+'\n')

for i in range (0, len(LONS)):
    lon = "{0:.4f}".format(float(LONS[i]))
    lat = "{0:.4f}".format(float(LATS[i]))
    maxmmi = "{0:.2f}".format(float(MMIS[i]))
    ammidum = "{0:.0f}".format(float(-1.0))
    ammi = numpy.where(MMImat[i,:,0] >= mmiwarnthreshold)[0]
    
    MMIprev = numpy.zeros([len(evolfiles),1])
    PGAprev = numpy.zeros([len(evolfiles),1])   
    for j in range (0, ttot):

        feo.write(str(len(evolfiles))+' '+str(SITE[i])+' '+lon+' '+lat+' '+str(j)+' ')
        for k in range (0, len(evolfiles)):
            if (math.isnan(MMImat[i,j,k]) == True):
                mmi = "{0:.4f}".format(float(MMIprev[k,0]))
                pga = "{0:.4f}".format(float(PGAprev[k,0]))
            else:
                mmi = "{0:.4f}".format(float(MMImat[i,j,k]))
                pga = "{0:.4f}".format(float(PGAmat[i,j,k]))
                MMIprev[k,0] = MMImat[i,j,k]
                PGAprev[k,0] = PGAmat[i,j,k]
                
            feo.write(str(algorithm[k])+' '+mmi+' '+pga+' ')
        if len(ammi) > 0:
            feo.write(str(ammi[0])+' '+maxmmi+' '+'\n')
        else:
            tstheory = "{0:.0f}".format(float(math.floor(ts[i,0])))
            feo.write(tstheory+' '+maxmmi+' '+'\n')

feo.close()
    
 

###########################################################
#Create GMT plots if turned on
###########################################################
plotson = props.getplotson()
if (plotson == 'yes'):
    GMCOMP_gmtplots.makeplots(algorithm,len(algorithm)-1,props)

