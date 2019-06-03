#!/usr/bin/python
import math
import GMCOMP_coord_tools
import numpy
import xml.etree.ElementTree as ET


def cireader(coreinfofile):
    tree = ET.parse(coreinfofile)
    root = tree.getroot()

    ci_id=[]
    ci_mag=[]
    ci_mag_uncer=[]
    ci_lat=[]
    ci_lat_uncer=[]
    ci_lon=[]
    ci_lon_uncer=[]
    ci_dep=[]
    ci_dep_uncer=[]
    ci_ot=[]
    ci_lik=[]
    
    ci_alg_vers = root.attrib['alg_vers']
    ci_category = root.attrib['category']
    ci_instance = root.attrib['instance']
    ci_message_type = root.attrib['message_type']
    ci_orig_sys = root.attrib['orig_sys']
    ci_timestamp = root.attrib['timestamp']
    ci_version = root.attrib['version']

    for child in root.iter():
        ca = child.tag
        if (ca == 'core_info'):
            ci_id = child.attrib['id']

        if (ca == 'mag'):
            ci_mag = float(child.text)

        if (ca == 'mag_uncer'):
            ci_mag_uncer = float(child.text)
            
        if (ca == 'lat'):
            ci_lat = float(child.text)
            
        if (ca == 'lat_uncer'):
            ci_lat_uncer = float(child.text)
            
        if (ca == 'lon'):
            ci_lon = float(child.text)
            
        if (ca == 'lon_uncer'):
            ci_lon_uncer = float(child.text)
            
        if (ca == 'depth'):
            ci_dep = float(child.text)
            
        if (ca == 'depth_uncer'):
            ci_dep_uncer = float(child.text)

        if (ca == 'orig_time'):
            ci_ot = child.text
            
        if (ca == 'likelihood'):
            ci_lik = float(child.text)
    mess_overview = [ci_alg_vers, ci_category, ci_instance, ci_message_type, ci_orig_sys, ci_timestamp, ci_version]
    mess_details = [ci_id, ci_mag, ci_mag_uncer, ci_lat, ci_lat_uncer, ci_lon, ci_lon_uncer, ci_dep, ci_dep_uncer, ci_ot, ci_lik]
            

    return(mess_overview, mess_details)

def ffreader(fffile):
    tree = ET.parse(fffile)
    root = tree.getroot()
    
    #predefining variables in case one is missing from xml file
    ff_id=[]
    ff_mag=[]
    ff_mag_uncer=[]
    ff_lat=[]
    ff_lat_uncer=[]
    ff_lon=[]
    ff_lon_uncer=[]
    ff_dep=[]
    ff_dep_uncer=[]
    ff_ot=[]
    ff_lik=[]
    lon_trans_unc=[]
    lat_trans_unc=[]
    total_len_unc=[]
    str_unc=[]
    dip_unc=[]

    ff_alg_vers = root.attrib['alg_vers']
    ff_category = root.attrib['category']
    ff_instance = root.attrib['instance']
    ff_message_type = root.attrib['message_type']
    ff_orig_sys = root.attrib['orig_sys']
    ff_timestamp = root.attrib['timestamp']
    ff_version = root.attrib['version']

    LONS=[]
    LATS=[]
    DEPS=[]
    SS=[]
    DS=[]
    SSunc=[]
    DSunc=[]
    fftag=0 #since coreinfo uses lon/lat as well, need to tag once ff message starts
    k=0
    numvert=10 #dummy variable for number of vertices before definition, should be 3 (triangles) or 4 (rectangles)
    
    for child in root.iter():
        ca = child.tag

        if (ca == 'core_info'):
            ff_id = child.attrib['id']
            
        if (ca == 'finite_fault'):
            fftag = 1

        if (ca == 'mag'):
            ff_mag = float(child.text)

        if (ca == 'mag_uncer'):
            ff_mag_uncer = float(child.text)
            
        if (ca == 'lat'):
            ff_lat = float(child.text)
            
        if (ca == 'lat_uncer'):
            ff_lat_uncer = float(child.text)
            
        if (ca == 'lon'):
            ff_lon = float(child.text)
            
        if (ca == 'lon_uncer'):
            ff_lon_uncer = float(child.text)
            
        if (ca == 'depth'):
            ff_dep = float(child.text)
            
        if (ca == 'depth_uncer'):
            ff_dep_uncer = float(child.text)

        if (ca == 'orig_time'):
            ff_ot = child.text
            
        if (ca == 'likelihood'):
            ff_lik = child.text

        if (ca == 'vertices'):
            numvert = int(child.attrib['number'])
            lons = []
            lats = []
            deps = []
            
        if (ca == 'lon' and fftag == 1):
            lons.append(float(child.text))
            
        if (ca == 'lat' and fftag == 1):
            lats.append(float(child.text))
            
        if (ca == 'depth' and fftag == 1):
            units = child.attrib['units']
            if (units == 'km'):
                gain = 1000
            if (units == 'm'):
                gain = 1
            deps.append(float(child.text)*gain)
            k=k+1
        
        if (fftag == 1):
            if (k == numvert):
                k=0
                DEPS.append(deps)
                LATS.append(lats)
                LONS.append(lons)
                
        if (ca == 'ss' and fftag == 1):
            units = child.attrib['units']
            if (units == 'm'):
                gain = 1
            if (units == 'cm'):
                gain = 0.01
            if (units == 'mm'):
                gain = 0.001
            SS.append(float(child.text)*gain)
            
        if (ca == 'ss_uncer' and fftag == 1):
            if (units == 'm'):
                gain = 1
            if (units == 'cm'):
                gain = 0.01
            if (units == 'mm'):
                gain = 0.001
            SSunc.append(float(child.text)*gain)
            
        if (ca == 'ds' and fftag == 1):
            if (units == 'm'):
                gain = 1
            if (units == 'cm'):
                gain = 0.01
            if (units == 'mm'):
                gain = 0.001
            DS.append(float(child.text)*gain)
            
        if (ca == 'ds_uncer' and fftag == 1):
            if (units == 'm'):
                gain = 1
            if (units == 'cm'):
                gain = 0.01
            if (units == 'mm'):
                gain = 0.001
            DSunc.append(float(child.text)*gain)
            
        #Global Uncertainty section           
        if (ca == 'lon_trans'):
            lon_trans_unc = float(child.text)
        if (ca == 'lat_trans'):
            lat_trans_unc = float(child.text)
        if (ca == 'total_len'):
            total_len_unc = float(child.text)
        if (str_unc == 'strike'):
            str_unc = float(child.text)
        if (ca == 'dip'):
            dip_unc = float(child.text)
                       
    DEPS = numpy.asarray(DEPS)
    LATS = numpy.asarray(LATS)
    LONS = numpy.asarray(LONS)
    SS = numpy.asarray(SS)
    DS = numpy.asarray(DS)
    
    mess_overview = [ff_alg_vers, ff_category, ff_instance, ff_message_type, ff_orig_sys, ff_timestamp, ff_version]
    mess_details = [ff_id, ff_mag, ff_mag_uncer, ff_lat, ff_lat_uncer, ff_lon, ff_lon_uncer, ff_dep, ff_dep_uncer, ff_ot, ff_lik, numvert]
    globunc_details = [lon_trans_unc, lat_trans_unc, total_len_unc, str_unc, dip_unc]

    return(mess_overview, mess_details, globunc_details, LONS, LATS, DEPS, SS, DS)

#This is for reading line source files such as from Finder
def ffreader_line(fffile):
    tree = ET.parse(fffile)
    root = tree.getroot()
    
    #predefining variables in case one is missing from xml file
    ff_id=[]
    ff_mag=[]
    ff_mag_uncer=[]
    ff_lat=[]
    ff_lat_uncer=[]
    ff_lon=[]
    ff_lon_uncer=[]
    ff_dep=[]
    ff_dep_uncer=[]
    ff_ot=[]
    ff_lik=[]

    ff_alg_vers = root.attrib['alg_vers']
    ff_category = root.attrib['category']
    ff_instance = root.attrib['instance']
    ff_message_type = root.attrib['message_type']
    ff_orig_sys = root.attrib['orig_sys']
    ff_timestamp = root.attrib['timestamp']
    ff_version = root.attrib['version']

    LONS=[]
    LATS=[]
    DEPS=[]

    fftag=0 #since coreinfo uses lon/lat as well, need to tag once ff message starts
    k=0
    
    for child in root.iter():
        ca = child.tag
        if (ca == 'core_info'):
            ff_id = child.attrib['id']           
        if (ca == 'finite_fault'):
            fftag = 1
        if (ca == 'mag'):
            ff_mag = float(child.text)
        if (ca == 'mag_uncer'):
            ff_mag_uncer = float(child.text)           
        if (ca == 'lat'):
            ff_lat = float(child.text)           
        if (ca == 'lat_uncer'):
            ff_lat_uncer = float(child.text)           
        if (ca == 'lon'):
            ff_lon = float(child.text)            
        if (ca == 'lon_uncer'):
            ff_lon_uncer = float(child.text)            
        if (ca == 'depth'):
            ff_dep = float(child.text)           
        if (ca == 'depth_uncer'):
            ff_dep_uncer = float(child.text)
        if (ca == 'orig_time'):
            ff_ot = child.text           
        if (ca == 'likelihood'):
            ff_lik = child.text

        if (ca == 'vertices'):
            lons = []
            lats = []
            deps = []
            
        if (ca == 'lon' and fftag == 1):
            LONS.append(float(child.text))
            
        if (ca == 'lat' and fftag == 1):
            LATS.append(float(child.text))
            
        if (ca == 'depth' and fftag == 1):
            units = child.attrib['units']
            if (units == 'km'):
                gain = 1000
            if (units == 'm'):
                gain = 1
            DEPS.append(float(child.text)*gain)
        if (ca == 'gm_info'):
            fftag = 2
        
##        if (fftag == 1):
##            DEPS.append(deps)
##            LATS.append(lats)
##            LONS.append(lons)
                       
    DEPS = numpy.asarray(DEPS)
    LATS = numpy.asarray(LATS)
    LONS = numpy.asarray(LONS)
    
    mess_overview = [ff_alg_vers, ff_category, ff_instance, ff_message_type, ff_orig_sys, ff_timestamp, ff_version]
    mess_details = [ff_id, ff_mag, ff_mag_uncer, ff_lat, ff_lat_uncer, ff_lon, ff_lon_uncer, ff_dep, ff_dep_uncer, ff_ot, ff_lik]


    return(mess_overview, mess_details, LONS, LATS, DEPS)

def eqinfo2gmreader(eqinfofile):
    tree = ET.parse(eqinfofile)
    root = tree.getroot()

    eqi_id=[]
    eqi_mag=[]
    eqi_mag_uncer=[]
    eqi_lat=[]
    eqi_lat_uncer=[]
    eqi_lon=[]
    eqi_lon_uncer=[]
    eqi_dep=[]
    eqi_dep_uncer=[]
    eqi_ot=[]
    eqi_lik=[]
    
    eqi_alg_vers = root.attrib['alg_vers']
    eqi_category = root.attrib['category']
    eqi_instance = root.attrib['instance']
    eqi_message_type = root.attrib['message_type']
    eqi_orig_sys = root.attrib['orig_sys']
    eqi_timestamp = root.attrib['timestamp']
    eqi_version = root.attrib['version']

    map_lon=list()
    map_lat=list()
    map_mmi=list()
    map_pga=list()

    contour_lon=list()
    contour_lat=list()
    contour_mmi=list()
    contour_pga=list()

    for child in root.iter():
        ca = child.tag
        if (ca == 'core_info'):
            eqi_id = child.attrib['id']

        if (ca == 'mag'):
            eqi_mag = float(child.text)

        if (ca == 'mag_uncer'):
            eqi_mag_uncer = float(child.text)
            
        if (ca == 'lat'):
            eqi_lat = float(child.text)
            
        if (ca == 'lat_uncer'):
            eqi_lat_uncer = float(child.text)
            
        if (ca == 'lon'):
            eqi_lon = float(child.text)
            
        if (ca == 'lon_uncer'):
            eqi_lon_uncer = float(child.text)
            
        if (ca == 'depth'):
            eqi_dep = float(child.text)
            
        if (ca == 'depth_uncer'):
            eqi_dep_uncer = float(child.text)

        if (ca == 'orig_time'):
            eqi_ot = child.text
            
        if (ca == 'likelihood'):
            eqi_lik = float(child.text) 
                   
        if (ca == 'grid_data'):
            line = child.text
            eqinfotype='map'
            data = line.splitlines()
            dataarray = numpy.asarray(data)
            len1=len(dataarray)
            for i in range (0, int(len1)):
                if (len(dataarray[i]) > 0):
                    dataline = dataarray[i].split()
                    map_lon.append(float(dataline[0]))
                    map_lat.append(float(dataline[1]))
                    map_pga.append(float(dataline[2]))
                    map_mmi.append(float(dataline[4]))

    if (eqinfotype == 'map'):
        gmout = [map_lon, map_lat, map_pga, map_mmi]
            
    mess_overview = [eqi_alg_vers, eqi_category, eqi_instance, eqi_message_type, eqi_orig_sys, eqi_timestamp, eqi_version]
    mess_details = [eqi_id, eqi_mag, eqi_mag_uncer, eqi_lat, eqi_lat_uncer, eqi_lon, eqi_lon_uncer, eqi_dep, eqi_dep_uncer, eqi_ot, eqi_lik]

    return(mess_overview, mess_details, gmout, eqinfotype)


def shakemapreader(shakemapfile):
    tree = ET.parse(shakemapfile)
    root = tree.getroot()

    sm_lon=list()
    sm_lat=list()
    sm_mmi=list()
    sm_pga=list()

    for child in root.iter():
        ca = child.tag
        if (ca == '{http://earthquake.usgs.gov/eqcenter/shakemap}grid_data'):
            line = child.text
            data = line.splitlines()
            dataarray = numpy.asarray(data)
            len1=len(dataarray)

            for i in range (0, int(len1)):
                if (len(dataarray[i]) > 0):
                    dataline = dataarray[i].split()
                    sm_lon.append(float(dataline[0]))
                    sm_lat.append(float(dataline[1]))
                    sm_pga.append(float(dataline[2]))
                    sm_mmi.append(float(dataline[4]))

    return(sm_lon, sm_lat, sm_pga, sm_mmi)

##[mess_overview, mess_details, LONS, LATS, DEPS] = ffreader_line('xmlfiles/finder_napa/south_napa_finder_message_version17.xml')

##[mess_overview, mess_details, gmout] = eqinfo2gmreader('xmlfiles/eqinfo2gm_napa/south_napa_eqinfo2gm_map_message_version0.xml')
##print (gmout[3])   
#shakemapreader('xmlfiles/shakemap_napa/shakemap_napa.xml')
