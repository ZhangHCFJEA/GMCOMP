#!/usr/bin/python
import math
import numpy
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
#######################################################
#This subroutine computes the theoretical P and S wave
#traveltime from an earthquake source to a receiver.
#Output is tp and ts
#Written by Brendan Crowell, University of Washington, 2017
#By default, I use the ak135 1-D velocity model
#To use other models, see availability at
#https://docs.obspy.org/packages/obspy.taup.html
#and change ak135 to model of your choice below
#
#This code will take in a vector of receivers, be sure
#to put into proper numpy format
#######################################################
#Variable Definitions
#######################################################
#eqlat is the latitude of the earthquake
#eqlon is the longitude of the earthquake
#eqdep is the depth of the earthquake, in km
#reclat is the latitude of the receiver
#reclon is the longitude of the receiver
#######################################################
def tt(eqlon,eqlat,eqdep,reclon,reclat):
    model = TauPyModel(model="ak135")
    [lenrec,aa] = numpy.shape(reclon)
    tp = 999999.9*numpy.ones([lenrec,1])
    ts = 999999.9*numpy.ones([lenrec,1])
    for i in range (0, lenrec):
        thetadist = locations2degrees(reclat[i], reclon[i], eqlat, eqlon)

        if thetadist < 7.0:
            parrival = model.get_travel_times(source_depth_in_km=eqdep,distance_in_degree=thetadist,phase_list=["p"])
            sarrival = model.get_travel_times(source_depth_in_km=eqdep,distance_in_degree=thetadist,phase_list=["s"])
            if len(parrival) > 0:
                parr = parrival[0]
                tp[i,0] = parr.time
            if len(parrival) == 0:
                parrival2 = model.get_travel_times(source_depth_in_km=eqdep,distance_in_degree=thetadist,phase_list=["P"])
                if len(parrival2) > 0:
                    parr = parrival2[0]
                    tp[i,0] = parr.time

            if len(sarrival) > 0:
                sarr = sarrival[0]
                ts[i,0] = sarr.time
            if len(sarrival) == 0:
                sarrival2 = model.get_travel_times(source_depth_in_km=eqdep,distance_in_degree=thetadist,phase_list=["S"])
                if len(sarrival2) > 0:
                    sarr = sarrival2[0]
                    ts[i,0] = sarr.time
    return(tp,ts)


