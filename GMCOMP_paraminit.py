#!/usr/bin/python

class Properties:
        def __init__(self,propfilename):
                self.dict={}
                infile=open(propfilename,'r')
                for line in infile:
                        if (line) and (line[0]!='#') and ('=' in line):
                                (key,val)=line.split('=')
                                self.dict[key]=val.strip()
                return

        ###########################################################
        #Event properties and basic logistics
        ###########################################################

        def getttot(self):
                if 'ttot' in self.dict:
                        return self.dict['ttot']

        def getotyr(self):
                if 'ot_yr' in self.dict:
                        return self.dict['ot_yr']
        
        def getotmo(self):
                if 'ot_mo' in self.dict:
                        return self.dict['ot_mo']
        
        def getotdy(self):
                if 'ot_dy' in self.dict:
                        return self.dict['ot_dy']
        
        def getothr(self):
                if 'ot_hr' in self.dict:
                        return self.dict['ot_hr']
        
        def getotmn(self):
                if 'ot_mn' in self.dict:
                        return self.dict['ot_mn']
        
        def getotsc(self):
                if 'ot_sc' in self.dict:
                        return self.dict['ot_sc']

        def getaccdatadir(self): 
                if 'accdatadir' in self.dict:
                        return self.dict['accdatadir']

        def getoutdir(self): 
                if 'outputdir' in self.dict:
                        return self.dict['outputdir']

        def getxmldir(self): 
                if 'xmldir' in self.dict:
                        return self.dict['xmldir']
                
        def getaccchanfile(self): 
                if 'accchanfile' in self.dict:
                        return self.dict['accchanfile']

        def geteq(self): 
                if 'eq' in self.dict:
                        return self.dict['eq']
                
        def geteqlon(self): 
                if 'eqlon' in self.dict:
                        return self.dict['eqlon']
                
        def geteqlat(self): 
                if 'eqlat' in self.dict:
                        return self.dict['eqlat']
                
        def geteqdep(self): 
                if 'eqdep' in self.dict:
                        return self.dict['eqdep']

        def getrunmseed(self): 
                if 'runmseed' in self.dict:
                        return self.dict['runmseed']

        ###########################################################
        #Ground motion model properties
        ###########################################################
               
        def getmmicomp(self): 
                if 'mmicomp' in self.dict:
                        return self.dict['mmicomp']

        def getgmpe(self): 
                if 'gmpe' in self.dict:
                        return self.dict['gmpe']

        def getgmice(self): 
                if 'gmice' in self.dict:
                        return self.dict['gmice']

        def getvs30file(self): 
                if 'vs30file' in self.dict:
                        return self.dict['vs30file']

        def getmmithreshold(self): 
                if 'mmithreshold' in self.dict:
                        return self.dict['mmithreshold']
                
        def getmmiwarnthreshold(self): 
                if 'mmiwarnthreshold' in self.dict:
                        return self.dict['mmiwarnthreshold']

        def getrruporrp(self): 
                if 'rruporrp' in self.dict:
                        return self.dict['rruporrp']

        def getpvalue(self): 
                if 'p' in self.dict:
                        return self.dict['p']
                
        ###########################################################
        #Finite fault message properties
        ###########################################################

        def getnumffmessages(self): 
                if 'numffmessages' in self.dict:
                        return self.dict['numffmessages']
                
        def getffmessages(self):
                if 'numffmessages' in self.dict:
                        nmess = int(self.dict['numffmessages'])
                ffmessages = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'ffmessage' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                ffmessages.append(fname)
                return(ffmessages)

        def getffalgs(self):
                if 'numffmessages' in self.dict:
                        nmess = int(self.dict['numffmessages'])
                ffalgs = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'ffalg' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                ffalgs.append(fname)
                return(ffalgs)
        
        ###########################################################
        #Line source message properties
        ###########################################################

        def getnumlsmessages(self): 
                if 'numlsmessages' in self.dict:
                        return self.dict['numlsmessages']
                
        def getlsmessages(self):
                if 'numlsmessages' in self.dict:
                        nmess = int(self.dict['numlsmessages'])
                lsmessages = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'lsmessage' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                lsmessages.append(fname)
                return(lsmessages)

        def getlsalgs(self):
                if 'numlsmessages' in self.dict:
                        nmess = int(self.dict['numlsmessages'])
                lsalgs = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'lsalg' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                lsalgs.append(fname)
                return(lsalgs)

        ###########################################################
        #Coreinfo message properties (for just point source algs)
        ###########################################################

        def getnumcoremessages(self): 
                if 'numcoremessages' in self.dict:
                        return self.dict['numcoremessages']
                
        def getcoremessages(self):
                if 'numcoremessages' in self.dict:
                        nmess = int(self.dict['numcoremessages'])
                coremessages = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'coremessage' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                coremessages.append(fname)
                return(coremessages)

        def getcorealgs(self):
                if 'numcoremessages' in self.dict:
                        nmess = int(self.dict['numcoremessages'])
                corealgs = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'corealg' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                corealgs.append(fname)
                return(corealgs)
        
        ###########################################################
        #Eqinfo2gm message properties
        ###########################################################
        
        def getnumeqinfo2gmmessages(self): 
                if 'numeqinfo2gmmessages' in self.dict:
                        return self.dict['numeqinfo2gmmessages']
                
        def geteqinfo2gmmessages(self):
                if 'numeqinfo2gmmessages' in self.dict:
                        nmess = int(self.dict['numeqinfo2gmmessages'])
                eqinfo2gmmessages = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'eqinfo2gmmessage' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                eqinfo2gmmessages.append(fname)
                return(eqinfo2gmmessages)

        def geteqinfo2gmalgs(self):
                if 'numeqinfo2gmmessages' in self.dict:
                        nmess = int(self.dict['numeqinfo2gmmessages'])
                eqinfoalgs = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'eqinfo2gmalg' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                eqinfoalgs.append(fname)
                return(eqinfoalgs)

        ###########################################################
        #SRCMOD fsp file properties. Separate folder for SRCMOD
        ###########################################################

        def getnumsrcmodfiles(self): 
                if 'numsrcmodfiles' in self.dict:
                        return self.dict['numsrcmodfiles']
                
        def getsrcmodfiles(self):
                if 'numsrcmodfiles' in self.dict:
                        nmess = int(self.dict['numsrcmodfiles'])
                srcmodfiles = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'srcmodfile' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                srcmodfiles.append(fname)
                return(srcmodfiles)

        def getsrcmodalgs(self):
                if 'numsrcmodfiles' in self.dict:
                        nmess = int(self.dict['numsrcmodfiles'])
                srcmodalgs = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'srcmodalg' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                srcmodalgs.append(fname)
                return(srcmodalgs)

        def getsrcmoddir(self): 
                if 'srcmoddir' in self.dict:
                        return self.dict['srcmoddir']
                
        ###########################################################
        #ShakeMap XML File Properties
        ###########################################################

        def getnumshakemapfiles(self): 
                if 'numshakemapfiles' in self.dict:
                        return self.dict['numshakemapfiles']
                
        def getshakemapfiles(self):
                if 'numshakemapfiles' in self.dict:
                        nmess = int(self.dict['numshakemapfiles'])
                shakemapfiles = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'shakemapfile' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                shakemapfiles.append(fname)
                return(shakemapfiles)
        
        def getshakemapalgs(self):
                if 'numshakemapfiles' in self.dict:
                        nmess = int(self.dict['numshakemapfiles'])
                shakemapalgs = list()
                for i in range(0,int(nmess)+1):
                        fmess = 'shakemapalg' + str(i)
                        if fmess in self.dict:
                                fname = self.dict[fmess]
                                shakemapalgs.append(fname)
                return(shakemapalgs)
        
        ###########################################################
        #Deformation Forward Modeling section
        ###########################################################
        
        def getdefmodon(self): 
                if 'defmodon' in self.dict:
                        return self.dict['defmodon']

        def getgridlat(self): 
                if 'gridlat' in self.dict:
                        return self.dict['gridlat']
                
        def getgridlon(self): 
                if 'gridlon' in self.dict:
                        return self.dict['gridlon']
                
        def getgridlatcen(self): 
                if 'gridlatcen' in self.dict:
                        return self.dict['gridlatcen']
                
        def getgridloncen(self): 
                if 'gridloncen' in self.dict:
                        return self.dict['gridloncen']
                
        def getgridnum(self): 
                if 'gridnum' in self.dict:
                        return self.dict['gridnum']

        ###########################################################
        #Plotting Parameters
        ###########################################################

        def getplotson(self): 
                if 'plotson' in self.dict:
                        return self.dict['plotson']

        def getgmtpath(self): 
                if 'gmtpath' in self.dict:
                        return self.dict['gmtpath']

        def getgmtversion(self): 
                if 'gmtversion' in self.dict:
                        return self.dict['gmtversion']

        ###########################################################
        #Video Plotting Parameters
        ###########################################################
                
        def getvidson(self): 
                if 'vidson' in self.dict:
                        return self.dict['vidson']

        def getvidoutdir(self): 
                if 'vidoutdir' in self.dict:
                        return self.dict['vidoutdir']
                
        def getvidlength(self): 
                if 'vidlength' in self.dict:
                        return self.dict['vidlength']




