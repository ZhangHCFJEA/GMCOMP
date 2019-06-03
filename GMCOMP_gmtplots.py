#!/usr/bin/python
import math
import numpy
import os
from GMCOMP_paraminit import Properties
###########################################################
#GMCOMP_gmtplots.py
#This takes in the main mmi evolution file with all algorithms
#and creates several different types of plots
#Last modified, May 31, 2019, Brendan Crowell, UW
###########################################################

def makeplots(alg,numalgs,props):
    eqlon = float(props.geteqlon())
    eqlat = float(props.geteqlat())
    eqdep = float(props.geteqdep())
    evolfile = props.getoutdir()+props.geteq()+'_'+props.getgmpe()+'_'+props.getgmice()+'_'+props.getmmicomp()+'_'+props.getrruporrp()+'_evolution.txt'

    analyzemmidata(evolfile, numalgs, props)

    if (props.getgmtversion() == '5'):
        #MMI Bias figure
        outfile1 = props.getoutdir() + props.geteq() + '_mmibias_density_' + props.getgmpe() + '_' + props.getgmice()+'_'+props.getmmicomp()+'_'+props.getrruporrp()
        print ('Creating figure', outfile1+'.eps')
        plotfile1 = props.getoutdir()+'plot_mmibias'
        f1 = open(plotfile1,'w')
        f1.write('#!/bin/bash'+'\n')
        f1.write('export PATH='+props.getgmtpath()+':$PATH'+'\n')
        f1.write('  NAME='+outfile1+'\n')
        f1.write('rm $NAME.eps'+'\n')
        f1.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
        f1.write('gmt set PS_PAGE_ORIENTATION landscape'+'\n')
        f1.write('gmt set PS_MEDIA A2'+'\n')
        f1.write('gmt set MAP_TITLE_OFFSET 0.25c'+'\n')
        f1.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
        f1.write('gmt set FONT_TITLE 16p,Helvetica'+'\n')
        f1.write('makecpt -Chaxby -T0/0.5/0.01 -D -I > tmp.cpt'+'\n')
        f1.write('psbasemap -R-50/50/-5/5 -JX-6i/3i -Ba25f5:"Warning Time (s)":/a5f1:"MMI Bias (predict-true)"::."'+props.geteq()+', '+alg[1]+'":WeSn -Y13i -K > $NAME.eps'+'\n')
        f1.write('psxy '+props.getoutdir()+props.geteq()+'_mmi_density0.txt -Ctmp.cpt -R -J -P -L -O -K >> $NAME.eps'+'\n')
        n=0
        for i in range(1,numalgs):
            if ((i-n*3)%3 == 0):
                f1.write('psxy '+props.getoutdir()+props.geteq()+'_mmi_density'+str(i)+'.txt -Ctmp.cpt -Ba25f5:"Warning Time (s)":/a5f1:"MMI Bias (predict-true)"::."'+props.geteq()+', '+alg[i+1]+'":WeSn -R -J -P -L -O -K -Y-4.25i -X-14i >> $NAME.eps'+'\n')
                n=n+1
            if ((i-n*3)%3 == 1):
                f1.write('psxy '+props.getoutdir()+props.geteq()+'_mmi_density'+str(i)+'.txt -Ctmp.cpt -Ba25f5:"Warning Time (s)":/a5f1:"MMI Bias (predict-true)"::."'+props.geteq()+', '+alg[i+1]+'":WeSn -R -J -P -L -O -K -X7i >> $NAME.eps'+'\n')
            if ((i-n*3)%3 == 2):
                f1.write('psxy '+props.getoutdir()+props.geteq()+'_mmi_density'+str(i)+'.txt -Ctmp.cpt -Ba25f5:"Warning Time (s)":/a5f1:"MMI Bias (predict-true)"::."'+props.geteq()+', '+alg[i+1]+'":WeSn -R -J -P -L -O -K -X7i >> $NAME.eps'+'\n')

        if (numalgs%3 == 1):
            xoff='7'
            yoff='3.5'
        if (numalgs%3 == 2):
            xoff='0'
            yoff='-1'
        if (numalgs%3 == 0):
            xoff='-7'
            yoff='-1'
        f1.write('psscale -O -D3i/0i/6i/0.5ih  -B0.5f0.1:"Probability":/:: -Ctmp.cpt -X'+xoff+'i -Y'+yoff+'i >> $NAME.eps'+'\n') 
        #f1.write('psconvert $NAME.eps -A -TE'+'\n')
        f1.write('exit 0'+'\n')
        f1.close()

        os.system('chmod 777' + ' ' + plotfile1)
        os.system('./'+plotfile1)

        #False positive, CDF figure
        outfile2 = props.getoutdir() + props.geteq() + '_cdf_fp_' + props.getgmpe() + '_' + props.getgmice()+'_'+props.getmmicomp()+'_'+props.getrruporrp()
        plotfile2 = props.getoutdir()+'plot_cdf_fp'
        print ('Creating figure', outfile2+'.eps')
        f2 = open(plotfile2,'w')
        f2.write('#!/bin/bash'+'\n')
        f2.write('export PATH='+props.getgmtpath()+':$PATH'+'\n')
        f2.write('  NAME='+outfile2+'\n')
        f2.write('rm $NAME.eps'+'\n')
        f2.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
        f2.write('gmt set PS_PAGE_ORIENTATION landscape'+'\n')
        f2.write('gmt set PS_MEDIA A2'+'\n')
        f2.write('gmt set MAP_TITLE_OFFSET 0.25c'+'\n')
        f2.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
        f2.write('gmt set FONT_TITLE 16p,Helvetica'+'\n')
        f2.write('psbasemap -R-50/50/0/1 -JX-6i/3i -Ba25f5:"Warning Time (s)":/a0.5f0.1:"CDF of Warned Stations"::."'+props.geteq()+', '+alg[1]+'":WeSn -Y13i -K > $NAME.eps'+'\n')
        f2.write('awk \'{print $1, $2}\' '+props.getoutdir()+props.geteq()+'_mmi_cdf0.txt |psxy -W3 -R -J -O -K >> $NAME.eps'+'\n')
        f2.write('gmt set FONT_LABEL 16p,Helvetica,red'+'\n')
        f2.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,red'+'\n')
        f2.write('awk \'{print $1, $3}\' '+props.getoutdir()+props.geteq()+'_mmi_cdf0.txt |psxy -W3,red -R -J -O -K -Ba25f5::/a0.5f0.1:"False Positive Rate":E >> $NAME.eps'+'\n')

        n=0
        for i in range(1,numalgs):
            if ((i-n*3)%3 == 0):
                f2.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
                f2.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
                f2.write('awk \'{print $1, $2}\' '+props.getoutdir()+props.geteq()+'_mmi_cdf'+str(i)+'.txt |psxy -W3 -R -J -O -K -Ba25f5:"Warning Time (s)":/a0.5f0.1:"CDF of Warned Stations"::."'+props.geteq()+', '+alg[i+1]+'":WeSn -Y-4.25i -X-15i >> $NAME.eps'+'\n')
                f2.write('gmt set FONT_LABEL 16p,Helvetica,red'+'\n')
                f2.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,red'+'\n')
                f2.write('awk \'{print $1, $3}\' '+props.getoutdir()+props.geteq()+'_mmi_cdf'+str(i)+'.txt |psxy -W3,red -R -J -O -K -Ba25f5::/a0.5f0.1:"False Positive Rate":E >> $NAME.eps'+'\n')
                n=n+1
            if ((i-n*3)%3 == 1):
                f2.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
                f2.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
                f2.write('awk \'{print $1, $2}\' '+props.getoutdir()+props.geteq()+'_mmi_cdf'+str(i)+'.txt |psxy -W3 -R -J -O -K -Ba25f5:"Warning Time (s)":/a0.5f0.1:"CDF of Warned Stations"::."'+props.geteq()+', '+alg[i+1]+'":WeSn -X7.5i >> $NAME.eps'+'\n')
                f2.write('gmt set FONT_LABEL 16p,Helvetica,red'+'\n')
                f2.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,red'+'\n')
                f2.write('awk \'{print $1, $3}\' '+props.getoutdir()+props.geteq()+'_mmi_cdf'+str(i)+'.txt |psxy -W3,red -R -J -O -K -Ba25f5::/a0.5f0.1:"False Positive Rate":E >> $NAME.eps'+'\n')
            if ((i-n*3)%3 == 2):
                f2.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
                f2.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
                f2.write('awk \'{print $1, $2}\' '+props.getoutdir()+props.geteq()+'_mmi_cdf'+str(i)+'.txt |psxy -W3 -R -J -O -K -Ba25f5:"Warning Time (s)":/a0.5f0.1:"CDF of Warned Stations"::."'+props.geteq()+', '+alg[i+1]+'":WeSn -X7.5i >> $NAME.eps'+'\n')
                f2.write('gmt set FONT_LABEL 16p,Helvetica,red'+'\n')
                f2.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,red'+'\n')
                f2.write('awk \'{print $1, $3}\' '+props.getoutdir()+props.geteq()+'_mmi_cdf'+str(i)+'.txt |psxy -W3,red -R -J -O -K -Ba25f5::/a0.5f0.1:"False Positive Rate":E >> $NAME.eps'+'\n')
        f2.write('exit 0'+'\n')
        f2.close()

        os.system('chmod 777' + ' ' + plotfile2)
        os.system('./'+plotfile2)


        #quad plot
        outfile3 = props.getoutdir() + props.geteq() + '_quad_' + props.getgmpe() + '_' + props.getgmice()+'_'+props.getmmicomp()+'_'+props.getrruporrp()
        plotfile3 = props.getoutdir()+'plot_quad'
        print ('Creating figure', outfile3+'.eps')
        f3 = open(plotfile3,'w')
        f3.write('#!/bin/bash'+'\n')
        f3.write('export PATH='+props.getgmtpath()+':$PATH'+'\n')
        f3.write('  NAME='+outfile3+'\n')
        f3.write('rm $NAME.eps'+'\n')
        f3.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
        f3.write('gmt set PS_PAGE_ORIENTATION landscape'+'\n')
        f3.write('gmt set PS_MEDIA A2'+'\n')
        f3.write('gmt set MAP_TITLE_OFFSET 0.25c'+'\n')
        f3.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
        f3.write('gmt set FONT_TITLE 16p,Helvetica'+'\n')
        f3.write('makecpt -Chaxby -T0/40/5 > tmp.cpt'+'\n')
        f3.write('psbasemap -R0/10/0/10 -JX4i/4i -Ba2g2:"Observed MMI":/a2g2:"Predicted MMI"::."'+props.geteq()+', '+alg[1]+'":WeSn -Y12i -K > $NAME.eps'+'\n')
        f3.write('awk \'{print $5, $4, -$6}\' '+props.getoutdir()+props.geteq()+'_mmi_quad0.txt |psxy -Sc0.2 -W -R -J -O -K -Ctmp.cpt >> $NAME.eps'+'\n')
        f3.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
        f3.write('0 '+props.getmmiwarnthreshold()+'\n')
        f3.write('10 '+props.getmmiwarnthreshold()+'\n')
        f3.write('EOF'+'\n')
        f3.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
        f3.write(props.getmmiwarnthreshold()+' 0'+'\n')
        f3.write(props.getmmiwarnthreshold()+' 10'+'\n')
        f3.write('EOF'+'\n')
        f3.write('echo "0.5 0.5 14,1,black 0 LM TN" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
        f3.write('echo "0.5 9.5 14,1,black 0 LM FP" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
        f3.write('echo "9.0 9.5 14,1,black 0 LM TP" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
        f3.write('echo "9.0 0.5 14,1,black 0 LM FN" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')

        n=0
        for i in range(1,numalgs):
            if ((i-n*3)%3 == 0):
                f3.write('awk \'{print $5, $4, -$6}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -Ctmp.cpt -R -J -O -K -Ba2g2:"Observed MMI":/a2g2:"Predicted MMI"::."'+props.geteq()+', '+alg[i+1]+'":WeSn  -Y-5.5i -X-11i >> $NAME.eps'+'\n')
                f3.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f3.write('0 '+props.getmmiwarnthreshold()+'\n')
                f3.write('10 '+props.getmmiwarnthreshold()+'\n')
                f3.write('EOF'+'\n')
                f3.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f3.write(props.getmmiwarnthreshold()+' 0'+'\n')
                f3.write(props.getmmiwarnthreshold()+' 10'+'\n')
                f3.write('EOF'+'\n')
                f3.write('echo "0.5 0.5 14,1,black 0 LM TN" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                f3.write('echo "0.5 9.5 14,1,black 0 LM FP" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                f3.write('echo "9.0 9.5 14,1,black 0 LM TP" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                f3.write('echo "9.0 0.5 14,1,black 0 LM FN" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                n=n+1
            if ((i-n*3)%3 == 1):
                f3.write('awk \'{print $5, $4, -$6}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -Ctmp.cpt -R -J -O -K -Ba2g2:"Observed MMI":/a2g2:"Predicted MMI"::."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X5.5i >> $NAME.eps'+'\n')
                f3.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f3.write('0 '+props.getmmiwarnthreshold()+'\n')
                f3.write('10 '+props.getmmiwarnthreshold()+'\n')
                f3.write('EOF'+'\n')
                f3.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f3.write(props.getmmiwarnthreshold()+' 0'+'\n')
                f3.write(props.getmmiwarnthreshold()+' 10'+'\n')
                f3.write('EOF'+'\n')
                f3.write('echo "0.5 0.5 14,1,black 0 LM TN" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                f3.write('echo "0.5 9.5 14,1,black 0 LM FP" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                f3.write('echo "9.0 9.5 14,1,black 0 LM TP" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                f3.write('echo "9.0 0.5 14,1,black 0 LM FN" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
            if ((i-n*3)%3 == 2):
                f3.write('awk \'{print $5, $4, -$6}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -Ctmp.cpt -R -J -O -K -Ba2g2:"Observed MMI":/a2g2:"Predicted MMI"::."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X5.5i >> $NAME.eps'+'\n')
                f3.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f3.write('0 '+props.getmmiwarnthreshold()+'\n')
                f3.write('10 '+props.getmmiwarnthreshold()+'\n')
                f3.write('EOF'+'\n')
                f3.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f3.write(props.getmmiwarnthreshold()+' 0'+'\n')
                f3.write(props.getmmiwarnthreshold()+' 10'+'\n')
                f3.write('EOF'+'\n')
                f3.write('echo "0.5 0.5 14,1,black 0 LM TN" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                f3.write('echo "0.5 9.5 14,1,black 0 LM FP" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                f3.write('echo "9.0 9.5 14,1,black 0 LM TP" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')
                f3.write('echo "9.0 0.5 14,1,black 0 LM FN" | pstext -R -J -O -K -F+f+a+j >> $NAME.eps'+'\n')

        f3.write('psscale -O -D4i/2i/4i/0.5i  -B10f5:"Warning Time (s)":/:: -Ctmp.cpt -X0.5i >> $NAME.eps'+'\n') 
        f3.write('exit 0'+'\n')
        f3.close()

        os.system('chmod 777' + ' ' + plotfile3)
        os.system('./'+plotfile3)

        #Warning Time vs observed MMI
        outfile4 = props.getoutdir() + props.geteq() + '_warntime_obsmmi_' + props.getgmpe() + '_' + props.getgmice()+'_'+props.getmmicomp()+'_'+props.getrruporrp()
        plotfile4 = props.getoutdir()+'plot_warntime'
        print ('Creating figure', outfile4+'.eps')
        f4 = open(plotfile4,'w')
        f4.write('#!/bin/bash'+'\n')
        f4.write('export PATH='+props.getgmtpath()+':$PATH'+'\n')
        f4.write('  NAME='+outfile4+'\n')
        f4.write('rm $NAME.eps'+'\n')
        f4.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
        f4.write('gmt set PS_PAGE_ORIENTATION landscape'+'\n')
        f4.write('gmt set PS_MEDIA A2'+'\n')
        f4.write('gmt set MAP_TITLE_OFFSET 0.25c'+'\n')
        f4.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
        f4.write('gmt set FONT_TITLE 16p,Helvetica'+'\n')
        f4.write('makecpt -Cpolar -T-3/3/0.5 > tmp.cpt'+'\n')
        f4.write('psbasemap -R-10/75/0/10 -JX4i/4i -Ba10g5:"Warning Time (s)":/a2g2:"Observed MMI"::."'+props.geteq()+', '+alg[1]+'":WeSn -Y12i -K > $NAME.eps'+'\n')
        f4.write('awk \'{print -$6, $5, $4-$5}\' '+props.getoutdir()+props.geteq()+'_mmi_quad0.txt |psxy -Sc0.2 -W -R -J -O -K -Ctmp.cpt >> $NAME.eps'+'\n')
        f4.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
        f4.write('-10 '+props.getmmiwarnthreshold()+'\n')
        f4.write('75 '+props.getmmiwarnthreshold()+'\n')
        f4.write('EOF'+'\n')

        n=0
        for i in range(1,numalgs):
            if ((i-n*3)%3 == 0):
                f4.write('awk \'{print -$6, $5, $4-$5}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -R -J -O -K -Ctmp.cpt -Ba10g5:"Warning Time (s)":/a2g2:"Observed MMI"::."'+props.geteq()+', '+alg[i+1]+'":WeSn  -Y-5.5i -X-11i >> $NAME.eps'+'\n')
                f4.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f4.write('-10 '+props.getmmiwarnthreshold()+'\n')
                f4.write('75 '+props.getmmiwarnthreshold()+'\n')
                f4.write('EOF'+'\n')
                n=n+1
            if ((i-n*3)%3 == 1):
                f4.write('awk \'{print -$6, $5, $4-$5}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -R -J -O -K -Ctmp.cpt -Ba10g5:"Warning Time (s)":/a2g2:"Observed MMI"::."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X5.5i >> $NAME.eps'+'\n')
                f4.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f4.write('-10 '+props.getmmiwarnthreshold()+'\n')
                f4.write('75 '+props.getmmiwarnthreshold()+'\n')
                f4.write('EOF'+'\n')
            if ((i-n*3)%3 == 2):
                f4.write('awk \'{print -$6, $5, $4-$5}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -R -J -O -K -Ctmp.cpt -Ba10g5:"Warning Time (s)":/a2g2:"Observed MMI"::."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X5.5i >> $NAME.eps'+'\n')
                f4.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f4.write('-10 '+props.getmmiwarnthreshold()+'\n')
                f4.write('75 '+props.getmmiwarnthreshold()+'\n')
                f4.write('EOF'+'\n')
        f4.write('psscale -O -D4i/2i/4i/0.5i  -B1f1:"MMI Bias (Pred.-Obs.)":/:: -Ctmp.cpt -X0.5i >> $NAME.eps'+'\n') 
        f4.write('exit 0'+'\n')
        f4.close()

        os.system('chmod 777' + ' ' + plotfile4)
        os.system('./'+plotfile4)

        #Cost Ratio
        outfile5 = props.getoutdir() + props.geteq() + '_costratio_' + props.getgmpe() + '_' + props.getgmice()+'_'+props.getmmicomp()+'_'+props.getrruporrp()
        plotfile5 = props.getoutdir()+'plot_costratio'
        print ('Creating figure', outfile5+'.eps')
        f5 = open(plotfile5,'w')
        f5.write('#!/bin/bash'+'\n')
        f5.write('export PATH='+props.getgmtpath()+':$PATH'+'\n')
        f5.write('  NAME='+outfile5+'\n')
        f5.write('rm $NAME.eps'+'\n')
        f5.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
        f5.write('gmt set PS_PAGE_ORIENTATION landscape'+'\n')
        f5.write('gmt set PS_MEDIA A2'+'\n')
        f5.write('gmt set MAP_TITLE_OFFSET 0.25c'+'\n')
        f5.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
        f5.write('gmt set FONT_TITLE 16p,Helvetica'+'\n')
        f5.write('psbasemap -R0/20/-0.5/1 -JX5i/4i -Ba5g5:"Cost Ratio, r":/a0.2g0.2:"Cost Savings Performance Metric, Q"::."'+props.geteq()+', '+alg[1]+'":WeSn -Y12i -K > $NAME.eps'+'\n')
        f5.write('awk \'{print $1, $2}\' '+props.getoutdir()+props.geteq()+'_mmi_costratio0.txt |psxy -W3 -R -J -O -K >> $NAME.eps'+'\n')

        n=0
        for i in range(1,numalgs):
            if ((i-n*3)%3 == 0):
                f5.write('awk \'{print $1, $2}\' '+props.getoutdir()+props.geteq()+'_mmi_costratio'+str(i)+'.txt |psxy -W3 -R -J -O -K -Ba5g5:"Cost Ratio, r":/a0.2g0.2:"Cost Savings Performance Metric, Q"::."'+props.geteq()+', '+alg[i+1]+'":WeSn  -Y-5.5i -X-13i >> $NAME.eps'+'\n')
                n=n+1
            if ((i-n*3)%3 == 1):
                f5.write('awk \'{print $1, $2}\' '+props.getoutdir()+props.geteq()+'_mmi_costratio'+str(i)+'.txt |psxy -W3 -R -J -O -K -Ba5g5:"Cost Ratio, r":/a0.2g0.2:"Cost Savings Performance Metric, Q"::."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X6.5i >> $NAME.eps'+'\n')
            if ((i-n*3)%3 == 2):
                f5.write('awk \'{print $1, $2}\' '+props.getoutdir()+props.geteq()+'_mmi_costratio'+str(i)+'.txt |psxy -W3 -R -J -O -K -Ba5g5:"Cost Ratio, r":/a0.2g0.2:"Cost Savings Performance Metric, Q"::."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X6.5i >> $NAME.eps'+'\n')

        f5.write('exit 0'+'\n')
        f5.close()

        os.system('chmod 777' + ' ' + plotfile5)
        os.system('./'+plotfile5)

        #Map of warning Time
        outfile6 = props.getoutdir() + props.geteq() + '_warntime_map_' + props.getgmpe() + '_' + props.getgmice()+'_'+props.getmmicomp()+'_'+props.getrruporrp()
        plotfile6 = props.getoutdir()+'plot_warntime_map'
        print ('Creating figure', outfile6+'.eps')
        f6 = open(plotfile6,'w')
        f6.write('#!/bin/bash'+'\n')
        f6.write('export PATH='+props.getgmtpath()+':$PATH'+'\n')
        f6.write('  NAME='+outfile6+'\n')
        f6.write('rm $NAME.eps'+'\n')
        f6.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
        f6.write('gmt set PS_PAGE_ORIENTATION landscape'+'\n')
        f6.write('gmt set PS_MEDIA A2'+'\n')
        f6.write('gmt set MAP_TITLE_OFFSET 0.15c'+'\n')
        f6.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
        f6.write('gmt set FONT_TITLE 16p,Helvetica'+'\n')
        f6.write('makecpt -Chaxby -T0/40/5 > tmp.cpt'+'\n')
        f6.write('psbasemap -R'+str(eqlon-1)+'/'+str(eqlon+1)+'/'+str(eqlat-1)+'/'+str(eqlat+1)+' -JM8 -Ba1f0.25:."'+props.geteq()+', '+alg[1]+'":WeSn -Y11.5i -K > $NAME.eps'+'\n')
        f6.write('pscoast -R -J -Na -Df -O -W -K  >> $NAME.eps'+'\n')
        f6.write('awk \'{print $2, $3, -$6}\' '+props.getoutdir()+props.geteq()+'_mmi_quad0.txt |psxy -Sc0.2 -W -R -J -O -K -Ctmp.cpt >> $NAME.eps'+'\n')

        n=0
        for i in range(1,numalgs):
            if ((i-n*3)%3 == 0):
                f6.write('pscoast -R -J -Na -Df -O -W -K -Ba1f0.25:."'+props.geteq()+', '+alg[i+1]+'":WeSn  -Y-5.0i -X-9i >> $NAME.eps'+'\n')
                f6.write('awk \'{print $2, $3, -$6}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -Ctmp.cpt -R -J -O -K >> $NAME.eps'+'\n')               
                n=n+1
            if ((i-n*3)%3 == 1):
                f6.write('pscoast -R -J -Na -Df -O -W -K -Ba1f0.25:."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X4.5i >> $NAME.eps'+'\n')
                f6.write('awk \'{print $2, $3, -$6}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -Ctmp.cpt -R -J -O -K >> $NAME.eps'+'\n')

            if ((i-n*3)%3 == 2):
                f6.write('pscoast -R -J -Na -Df -O -W -K -Ba1f0.25:."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X4.5i >> $NAME.eps'+'\n')
                f6.write('awk \'{print $2, $3, -$6}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -Ctmp.cpt -R -J -O -K >> $NAME.eps'+'\n')

        f6.write('psscale -O -D4i/2i/4i/0.5i  -B10f5:"Warning Time (s)":/:: -Ctmp.cpt >> $NAME.eps'+'\n') 
        f6.write('exit 0'+'\n')
        f6.close()

        os.system('chmod 777' + ' ' + plotfile6)
        os.system('./'+plotfile6)

        #Map of MMI Bias
        outfile7 = props.getoutdir() + props.geteq() + '_mmibias_map_' + props.getgmpe() + '_' + props.getgmice()+'_'+props.getmmicomp()+'_'+props.getrruporrp()
        plotfile7 = props.getoutdir()+'plot_mmibias_map'
        print ('Creating figure', outfile7+'.eps')
        f7 = open(plotfile7,'w')
        f7.write('#!/bin/bash'+'\n')
        f7.write('export PATH='+props.getgmtpath()+':$PATH'+'\n')
        f7.write('  NAME='+outfile7+'\n')
        f7.write('rm $NAME.eps'+'\n')
        f7.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
        f7.write('gmt set PS_PAGE_ORIENTATION landscape'+'\n')
        f7.write('gmt set PS_MEDIA A2'+'\n')
        f7.write('gmt set MAP_TITLE_OFFSET 0.15c'+'\n')
        f7.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
        f7.write('gmt set FONT_TITLE 16p,Helvetica'+'\n')
        f7.write('makecpt -Cpolar -T-3/3/0.5 > tmp.cpt'+'\n')
        f7.write('psbasemap -R'+str(eqlon-1)+'/'+str(eqlon+1)+'/'+str(eqlat-1)+'/'+str(eqlat+1)+' -JM8 -Ba1f0.25:."'+props.geteq()+', '+alg[1]+'":WeSn -Y11.5i -K > $NAME.eps'+'\n')
        f7.write('pscoast -R -J -Na -Df -O -W -K  >> $NAME.eps'+'\n')
        f7.write('awk \'{print $2, $3, $4-$5}\' '+props.getoutdir()+props.geteq()+'_mmi_quad0.txt |psxy -Sc0.2 -W -R -J -O -K -Ctmp.cpt >> $NAME.eps'+'\n')

        n=0
        for i in range(1,numalgs):
            if ((i-n*3)%3 == 0):
                f7.write('pscoast -R -J -Na -Df -O -W -K -Ba1f0.25:."'+props.geteq()+', '+alg[i+1]+'":WeSn  -Y-5.0i -X-9i >> $NAME.eps'+'\n')
                f7.write('awk \'{print $2, $3, $4-$5}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -Ctmp.cpt -R -J -O -K >> $NAME.eps'+'\n')               
                n=n+1
            if ((i-n*3)%3 == 1):
                f7.write('pscoast -R -J -Na -Df -O -W -K -Ba1f0.25:."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X4.5i >> $NAME.eps'+'\n')
                f7.write('awk \'{print $2, $3, $4-$5}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -Ctmp.cpt -R -J -O -K >> $NAME.eps'+'\n')

            if ((i-n*3)%3 == 2):
                f7.write('pscoast -R -J -Na -Df -O -W -K -Ba1f0.25:."'+props.geteq()+', '+alg[i+1]+'":WeSn  -X4.5i >> $NAME.eps'+'\n')
                f7.write('awk \'{print $2, $3, $4-$5}\' '+props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt |psxy -Sc0.2 -W -Ctmp.cpt -R -J -O -K >> $NAME.eps'+'\n')

        f7.write('psscale -O -D4i/2i/4i/0.5i  -B1f1:"MMI Bias (Pred.-Obs.)":/:: -Ctmp.cpt >> $NAME.eps'+'\n') 
        f7.write('exit 0'+'\n')
        f7.close()

        os.system('chmod 777' + ' ' + plotfile7)
        os.system('./'+plotfile7)


        #MMI Warntime Density
        outfile8 = props.getoutdir() + props.geteq() + '_mmi_warntime_density_' + props.getgmpe() + '_' + props.getgmice()+'_'+props.getmmicomp()+'_'+props.getrruporrp()
        print ('Creating figure', outfile8+'.eps')
        plotfile8 = props.getoutdir()+'plot_mmiwarntimedensity'
        f8 = open(plotfile8,'w')
        f8.write('#!/bin/bash'+'\n')
        f8.write('export PATH='+props.getgmtpath()+':$PATH'+'\n')
        f8.write('  NAME='+outfile8+'\n')
        f8.write('rm $NAME.eps'+'\n')
        f8.write('gmt set FONT_ANNOT_PRIMARY 16p,Helvetica,black'+'\n')
        f8.write('gmt set PS_PAGE_ORIENTATION landscape'+'\n')
        f8.write('gmt set PS_MEDIA A2'+'\n')
        f8.write('gmt set MAP_TITLE_OFFSET 0.25c'+'\n')
        f8.write('gmt set FONT_LABEL 16p,Helvetica,black'+'\n')
        f8.write('gmt set FONT_TITLE 16p,Helvetica'+'\n')
        f8.write('makecpt -Chaxby -T0/0.5/0.01 -D -I > tmp.cpt'+'\n')
        f8.write('psbasemap -R'+props.getmmiwarnthreshold()+'/10/-10/50 -JX4i/4i -Ba1f0.5:"MMI":/a10f1:"Warning Time (s)"::."'+props.geteq()+', '+alg[1]+'":WeSn -Y12i  -K > $NAME.eps'+'\n')
        f8.write('psxy '+props.getoutdir()+props.geteq()+'_mmi_warntimedensity0.txt -Ctmp.cpt -R -J -P -L -O -K >> $NAME.eps'+'\n')
        f8.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
        f8.write('0 0'+'\n')
        f8.write('10 0'+'\n')
        f8.write('EOF'+'\n')
        n=0
        for i in range(1,numalgs):
            if ((i-n*3)%3 == 0):
                f8.write('psxy '+props.getoutdir()+props.geteq()+'_mmi_warntimedensity'+str(i)+'.txt -Ctmp.cpt -Ba1f0.5:"MMI":/a10f1:"Warning Time (s)"::."'+props.geteq()+', '+alg[i+1]+'":WeSn -R -J -P -L -O -K -Y-5.5i -X-11i  >> $NAME.eps'+'\n')
                f8.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f8.write('0 0'+'\n')
                f8.write('10 0'+'\n')
                f8.write('EOF'+'\n')
                n=n+1
            if ((i-n*3)%3 == 1):
                f8.write('psxy '+props.getoutdir()+props.geteq()+'_mmi_warntimedensity'+str(i)+'.txt -Ctmp.cpt -Ba1f0.5:"MMI":/a10f1:"Warning Time (s)"::."'+props.geteq()+', '+alg[i+1]+'":WeSn -R -J -P -L -O -K -X5.5i >> $NAME.eps'+'\n')
                f8.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f8.write('0 0'+'\n')
                f8.write('10 0'+'\n')
                f8.write('EOF'+'\n')
            if ((i-n*3)%3 == 2):
                f8.write('psxy '+props.getoutdir()+props.geteq()+'_mmi_warntimedensity'+str(i)+'.txt -Ctmp.cpt -Ba1f0.5:"MMI":/a10f1:"Warning Time (s)"::."'+props.geteq()+', '+alg[i+1]+'":WeSn -R -J -P -L -O -K -X5.5i >> $NAME.eps'+'\n')
                f8.write('psxy -W2 -R -J -O -K << EOF >> $NAME.eps'+'\n')
                f8.write('0 0'+'\n')
                f8.write('10 0'+'\n')
                f8.write('EOF'+'\n')

        if (numalgs%3 == 1):
            xoff='7'
            yoff='3.5'
        if (numalgs%3 == 2):
            xoff='0'
            yoff='-1'
        if (numalgs%3 == 0):
            xoff='-7'
            yoff='-1'
        f8.write('psscale -O -D3i/0i/6i/0.5ih  -B0.5f0.1:"Probability":/:: -Ctmp.cpt -X'+xoff+'i -Y'+yoff+'i >> $NAME.eps'+'\n') 
        f8.write('exit 0'+'\n')
        f8.close()

        os.system('chmod 777' + ' ' + plotfile8)
        os.system('./'+plotfile8)

###########################################################
#This section takes in the MMI evolution file and creates
#files more suitable for plotting
###########################################################

def analyzemmidata(evolfile, numalgs, props):
    file_length=0
    with open(evolfile, 'rt') as g:
        rows = (line.split() for line in g)
        for grow in rows:
            if (grow[0].startswith('#')):
                pass
            else:
                file_length = file_length+1
    site = list()
    lon = numpy.nan*numpy.ones([file_length,1])
    lat = numpy.nan*numpy.ones([file_length,1])
    t = numpy.nan*numpy.ones([file_length,1])
    mmiext = numpy.nan*numpy.ones([file_length,1])
    maxmmi = numpy.nan*numpy.ones([file_length,1])

    mmi = numpy.nan*numpy.ones([file_length,numalgs])
    k=0
    with open(evolfile, 'rt') as g:
        rows = (line.split() for line in g)
        for grow in rows:
            if (grow[0].startswith('#')):
                pass
            else:
                site.append(grow[1])
                lon[k,0]=float(grow[2])
                lat[k,0]=float(grow[3])
                t[k,0]=float(grow[4])
                mmiext[k,0]=float(grow[numalgs*3+8])
                maxmmi[k,0]=float(grow[numalgs*3+9])
                for i in range(0,numalgs):
                    mmi[k,i] = float(grow[9+3*i])
                    
                k=k+1
    wtime=t-mmiext #warning time
    
    numsites = len(numpy.unique(lon))
    aa = numpy.where(maxmmi >= float(props.getmmiwarnthreshold()))[0]
    totper = len(aa)/len(mmiext)*numsites
    
    mmirange = numpy.linspace(-5,5,num=51)

    for i in range(0,numalgs):
        evoloutfile=props.getoutdir()+props.geteq()+'_mmi_density'+str(i)+'.txt'
        cdfoutfile=props.getoutdir()+props.geteq()+'_mmi_cdf'+str(i)+'.txt'
        qtreefile=props.getoutdir()+props.geteq()+'_mmi_quad'+str(i)+'.txt'
        costratiofile=props.getoutdir()+props.geteq()+'_mmi_costratio'+str(i)+'.txt'
        warntimedensityfile=props.getoutdir()+props.geteq()+'_mmi_warntimedensity'+str(i)+'.txt'
        
        fevol = open(evoloutfile,'w')
        fcdf = open(cdfoutfile,'w')
        fqt = open(qtreefile,'w')
        fcr = open(costratiofile,'w')
        fwtd = open(warntimedensityfile,'w')

        mmialg = mmi[:,i]
        mmialg = mmialg.reshape((len(maxmmi),1))

        mmibias = mmialg-maxmmi
        #cdf and mmi evol files
        for j in range(0,100):
            secsum=0
            overwarned=0
            a3 = numpy.where((wtime == j-50) & (maxmmi < float(props.getmmiwarnthreshold())) & (mmialg >= float(props.getmmiwarnthreshold())))[0]
            tp = numpy.where((wtime == j-50) & (maxmmi >= float(props.getmmiwarnthreshold())) & (mmialg >= float(props.getmmiwarnthreshold())))[0]
            overwarned = overwarned+len(a3)/numsites
            secsum = secsum+len(tp)/totper
            for k in range(1, len(mmirange)):
                a1 = numpy.where((wtime == j-50) & (mmibias >= mmirange[k-1]) & (mmibias < mmirange[k]) & (maxmmi >= float(props.getmmiwarnthreshold())))[0]
                a2 = numpy.where((wtime == j-50) & (mmibias >= mmirange[k-1]) & (mmibias < mmirange[k]) & (maxmmi >= float(props.getmmiwarnthreshold())) & (mmialg > 4))[0]               
                mmihigh = "{0:.2f}".format(float(mmirange[k]))
                mmilow = "{0:.2f}".format(float(mmirange[k-1]))
                zval =  "{0:.4f}".format(float(len(a1)/totper))
                fevol.write('>-Z'+zval+'\n')
                fevol.write(str((j-50)*-1)+' '+mmilow+'\n')
                fevol.write(str((j-50)*-1)+' '+mmihigh+'\n')
                fevol.write(str((j-49)*-1)+' '+mmihigh+'\n')
                fevol.write(str((j-49)*-1)+' '+mmilow+'\n')
                fevol.write(str((j-50)*-1)+' '+mmilow+'\n')
                

                #secsum=secsum+len(a2)/totper
            ssum =  "{0:.4f}".format(float(secsum))
            owarn = "{0:.4f}".format(float(overwarned))
            fcdf.write(str((j-50)*-1)+' '+ssum+' '+owarn+'\n')
            
        fevol.close()
        fcdf.close()

        truepos=list()
        trueneg=list()
        falsepos=list()
        falseneg=list()

        mpredlist=list()
        wtlist=list()

        #quad plot of first alert
        uniquelons = numpy.unique(lon)
        for us in range(0,len(uniquelons)):
            sectoconsider = int(props.getttot())
            asite = numpy.arange(sectoconsider)+us*sectoconsider
            if (mmiext[asite[0]] < 200):
                siteid = site[asite[0]]
                mmisite = mmialg[asite]
                
                tsite = t[asite]
                twarn = tsite-mmiext[asite]
                aqt = numpy.where(numpy.absolute(mmisite) > 0)[0]

                mmi_pred = mmisite[aqt[0]] #mmi predicted by algorithm at site
                mmi_obs = maxmmi[asite[0]]
                tw = twarn[aqt[0]]

                if ((float(mmi_obs) >= float(props.getmmiwarnthreshold())) & (float(mmi_pred) >= float(props.getmmiwarnthreshold()))):
                    mpredlist.append(mmi_obs)
                    if (tw < -50):
                        tw = -50
                    wtlist.append(tw)

                mmip = "{0:.2f}".format(float(mmi_pred))
                mmio = "{0:.2f}".format(float(mmi_obs))
                tw =  "{0:.2f}".format(float(tw))
                lon_site = "{0:.4f}".format(float(lon[asite[0]]))
                lat_site = "{0:.4f}".format(float(lat[asite[0]]))

                fqt.write(siteid+' '+lon_site+' '+lat_site+' '+mmip+' '+mmio+' '+tw+'\n')

                if ((mmi_pred >= float(props.getmmiwarnthreshold())) & (mmi_obs >= float(props.getmmiwarnthreshold()))):
                    truepos.append(mmi_pred)
                if ((mmi_pred >= float(props.getmmiwarnthreshold())) & (mmi_obs < float(props.getmmiwarnthreshold()))):
                    falsepos.append(mmi_pred)
                if ((mmi_pred < float(props.getmmiwarnthreshold())) & (mmi_obs >= float(props.getmmiwarnthreshold()))):
                    falseneg.append(mmi_pred)
                if ((mmi_pred < float(props.getmmiwarnthreshold())) & (mmi_obs < float(props.getmmiwarnthreshold()))):
                    trueneg.append(mmi_pred)

        rrange = numpy.linspace(1.1,20,num=250)
        for rr in range(0,len(rrange)):
            TP = len(truepos)
            TN = len(trueneg)
            FP = len(falsepos)
            FN = len(falseneg)
            
            cr = (TP - (FP + (FP+FN))/(rrange[rr]-1))/(TP+FN)
            crout = "{0:.2f}".format(float(cr))
            rout = "{0:.2f}".format(float(rrange[rr]))
            fcr.write(rout+' '+crout+'\n')
                
        fqt.close()
        fcr.close()

        numrange = (10-float(props.getmmiwarnthreshold()))*2+1
        mmirange2 = numpy.linspace(float(props.getmmiwarnthreshold()),10,num=numrange)

        mpredarray = numpy.asarray(mpredlist)
        wtarray = numpy.asarray(wtlist)

        for j in range(0,60):
            for k in range (1, len(mmirange2)):
                a1 = numpy.where((wtarray >= j-50) & (wtarray < j-48) & (mpredarray >= mmirange2[k-1]) & (mpredarray < mmirange2[k]) )[0]

                a2 = numpy.where((mpredarray >= mmirange2[k-1]) & (mpredarray < mmirange2[k]) )[0]

                a3 = numpy.where((wtarray > 0) & (mpredarray >= mmirange2[k-1]) & (mpredarray < mmirange2[k]) )[0]

                totper2 = len(a2)
                if (totper2 == 0):
                    zval =  "{0:.4f}".format(float(0.0))
                    zval2 = "{0:.4f}".format(float(0.0))
                else:
                    zval =  "{0:.4f}".format(float(len(a1)/totper2))
                    zval2 =  "{0:.4f}".format(float(len(a3)/totper2))
                    

                mmihigh = "{0:.2f}".format(float(mmirange2[k]))
                mmilow = "{0:.2f}".format(float(mmirange2[k-1]))
                
                fwtd.write('>-Z'+zval+'\n')
                fwtd.write(mmilow+' '+str((j-50)*-1)+'\n')
                fwtd.write(mmihigh+' '+str((j-50)*-1)+'\n')
                fwtd.write(mmihigh+' '+str((j-48)*-1)+'\n')
                fwtd.write(mmilow+' '+str((j-48)*-1)+'\n')
                fwtd.write(mmilow+' '+str((j-50)*-1)+'\n')


        fwtd.close()
        
    

