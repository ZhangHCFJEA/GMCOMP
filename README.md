# GMCOMP
<b>Comparisons of Ground Motion Predictions for EEW</b>

This code package is written to compare ground motion predictions made from different algorithms. The original intent is to compare the source models coming out of geodetic algorithms for ShakeAlert, the different seismic algorithms, and the ground truth (real PGA/PGV values recorded at seismic stations). 

The usage of this code package is fully dictated by the .props file. If the name of the properties file is changed, only the following line towards the top of <i>GMCOMP.py</i> needs to be modified:

<i>props = Properties('gmcomp.props')</i>

As such, maintaining the integrity of the .props file is key. The variables should have no spaces on a line such that:

<i>variable=value</i>

Since there are lots of variables and paths, my suggestion is to copy your original .props file and create a new one for each earthquake run you want to perform such that:

<i>earthquake1.props</i>, <i>earthquake2.props</i>, etc.

Output file names will be modified according to specific properties you use such as earthquake name, gmpe used, gmice used, definition of peak motion, and algorithm name, so you don't need to make a new properties file for small changes such as gmpe used since the output files will not overwrite others for the same earthquake.

<hr style="height:1px;border:none;color:#333;background-color:#333;" />

<b>Prerequisites</b>

In order to run GMCOMP, you will need several Python3 packages. On OSX, use <i>pip3</i> for all of these, on other operating systems it will be dependent on your package manager, but these are fairly standard packages.

<i>obspy</i>

<i>shapely</i>

<i>numpy</i>

For plotting, I have set up several GMT scripts. The most annoying thing with GMT is ensuring your path is correct. I have added the GMT path as a variable in the props file; this should be the same as in your .bashrc file. If your version of GMT is finnicky, simply turn off plotting in the props file.


<b>Ground Motion Models (GMMs):</b>

There are two types of GMMs within GMCOMP, Ground Motion Prediction Equations (GMPEs) and Ground Motion Intensity Conversion Equations (GMICEs). The two Python packages are aptly called <i>GMPE.py</i> and <i>GMICE.py</i>.

<i>GMPE.py</i> - There are currently 4 GMPEs within the main package that compute PGA or PGV. The name of the GMPE is stored in the properties dictionary (read from the .props file), and is called within GMPE.py through <i>props.getgmpe()</i>. To perform ground motion predictions, the following is called:

<i>Y = GMPE.Y(M, Rrup, Rjb, Rx, Ztor, Zhyp, VS30, Hw, W, dip, rake, mode, props)</i>

and the return is either the values of PGA or PGV at the locations specified. Here, the locations of interest are embedded within the different distance measures and site conditions.

The variables are:

<i>M</i> - magnitude (single value)

<i>Rrup</i> - Rupture distance (km, a numpy array, n by 1, where n is the number of sites/locations you are interested in)

<i>Rjb</i> - Joyner Boore distance (km, array)

<i>Rx</i> - Distance to the updip projection of the fault plane (km, array)

<i>Ztor</i> - Depth to the top of the rupture plane (km, array)

<i>Zhyp </i>- hypocentral depth (km, single value)

<i>VS30</i> - Shear wave velocity in upper 30 meters at site (m/s, array)

<i>Hw</i> - Hanging wall flag (1 for on hanging wall or 0 for not, array)

<i>W</i> - width of fault plane (km, single value)

<i>dip</i> - dip of the fault plane (degrees, single value)

<i>rake</i> - rake of the earthwauke (degrees, single value)

<i>mode</i> - 0 for PGA, 1 for PGV. PGV is not fully operational

<i>props</i> - the dictionary of properties.

The 4 different GMPEs are:

<i>cy08</i> - Chiou and Youngs (2008)

<i>cb14</i> - Campbell and Bozorgnia (2014)

<i>ba08</i> - Boore and Atkinson (2008)

<i>bradley13</i> - Bradley (2013)


<i> GMICE.py </i> - There are two GMICEs embedded within GMICE.py. They both can take in either the PGA or PGV values and convert to instrumental MMI. The particular GMICE equation used is embedded in the properties dictionary, accessible through props.getgmice(). The usage is the following:

<i>MMI = GMICE.MMI_Y(Y, mode, props)</i>

where <i>Y</i> is either the PGA or PGV values, <i>mode</i> is 0 for PGA, 1 for PGV, and <i>props</i> is the dictionary of properties. The two GMICE equations are as follows:

<i>wald99</i> - Wald et al. (1999)

<i>wgrw12</i> - Worden et al. (2012)


<b>GMCOMP_coord_tools.py</b>

The coord_tools package performs standard coordinate transformations. All subroutines are array capable. Not all of them are used, but the following subroutines are available:

GMCOMP_coord_tools.lla2ecef(lat,lon,alt) - lon,lat,alt to ITRF cartesian

GMCOMP_coord_tools.ecef2lla(x,y,z) - ITRF Cartesian to lon,lat,alt

GMCOMP_coord_tools.dxyz2dneu(dx,dy,dz,lat,lon) - displacement in xyz to neu

GMCOMP_coord_tools.ll2utm(lon,lat,lon0,lat0) - lon/lat to UTM

GMCOMP_coord_tools.utm2ll(UTMeasting,UTMnorthing,lon0,lat0) - UTM to lon/lat


<b>Plotting</b>

In the .props file, there are three variables related to plotting. <b><i>plotson</i></b> is either <i>yes</i> or <i>no</i> to turn on plotting, <b><i>gmtpath</b></i> controls the path to your GMT bin folder, and <b><i>gmtversion</b></i> gives the overarching version of GMT (right now only 5 is supported but I will write GMT4 scripts eventually).

On Macs, the easiest way to get GMT5 is through homebrew such that:

<i>brew install gmt@5</i>

This will place your GMT path in: 

<i>/usr/local/Cellar/gmt/"version"/bin</i>

The main plotting script is <i>GMCOMP_gmtplots.py</i> and this will be called at the end of <i>GMCOMP.py</i>. A number of GMT scripts will be created in your output directory (outputdir in .props), and you can modify these to change formatting of figures according to your needs. You can simply go to your output directory, modify, and rerun the following from the command line in the main GMCOMP directory (since there are path references that wouldn't work within the output directory):

<hr>

./plot_mmibias

This creates MMI Bias PDFs as a function of warning time. If the warning threshold is not exceeded, I use theoretical S wave traveltimes derived from the Obspy taup package. 

![MMI Bias Density Figure](../master/plots/Napa2014_mmibias_density_cy08_wald99_m_rp.jpg)

<hr>

./plot_cdf_fp

This creates a pseudo-cdf and false-positive chart. What I mean by pseudo-cdf is the plot is the cumulative number of warned stations as a function of warning time, however, that number can decrease as the model results are updated. This would assume optimal cancellations of alerts, which is not the current paradigm. 

![CDF False Positive Figure](../master/plots/Napa2014_cdf_fp_cy08_wald99_m_rp.jpg)

<hr>

./plot_quad

This creates a plot of observed versus predicted MMI. The reason I call it a quad plot is the four quadrants can be segregated into true positive (TP), true negative (TN), false positive (FP), and false negative (FN). The dots are color-coded by warning time, fully saturated black dots represent no or negative warning time. The vertical and horizontal lines are modified by the mmiwarnthreshold defined in the .props file.

![Quad Figure](../master/plots/Napa2014_quad_cy08_wald99_m_rp.jpg)

<hr>

./plot_warntime

This plots the observed MMI versus the warning time. The color coding is by MMI bias (predicted minus observed). The horizontal line is modified by the mmiwarnthreshold defined in the .props file.

![Warning Time MMI Figure](../master/plots/Napa2014_warntime_obsmmi_cy08_wald99_m_rp.jpg)

<hr>

./plot_costratio

This plots the cost savings performance metric, Q as a function of the cost ratio, r.

![Cost Ratio Figure](../master/plots/Napa2014_costratio_cy08_wald99_m_rp.jpg)

<hr>

./plot_mmibias_map

This plots the MMI bias of the first alert for an algorithm.

![MMI Bias Map](../master/plots/Napa2014_mmibias_map_cy08_wald99_m_rp.jpg)

<hr>

./plot_warntime_map

This plots the warning time map for a given algorithm based on the first alert times.

![Warning Time Map](../master/plots/Napa2014_warntime_map_cy08_wald99_m_rp.jpg)

<hr>

./plot_mmiwarntimedensity

This plots the probability of warning times based upon a given mmi band. I have set this to 0.5 MMI units between the warning threshold and MMI=10.

![Warning Time Density Plot](../master/plots/Napa2014_mmi_warntime_density_cy08_wald99_m_rp.jpg)





