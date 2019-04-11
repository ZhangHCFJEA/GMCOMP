# GMCOMP
Comparisons of Ground Motion Predictions for EEW

This code package is written to compare ground motion predictions made from different algorithms. The original intent is to compare the source models coming out of geodetic algorithms for ShakeAlert, the different seismic algorithms, and the ground truth (real PGA/PGV values recorded at seismic stations). 


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



