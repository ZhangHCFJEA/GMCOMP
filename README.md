# GMCOMP
Comparisons of Ground Motion Predictions for EEW

This code package is written to compare ground motion predictions made from different algorithms. The original intent is to compare the source models coming out of geodetic algorithms for ShakeAlert, the different seismic algorithms, and the ground truth (real PGA/PGV values recorded at seismic stations). 


<b>Ground Motion Models (GMMs):</b>

There are two types of GMMs within GMCOMP, Ground Motion Prediction Equations (GMPEs) and Ground Motion Intensity Conversion Equations (GMICEs). The two Python packages are aptly called <i>GMPE.py</i> and <i>GMICE.py</i>.

<i>GMPE.py</i> - There are currently 4 GMPEs within the main package that compute PGA or PGV. The name of the GMPE is stored in the properties dictionary (read from the .props file), and is called within GMPE.py through props.getgmpe(). To perform ground motion predictions, the following is called:

GMPE.Y(M, Rrup, Rjb, Rx, Ztor, Zhyp, VS30, Hw, W, dip, rake, mode, props)

and the return is either the values of PGA or PGV at the locations specified. Here, the locations of interest are embedded within the different distance measures and site conditions.

The variables are:
M - magnitude (single value)
Rrup - Rupture distance (km, a numpy array, n by 1, where n is the number of sites/locations you are interested in)
Rjb - Joyner Boore distance (km, array)
Rx - Distance to the updip projection of the fault plane (km, array)
Ztor - Depth to the top of the rupture plane (km, array)
Zhyp - hypocentral depth (km, single value)
VS30 - Shear wave velocity in upper 30 meters at site (m/s, array)
Hw - Hanging wall flag (1 for on hanging wall or 0 for not, array)
W - width of fault plane (km, single value)
dip - dip of the fault plane (degrees, single value)
rake - rake of the earthwauke (degrees, single value)
mode - 0 for PGA, 1 for PGV. PGV is not fully operational
props - the dictionary of properties.

The 4 different GMPEs are:

cy08 - Chiou and Youngs (2008)
cb14 - Campbell and Bozorgnia (2014)
ba08 - Boore and Atkinson (2008)
bradley13 - Bradley (2013)

