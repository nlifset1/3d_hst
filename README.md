# 3d_hst

**There are three groups of scripts: Aperture Radius, Nth Nearest, and Colors.**

##Nth Nearest

###Tabler
This script takes in the 3dhst all fields data file and creates a new file called *values.dat* with info relevant to the nth nearest number density calculation.
###Mass_Density_Z
This script takes in the 3dhst all fields file as well as the *values.dat* file and then creates a plot which shows galaxy number density per galaxy mass in four redshift bins.

##Aperture Radius

###Tabler_Radius
This script takes in the 3dhst all fields file and creates a new file called *values_R.dat* which has data on the galaxy number density and random background density per radius for central galaxies.
###Tabler_Radius2
This script uses the 3dhst all fields catalog as well as *values_R.dat* and creates a new file called *values_R2* which has the same data as *values_R.dat* but slit into bins based on satellite mass.
###Z_error
This script uses the 3dhst all fields catalog as well as *values_R.dat* and creates a new file called *z_errors.dat* which has error bar data for Radius_Density _Z
###Radius_Density _Z
This script uses the 3dhst all fields catalog, *values_r.dat*, and *z_errors.dat* and creates a plot the shows galaxy number density per aperture radius binned by central galaxy redshift.
###Lmass_error
This script uses the 3dhst all fields catalog as well as *values_R.dat* and creates a new file called *lmass_errors.dat* which has error bar data for Radius_Density _Mass
###Radius_Density _Mass
This script uses the 3dhst all fields catalog, *values_r.dat*, and *lmass_errors.dat* and creates a plot the shows galaxy number density per aperture radius binned by central galaxy mass.
###satmass_error
This script uses the 3dhst all fields catalog as well as *values_R2.dat* and creates a new file called *lmass_errors2.dat* which has error bar data for Radius_Density _Satmass.
###Radius_Density _Satmass
This script uses  *values_r2.dat* and *lmass_errors2.dat*, and creates a plot the shows galaxy number density per aperture radius binned by satellite galaxy mass.
###Radius_Density _ZTal
This script uses the 3dhst all fields catalog, *values_r.dat*, and *z_errors.dat* and creates a plot the shows galaxy number density per aperture radius binned by central galaxy redshift and includes Tomer Tal's 2013 data.
###Radius_Density _Compare
This script uses the 3dhst all fields catalog and *values_r.dat* and creates a plot that shows galaxy number density per aperture radius with and without the background subtraction.

##Colors

###concatenator
This script takes in all of the individual 3dhst field files and concatenates them into one, because the all fields file has no U,V, or J data. This creates the new file *3dhst_big.dat*.
###Color_Table _initial
This script uses *3dhst_big.dat* and creates a new file called *Color_values.dat* which has SFR color data for each galaxy.
###Color_Table
This script uses *color_values.dat* and creates a new file *values_color.dat* which has density data binned by sfr color.
###Color_Table _error
This script uses *values_color.dat* and creates a new file called *color_errors.dat* which has error bar data for Color_Panel.
###Color_Panel
This script uses *values_color.dat* and creates a plot of galaxy number density per aperture radius binned by central galaxy sfr color, satellite galaxy sfr color, and central galaxy redshift.
###Color_Table _conformity
This script uses *values_color.dat* and creates a new file called *color_tabled.dat* which has data on satellite quiescent percentage for each central galaxy binned by central galaxy sfr color.
###Galaxy_conformity1
This script uses *color_tabled.dat* and creates a plot that shows satellite quiescent percentage per aperture radius with centrals binned by sfr color and redshift.
###Galaxy_conformity2
This script uses *color_tabled.dat* and creates a plot that shows galaxy conformity per aperture radius with centrals binned by redshift.
###Color
This script uses the 3dhst all fields catalog as well as *3dhst_big.dat* and creates a U-V vs V-J plot with galaxies binned by mass.
###Color_Z
This script uses the 3dhst all fields catalog as well as *3dhst_big.dat* and creates a U-V vs V-J plot with galaxies binned by redshift.
