# MahaliPlotting
This will hold the sensor fussion plotting code for the Mahali project


# Requirements
This runs on Python 2.7.9. The packages required include
* numpy
* scipy
* pytables
* matplotlib
* basemap
* [GeoDataPython](https://github.com/jswoboda/GeoDataPython)

It is suggested that the user install the anaconda distribution of python. This software has be tested on Linux and Mac machines only.

# Use

Currently the plotting is done through the python file plotdata.py. This fill will use the GeoDataPython module to interpolate and plot images of all sky data and plot vertical TEC in a scatter plot at pierce points. The TEC data is read in from ionofiles. The all sky data is from [Don Hampton](https://amisr.asf.alaska.edu/PKR/DASC/RAW/) at UAF and saved as FITS files. 

The user can run this code from the command line but first they must download the data into directries that can be accessed. The user can then run the code as so

	$python plotdata.py -a <all skypath> -w <wavelength>, -i <ionofile dir>, -t <time interval>,  -p <plotdirectory> -r <type y to reinterpolate all sky data> 
	
Or as in this example

	$python plotdata.py -a ~/DATA/Mahali/allsky -i ~/DATA/Mahali/GPS_Processed_IonoFiles/site_280_2015/ -w 558 -t 2 -p ~/Documents/python/MahaliPlotting/plots10172015
	
Future updates will allow the user to download all sky data and give a method to select times using the command line.