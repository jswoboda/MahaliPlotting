==============
MahaliPlotting
==============
This holds sensor fusion plotting code for the Mahali project. 
Python 3 or 2.7.

.. contents::

Requirements
============
One time Python environment setup, IF you don't already have Python::

    ./reference/setuppy.sh

Setup of program (can be issued for existing Python install)::

    ./setup.sh

It is suggested that the user install the Miniconda/Anaconda distribution of python. 

Get DASC raw data
=================
https://github.com/scivision/dascutils

Use
===
Currently the plotting is done through the python file plotdata.py. This fill will use the GeoDataPython module to interpolate and plot images of all sky data and plot vertical TEC in a scatter plot at pierce points. The TEC data is read in from ionofiles. The all sky data is from [Don Hampton](https://amisr.asf.alaska.edu/PKR/DASC/RAW/) at UAF and saved as FITS files. 

The user can run this code from the command line but first they must download the data into directries that can be accessed. The user can then run the code as::

	python plotdata.py -a <all skypath> -w <wavelength>, -i <ionofile dir>, -t <time interval>,  -p <plotdirectory> -r <type y to reinterpolate all sky data> 
	
Or as in this example::

	python plotdata.py -a ~/DATA/Mahali/allsky -i ~/DATA/Mahali/GPS_Processed_IonoFiles/site_280_2015/ -w 558 -t 2 -p ~/Documents/python/MahaliPlotting/plots10172015
	
Future updates will allow the user to download all sky data and give a method to select times using the command line.
