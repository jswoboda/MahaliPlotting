#!/usr/bin/env python
"""  

"""
import os
from GeoData.filescraping import datedwebsite




def main():
    dirloc = 'https://amisr.asf.alaska.edu/PKR/DASC/RAW/'
    
    timedist = [1444608000,1444694400]
    basedir = os.path.expanduser('~/Documents/python/mahali/allsky/')
    datedwebsite(dirloc,timedist,basedir)
if __name__ == "__main__":
     main()