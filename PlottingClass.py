#!/usr/bin/env python

from MahaliPlotting import runPlotClass2

if __name__== '__main__':

    from argparse import ArgumentParser
    descr = '''
             This script will run Mahali Plotting software. The user needs to
             specify the locations of the different types of data and the time
             limits if they don't want all of the data processed. This specific
             module runs everything in a class struture so it can be used with
             other modules easier.
            '''
    p = ArgumentParser(description=descr)
    p.add_argument("ini",help='Config file name.')
    p.add_argument("radar",help='Radar hdf5 file')
    p.add_argument("asky",help='The allsky data directory.')
    p.add_argument("-i","--ifile",help='The directory that holds all of the TEC data in ionofile formats.')
    p.add_argument("-p", "--pdir",help='plot output directory',default='.')
    p = p.parse_args()

    runPlotClass2(p.ini, p.pdir, p.ifile, p.asky, p.radar)
