#!/bin/sh

#----- setup environ

py2path=$(dirname $(which python2))

echo $py2path

${py2path}/conda install --file requirements.txt

${py2path}/pip install https://github.com/jswoboda/GeoDataPython/tarball/master#egg=GeoDataPython
