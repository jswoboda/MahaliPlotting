#!/bin/sh

#----- setup environ

py2path=$(dirname $(which python2))

echo $py2path

${py2path}/conda install --file requirements.txt
(
cd ..

git clone https://github.com/jswoboda/GeoDataPython
cd GeoDataPython
python2 setup.py develop
)
