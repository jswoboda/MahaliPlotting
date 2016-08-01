#!/bin/sh

# we use this instead of setup.py so that GeoDataPython is also installed in develop mode.

conda install --file requirements.txt
pip install pathlib2

conda install -c menpo mayavi

(
cd ..

git clone https://github.com/jswoboda/GeoDataPython
cd GeoDataPython
git pull #in case it was already installed
python setup.py develop
)
