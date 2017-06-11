#!/bin/bash

# we use this instead of setup.py so that GeoDataPython is also installed in develop mode.

(
cd ..

git clone https://github.com/jswoboda/GeoDataPython
cd GeoDataPython
git pull #in case it was already installed
python setup.py develop
)

python setup.py develop
