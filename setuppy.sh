#!/bin/sh

wget -nc http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda2
export PATH=$HOME/miniconda2/bin:$PATH

