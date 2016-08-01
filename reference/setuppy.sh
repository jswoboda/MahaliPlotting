#!/bin/bash

wget -nc http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda3
export PATH=$HOME/miniconda3/bin:$PATH
echo 'export PATH=$HOME/miniconda3/bin:$PATH' >> ~/.bashrc
