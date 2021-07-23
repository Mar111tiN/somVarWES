#!/bin/sh


DATA="/fast/users/szyskam_c/scratch/snaketest/data"

# github repos
EB="https://github.com/Mar111tiN/ebscore.git"
HDR="https://github.com/Mar111tiN/HDRdetect.git"
CNV="https://github.com/Mar111tiN/myCNV.git"
PRM="https://github.com/Mar111tiN/primertools.git"

# load the packages from github
cd scripts && \
git clone $EB && \
git clone $HDR && \
git clone $CNV && \
git clone $PRM && \
cd ..

# create the snake-env environment (if it does not exist yet)
conda env create -f env/snake-env.yml -n snake-env

 