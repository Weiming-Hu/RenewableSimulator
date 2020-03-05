# Evergreen: A Renewable Energy Simulator

## Introduction

`evergreen` is implemented with [the System Advisory Model (SAM)](https://sam.nrel.gov/). SAM is a performance and financial model designed to facilitate decision making for people involved in the renewable energy industry.

This simulator focuses on solar photovoltaic power and wind power simulation. It uses ensemble weather predictions from [the Analog Ensemble package](https://weiming-hu.github.io/AnalogsEnsemble/).

## Installation

The following script walks you through the steps of installing `evergreen`.

```
# If you are on a cluster, load the dependencies.
# This is tested on XSEDE Stampede2
#
module load boost netcdf

# Download the latest version of source files
wget https://github.com/Weiming-Hu/RenewableSimulator/archive/master.zip

# Extract contents
unzip master.zip

# Carry an out-of-tree build
cd RenewableSimulator-master/
mkdir build && cd build

# Note that the CMAKE_PREFIX_PATH is where the AnEn and AnEnIO libraries are
CC=icc CXX=icpc cmake -DCMAKE_PREFIX_PATH=~/packages/release/lib -DCMAKE_INSTALL_PREFIX=~/packages/release -DCMAKE_INSTALL_RPATH="`echo ~`/packages/release/lib;$TACC_BOOST_LIB;$TACC_NETCDF_LIB" ..

# Build
make

# Install
make install
```

## Feedbacks

We appreciate collaborations and feedbacks from users. Please contact maintainer [Weiming Hu](http://weiming.ddns.net) through [weiming@psu.edu](weiming@psu.edu), or create tickets if you have any problems.

Thank you!

```
# "`-''-/").___..--''"`-._
#  (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
#  (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#    _ ..`--'_..-_/  /--'_.' ,'
#  (il),-''  (li),'  ((!.-'
# 
# Authors: 
#     Weiming Hu <weiming@psu.edu>
#     Guido Cervone <cervone@psu.edu>
#
# Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
# Department of Geography and Institute for CyberScience
# The Pennsylvania State University
```
