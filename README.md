# A Weather-Coupled Renewable Energy Simulator

## Introduction

The renewable simulator is implemented with [the System Advisory Model (SAM)](https://sam.nrel.gov/). SAM is a performance and financial model designed to facilitate decision making for people involved in the renewable energy industry.

This simulator focuses on solar photovoltaic power and wind power simulation. It uses weather predictions from [the Analog Ensemble package](https://github.com/Weiming-Hu/AnalogsEnsemble) and prepare the data so that they can be fed into SAM functions. It also provides functions to view the different model available in SAM and their requirements.

## Installation

Please clone or download the repository. The program has the following requirements:

- Python 2.7
- Python package `netCDF4`

After checking the requirements, run the following check commands in a terminal:

```
$ python runme.py 
[pvwattsv5] Annual ennergy: 6791.42236328
[pvwattsv5_1ts] AC power: 39.2786483765; DC power: 67.7225570679
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
