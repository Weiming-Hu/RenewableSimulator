# Evergreen: A Renewable Energy Simulator

## Introduction

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

This project is designed for assessing prediction uncertainty for renewable energy production. It builds on top of [Analog Ensemble](https://weiming-hu.github.io/AnalogsEnsemble/) and [pvlib](https://pvlib-python.readthedocs.io/en/stable/). It currently supports running with data from Analog Ensemble and from SURFRAD.

It is implemented in Python 3.

## Usage

For your reference, [here](https://github.com/Weiming-Hu/RenewableSimulator/issues/2) provides a collection of scripts to prepare the environment on HPC platforms.

Once you have the environment set up, you can use the following code to see the available options:

```shell script
python runner_pv_anen.py -h

python runner_pv_surfrad.py -h
```

## Profiling

### yappi

```shell script
# Run the program through the profiler
python runner_pv.py --profile --profiler yappi

# On Mac OS. kcachegrind on Linux
qcachegrind yappi_2020-04-12-18-08-48_rank-0.log
```

### pyinstrument

```shell script
# Run the program through the profiler
python runner_pv.py --profile --profiler pyinstrument
```

### line_profiler

```shell script
# Run the program through the profiler
kernprof -l runner_pv.py --profile --profiler line_profiler

# Generate text output
python -m line_profiler evergreen.py.lprof
```

### Simple Clock

If you are ready for production but still want to have a general idea of how much time was spent, the simple clock would be a good solution. It just uses `time.time()` to get current time so its overhead is pretty low.

```shell script
python runner_pv.py --profile
```

## Feedback

We appreciate collaborations and feedback from users. Please contact maintainer [Weiming Hu](http://weiming.ddns.net) or create tickets if you have any problems.

Thank you!

