# "`-''-/").___..--''"`-._
#  (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
#  (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#    _ ..`--'_..-_/  /--'_.' ,'
#  (il),-''  (li),'  ((!.-'
#
# Author:
#
#     Weiming Hu <weiming@psu.edu>
#
# Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
# Department of Geography and Institute for CyberScience
# The Pennsylvania State University
#

# Built-in python modules
import argparse

import pandas as pd
from pvlib import location, atmosphere, temperature, pvsystem, irradiance


if __name__ == '__main__':

    welcome_msg = "Test get_solarposition"
    parser = argparse.ArgumentParser(description=welcome_msg)
    parser.add_argument('--iterations', type=int, default='1000')
    parser.add_argument('--method', default='nrel_numpy')

    args = parser.parse_args()

    for it in range(args.iterations):
        current_location = location.Location(latitude=40, longitude=-70)
        current_posix = 1514782800
        current_time = pd.Timestamp(current_posix, tz="UTC", unit='s')
        solar_position = current_location.get_solarposition(current_time, method=args.method)
