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

# External packages
import os
from netCDF4 import Dataset
from PySAM import Pvwattsv51ts as pv

# Internal packages
from Scenarios import Scenarios


if __name__ == '__main__':

    # welcome_msg = "evergreen -- our pursuit of a more sustainable future\n" + \
    #     "Developed by Weiming Hu [weiming-hu.github.io]\n"
    #
    # print(welcome_msg)
    #
    # # These are my input variables
    # anen_file =  os.path.expanduser("~/data/analogs_PA.nc")
    #
    # if not os.path.isfile(anen_file):
    #     raise Exception("File {} not exists".format(anen_file))
    #
    #
    # """
    # Set up scenarios to simulate
    # """
    #
    # # Define multiple scenarios to simulate.
    # # Scenarios will be created using permutation of values from all keys.
    # #
    # pvwatts_scenarios = Scenarios();
    # pvwatts_scenarios["system_capacity"] = [1, 2, 3, 4];    # System size (DC nameplate) in kW
    # pvwatts_scenarios["module_type"] = [0, 1, 2];           # 0 for standard, 1 for premium, and 2 for thin film
    # pvwatts_scenarios["losses"] = [0, 20, 40, 80, 99];      # Percent of total system loss within [-5, 99]
    #
    # # Define fixed configuration
    # pvwatts_config = {
    #     "time_step": 1,                                     # SSC default
    #     "dc_ac_ratio": 1.1,                                 # SSC default
    #     "inv_eff": 96,                                      # SSC default
    #     "array_type": 0,                                    # Fixed panel
    #     "gcr": 0.4,                                         # SSC default ground coverage ratio
    #     "tilt": 0,                                          # Horizontal
    # }
    #
    # pvwatts_scenarios.print()
    #
    #
    # """
    # Read input from an NetCDF file
    # """
    # rootgrp = Dataset(anen_file, "r")
    #
    # rootgrp.close()
    #
    #
    # """
    # Run batch simulations
    # """
    # index = 1
    #
    # current_scenario = pvwatts_scenarios.get_scenario(index)

    mydict = {
        "PVWatts": {
            "beam": 100,
            "diffuse": 20,
            "poa": 100,
            "alb": 0.13,
            "wspd": 1.14,
            "tamb": 4.75,
            "tcell": 4.75,

            # "time_step": 1,

            "year": 2018,
            "month": 12,
            "day": 14,
            "hour": 8,
            "minute": 0,

            "lat": 39,
            "lon": -80,
            "tz": -5,
        },

        "SystemDesign": {
            "array_type": 0,
            "azimuth": 180,
            # "dc_ac_ratio": 0.9,
            # "gcr": 0.4,
            # "inv_eff": 99,
            "losses": 0,
            # "module_type": 0,
            "system_capacity": 10,
            "tilt": 0,
        }
    }

    pvwatts_module = pv.new()

    pvwatts_module.assign(mydict)

    pvwatts_module.execute(1)
