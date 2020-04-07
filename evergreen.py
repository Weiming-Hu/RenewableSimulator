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
import os

# Scientific python add-ons
import pandas as pd
from netCDF4 import Dataset
from pvlib import pvsystem, location, irradiance, atmosphere

# Self hosted modules
from Scenarios import Scenarios

# Temporary modules that is used during development only


if __name__ == '__main__':

    welcome_msg = "evergreen -- our pursuit of a more sustainable future\n" + \
        "Developed by Weiming Hu [weiming-hu.github.io]\n"

    print(welcome_msg)

    # These are my input variables
    anen_file =  os.path.expanduser("~/data/analogs_PA.nc")

    variable_dict = {
        "ghi": "DownwardShortwaveRadiation",
        "alb": "Albedo",
        "tamb": "Temperature_2m",
        "wspd": "WindSpeed_10m",
        "lon": "Xs",
        "lat": "Ys",
        "date": "test_times",
        "flt": "FLTs",
    }


    if not os.path.isfile(anen_file):
        raise Exception("File {} not exists".format(anen_file))





    nc = Dataset(anen_file, "r")

    # Read variables
    nc_vars = {key:nc.variables.get(value) for (key, value) in variable_dict.items()}


    station_index = 100
    test_day_index = 100
    lead_time_index = 16
    member_index = 0

    latitude = nc_vars["lat"][station_index].data
    longitude = nc_vars["lon"][station_index] - 360

    surface_tilt = 0
    surface_azimuth = 180
    albedo = nc_vars["alb"][member_index, lead_time_index, test_day_index, station_index] / 100


    ghi = nc_vars["ghi"][member_index, lead_time_index, test_day_index, station_index].data


    # Solar position
    posix_timestampe = nc_vars["date"][test_day_index] + nc_vars["flt"][lead_time_index]
    current_time = pd.Timestamp(posix_timestampe, tz = "UTC", unit = 's')
    current_location = location.Location(latitude = latitude, longitude = nc_vars["lon"][station_index])

    solar_position = current_location.get_solarposition(current_time)

    # Extraterrestrial DNI
    dni_extra = irradiance.get_extra_radiation(current_time)

    # Air mass
    air_mass = atmosphere.get_relative_airmass(solar_position["apparent_zenith"])

    # POA sky diffuse
    dni_dict = irradiance.disc(ghi, solar_position["zenith"], current_time)
    poa_sky_diffuse = irradiance.haydavies(surface_tilt, surface_azimuth, ghi, dni_dict["dni"], dni_extra, solar_position["apparent_zenith"], solar_position["azimuth"])

    # POA ground diffuse
    poa_ground_diffuse = irradiance.get_ground_diffuse(surface_tilt, ghi, albedo)

    # AOI
    aoi = irradiance.aoi(surface_tilt, surface_azimuth, solar_position["apparent_zenith"], solar_position["azimuth"])

    # POA total
    poa_irradiance = irradiance.poa_components(aoi, dni_dict["dni"], poa_sky_diffuse, poa_ground_diffuse)

    # Cell and module temperature
    wspd = nc_vars["wspd"][member_index, lead_time_index, test_day_index, station_index].data
    ambient_temperature = nc_vars["tamb"][member_index, lead_time_index, test_day_index, station_index] - 273.15
    cell_temperature = pvsystem.temperature.sapm_cell(poa_irradiance['poa_global'], ambient_temperature, wspd, -3.47, -0.0594, 3)

    # Select module
    module = pvsystem.retrieve_sam("SandiaMod").Silevo_Triex_U300_Black__2014_

    # DC power
    effective_irradiance = pvsystem.sapm_effective_irradiance(poa_irradiance.poa_direct, poa_irradiance.poa_diffuse, air_mass, aoi, module)

    nc.close()

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




