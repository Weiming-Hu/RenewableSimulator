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
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from pvlib import pvsystem, location, irradiance, atmosphere, temperature

# Self hosted modules
from Scenarios import Scenarios


def run_pv_simulations(nc_file, variable_dict, scenarios):

    # Sanity checks
    if not isinstance(scenarios, Scenarios):
        raise Exception("Input scenarios should an object of the class Scenarios")

    if not os.path.isfile(nc_file):
        raise Exception("File {} not exists".format(nc_file))

    # Open the NetCDF file
    nc = Dataset(nc_file, "a")

    # Read variables as a dictionary
    nc_vars = {key:nc.variables.get(value) for (key, value) in variable_dict.items()}

    # Sanity check to make sure that all variables are present and read successfully
    if not all(nc_vars.values()):
        raise Exception("Some variables are not found from the input file ({})".format(nc_file))

    # Determine the dimensions of the problem to simulate
    num_stations = nc.dimensions['num_stations'].size
    num_days = nc.dimensions['num_test_times'].size
    num_lead_times = nc.dimensions['num_flts'].size
    num_analogs = nc.dimensions['num_analogs'].size

    # Create our array to store results
    nc_output_group = nc.createGroup("PV_simulation")
    dc = nc_output_group.createVariable("DC", "f8", ("num_analogs", "num_flts", "num_test_times", "num_stations"))

    # TODO: What to do with these configurations
    surface_tilt = 0
    surface_azimuth = 180
    tcell_model_parameters = temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"]["open_rack_glass_glass"]
    pv_module = pvsystem.retrieve_sam("SandiaMod").Silevo_Triex_U300_Black__2014_

    # Batch run simulations
    for station_index in range(num_stations):

        # Determine the current location
        latitude = nc_vars["lat"][station_index].data
        longitude = nc_vars["lon"][station_index].data

        current_location = location.Location(latitude = latitude, longitude = longitude)

        for day_index in range(num_days):
            for lead_time_index in range(num_lead_times):

                # Determine the current time
                current_posix = nc_vars["date"][day_index] + nc_vars["flt"][lead_time_index]
                current_time = pd.Timestamp(current_posix, tz="UTC", unit='s')

                # Calculate sun position
                solar_position = current_location.get_solarposition(current_time)

                # Calculate extraterrestrial DNI
                dni_extra = irradiance.get_extra_radiation(current_time)

                # Calculate air mass
                air_mass = atmosphere.get_relative_airmass(solar_position["apparent_zenith"])

                for analog_index in range(num_analogs):

                    # Extra weather forecasts
                    ghi = nc_vars["ghi"][analog_index, lead_time_index, day_index, station_index].data
                    albedo = nc_vars["alb"][analog_index, lead_time_index, day_index, station_index] / 100
                    wspd = nc_vars["wspd"][analog_index, lead_time_index, day_index, station_index].data
                    tamb = nc_vars["tamb"][analog_index, lead_time_index, day_index, station_index] - 273.15

                    # Decompose DNI from GHI
                    dni_dict = irradiance.disc(ghi, solar_position["zenith"], current_time)
                    dni = dni_dict["dni"]

                    # Calculate POA sky diffuse
                    poa_sky_diffuse = irradiance.haydavies(
                        surface_tilt, surface_azimuth, ghi, dni, dni_extra,
                        solar_position["apparent_zenith"], solar_position["azimuth"])

                    # Calculate POA ground diffuse
                    poa_ground_diffuse = irradiance.get_ground_diffuse(surface_tilt, ghi, albedo)

                    # Calculate angle of incidence
                    aoi = irradiance.aoi(surface_tilt, surface_azimuth, solar_position["apparent_zenith"],
                                         solar_position["azimuth"])

                    # Calculate POA total
                    poa_irradiance = irradiance.poa_components(aoi, dni, poa_sky_diffuse, poa_ground_diffuse)

                    # Calculate cell temperature
                    tcell = pvsystem.temperature.sapm_cell(
                        poa_irradiance['poa_global'], tamb, wspd,
                        tcell_model_parameters['a'], tcell_model_parameters['b'], tcell_model_parameters["deltaT"])

                    # Calculate DC power
                    effective_irradiance = pvsystem.sapm_effective_irradiance(
                        poa_irradiance.poa_direct, poa_irradiance.poa_diffuse, air_mass, aoi, pv_module)

                    # Assign results
                    dc[analog_index, lead_time_index, day_index, station_index] = effective_irradiance



    nc.close()


if __name__ == '__main__':

    welcome_msg = "evergreen -- our pursuit of a more sustainable future\n" + \
        "Developed by Weiming Hu [weiming-hu.github.io]\n"

    print(welcome_msg)

    # This is my input NetCDF file generated from Analog Ensemble
    nc_file =  os.path.expanduser("~/data/analogs_PA.nc")

    # This is a dictionary of variables names in the PV simulator (keys) and the Analog Ensemble weather file (values)
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


    # Define multiple scenarios to simulate.
    scenarios = Scenarios();

    # TODO: what are some scenarios to try

    # Run the simulator
    run_pv_simulations(nc_file, variable_dict, scenarios)

    # Finishing
    print("PV simulation is complete!")
