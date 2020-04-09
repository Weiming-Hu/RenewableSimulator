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

import pandas as pd

from os import listdir, path
from progress.bar import IncrementalBar
from pvlib import pvsystem, irradiance, iotools


def read_hourly_SURFRAD(folder, progress=True):

    # List all data files in the folder
    folder = path.expanduser(folder)
    all_files = [path.join(folder, file) for file in listdir(folder)]

    # Initialize variables to hold data
    df_list = []
    meta = None

    if progress:
        print("Reading data from folder {}".format(folder))

    # Initialize a progress bar
    pbar = IncrementalBar("Reading files", max=len(all_files))
    pbar.suffix = '%(percent).1f%% - %(eta)ds'

    # Read each daily file as a data frame and append it to the data frame list
    for file in all_files:

        if progress:
            pbar.next()

        # Read file
        daily_surfrad = iotools.read_surfrad(file)

        # Parse data
        records = daily_surfrad[0]

        if meta is None:
            meta = daily_surfrad[1]

        # Subset only hourly rows
        records = records[records['minute'] == 0]

        # Append hourly data to the list
        df_list.append(records)

    if progress:
        pbar.finish()

    # Row bind all data frames to get the hourly data for the entire year
    hourly_data = pd.concat(df_list)

    return hourly_data, meta


def simulate_single_instance(ghi, dni_extra, tamb, wspd, albedo, current_time, surface_tilt, surface_azimuth,
                             pv_module, air_mass, tcell_model_parameters, solar_position):
    """
    Simulates PV energy production for a single instance at a single location for one timestamp.
    """

    # Decompose DNI from GHI
    #
    # TODO: There are different separation models.
    #
    dni_dict = irradiance.disc(ghi, solar_position["zenith"], current_time)

    # dni_dict = irradiance.erbs(ghi, solar_position["zenith"], current_time)

    # transmittance = forecast_model.cloud_cover_to_transmittance_linear(20)
    # dni_dict = irradiance.liujordan(solar_position["zenith"], transmittance, air_mass)

    dni = dni_dict["dni"]

    # Calculate POA sky diffuse
    #
    # TODO: There are different models to estimate diffuse radiation.
    #
    poa_sky_diffuse = irradiance.haydavies(
        surface_tilt, surface_azimuth, ghi, dni, dni_extra,
        solar_position["apparent_zenith"], solar_position["azimuth"])

    # Calculate POA ground diffuse
    #
    # TODO: There are different ground surface types.
    #
    poa_ground_diffuse = irradiance.get_ground_diffuse(surface_tilt, ghi, albedo)

    # Calculate angle of incidence
    aoi = irradiance.aoi(surface_tilt, surface_azimuth, solar_position["apparent_zenith"], solar_position["azimuth"])

    # Calculate POA total
    poa_irradiance = irradiance.poa_components(aoi, dni, poa_sky_diffuse, poa_ground_diffuse)

    # Calculate cell temperature
    tcell = pvsystem.temperature.sapm_cell(
        poa_irradiance['poa_global'], tamb, wspd,
        tcell_model_parameters['a'], tcell_model_parameters['b'], tcell_model_parameters["deltaT"])

    # Calculate effective irradiance
    effective_irradiance = pvsystem.sapm_effective_irradiance(
        poa_irradiance.poa_direct, poa_irradiance.poa_diffuse, air_mass, aoi, pv_module)

    # Calculate power
    sapm_out = pvsystem.sapm(effective_irradiance, tcell, pv_module)

    return sapm_out
