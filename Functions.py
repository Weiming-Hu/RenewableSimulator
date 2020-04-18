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

import math
import numpy as np
import pandas as pd

from os import listdir, path
from progress.bar import IncrementalBar
from pvlib import pvsystem, irradiance, iotools, location, atmosphere


def get_start_index(total, num_procs, rank):
    """
    Get the start index of a given rank

    :param total: The total number of instances
    :param num_procs: The total number of processes
    :param rank: The current rank number
    :return: The start index of the instance for the current rank
    """

    if total < num_procs:
        raise Exception("You are probably wasting computing resources. You requested too many processes.")

    if rank >= num_procs:
        raise Exception("The rank index {} can not be greater or equal to the number of processes {}".format(
            rank, num_procs))

    return math.ceil(rank * total / num_procs)


def get_end_index(total, num_procs, rank):
    """
    Get the end index of a given rank

    :param total: The total number of instances
    :param num_procs: The total number of processes
    :param rank: The current rank number
    :return: The end index of the instance for the current rank. Note that this index is exclusive.
    """

    if total < num_procs:
        raise Exception("You are probably wasting computing resources. You requested too many processes.")

    if rank >= num_procs:
        raise Exception("The rank index {} can not be greater or equal to the number of processes {}".format(
            rank, num_procs))

    if rank == num_procs:
        return total

    return math.ceil((rank + 1) * total / num_procs)


def get_sub_total(total, num_procs, rank):
    """
    Get the number of instances for the given rank.
    :param total: The total number of instances
    :param num_procs: The total number of processes
    :param rank: the current rank number
    :return: The number of instances for the given rank
    """
    return get_end_index(total, num_procs, rank) - get_start_index(total, num_procs, rank) + 1


def read_hourly_surfrad(folder, progress=True):
    """
    Reads only the hourly data from the daily data files from the input folder. It is assumed that all files in the
    folder belong to the same location.

    :param folder: A data folder with SURFRAD daily data files
    :param progress: Whether to show a progress bar
    :return: Hourly data frame and the meta information
    """

    # List all data files in the folder
    folder = path.expanduser(folder)
    all_files = [path.join(folder, file) for file in listdir(folder)]

    # Initialize variables to hold data
    df_list = []
    meta = None

    if progress:
        print("Reading data from folder {} ...".format(folder))

    # Initialize a progress bar
    pbar = IncrementalBar("Reading files", max=len(all_files))
    pbar.suffix = '%(percent).1f%% - %(eta)ds'

    # Read each daily file as a data frame and append it to the data frame list
    for file in all_files:

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
            pbar.next()

    if progress:
        pbar.finish()

    # Row bind all data frames to get the hourly data for the entire year
    hourly_data = pd.concat(df_list)

    return hourly_data, meta


def simulate_sun_positions(days, lead_times, lats, lons, solar_position_method="nrel_numpy", silent=False):
    """
    Simulations from each input sun position and each input time.
    """
    assert (len(lats) == len(lons)), "Numbers of latitudes and longitudes are not consistent"

    # The dimension of the simulation
    num_days = len(days)
    num_lead_times = len(lead_times)
    num_stations = len(lons)

    # Initialize a progress bar
    if not silent:
        pbar = IncrementalBar("Sun position silumations", max=num_stations * num_days)
        pbar.suffix = '%(percent).1f%% - %(eta)ds'

    # Initialize output variables
    keys = ["dni_extra", "air_mass", "zenith", "apparent_zenith", "azimuth"]
    sky_dict = {key: np.zeros((num_lead_times, num_days, num_stations)) for key in keys}

    for station_index in range(num_stations):

        # Determine the current location
        current_location = location.Location(latitude=lats[station_index], longitude=lons[station_index])

        for day_index in range(num_days):
            for lead_time_index in range(num_lead_times):

                # Determine the current time
                current_posix = days[day_index] + lead_times[lead_time_index]
                current_time = pd.Timestamp(current_posix, tz="UTC", unit='s')

                # Calculate sun position
                solar_position = current_location.get_solarposition(current_time, method=solar_position_method)

                # Calculate extraterrestrial DNI
                sky_dict["dni_extra"][lead_time_index, day_index, station_index] = \
                    irradiance.get_extra_radiation(current_time)

                # Calculate air mass
                sky_dict["air_mass"][lead_time_index, day_index, station_index] =\
                    atmosphere.get_relative_airmass(solar_position["apparent_zenith"])

                # Store other keys
                sky_dict["zenith"][lead_time_index, day_index, station_index] = solar_position["zenith"]
                sky_dict["apparent_zenith"][lead_time_index, day_index, station_index] = solar_position["apparent_zenith"]
                sky_dict["azimuth"][lead_time_index, day_index, station_index] = solar_position["azimuth"]

            if not silent:
                pbar.next()

    if not silent:
        pbar.finish()

    return sky_dict


def simulate_power(ghi_arr, tamb_arr, wspd_arr, albedo_arr, days, lead_times, sky_dict,
                   surface_tilt, surface_azimuth, pv_module, tcell_model_parameters, silent):
    """
    Simulates PV energy production
    """

    # Determine the dimensions
    num_analogs = ghi_arr.shape[0]
    num_lead_times = ghi_arr.shape[1]
    num_days = ghi_arr.shape[2]
    num_stations = ghi_arr.shape[3]

    # Initialization
    p_mp = np.zeros((num_analogs, num_lead_times, num_days, num_stations))

    if not silent:
        pbar = IncrementalBar("Power scenario simulation", max=num_stations * num_days)
        pbar.suffix = '%(percent).1f%% - %(eta)ds'

    for station_index in range(num_stations):
        for day_index in range(num_days):
            for lead_time_index in range(num_lead_times):
                for analog_index in range(num_analogs):

                    ghi = ghi_arr[analog_index, lead_time_index, day_index, station_index]

                    if ghi != 0:
                        albedo = albedo_arr[analog_index, lead_time_index, day_index, station_index]
                        wspd = wspd_arr[analog_index, lead_time_index, day_index, station_index]
                        tamb = tamb_arr[analog_index, lead_time_index, day_index, station_index]
                        air_mass = sky_dict["air_mass"][lead_time_index, day_index, station_index]
                        dni_extra = sky_dict["dni_extra"][lead_time_index, day_index, station_index]
                        zenith = sky_dict["zenith"][lead_time_index, day_index, station_index]
                        apparent_zenith = sky_dict["apparent_zenith"][lead_time_index, day_index, station_index]
                        azimuth = sky_dict["azimuth"][lead_time_index, day_index, station_index]

                        # Determine the current time
                        current_posix = days[day_index] + lead_times[lead_time_index]
                        current_time = pd.Timestamp(current_posix, tz="UTC", unit='s')

                        # Decompose DNI from GHI
                        dni_dict = irradiance.disc(ghi, zenith, current_time)
                        dni = dni_dict["dni"]

                        # Calculate POA sky diffuse
                        poa_sky_diffuse = irradiance.haydavies(
                            surface_tilt, surface_azimuth, ghi, dni, dni_extra, apparent_zenith, azimuth)

                        # Calculate POA ground diffuse
                        poa_ground_diffuse = irradiance.get_ground_diffuse(surface_tilt, ghi, albedo)

                        # Calculate angle of incidence
                        aoi = irradiance.aoi(surface_tilt, surface_azimuth, apparent_zenith, azimuth)

                        # Calculate POA total
                        poa_irradiance = irradiance.poa_components(aoi, dni, poa_sky_diffuse, poa_ground_diffuse)

                        # Calculate cell temperature
                        tcell = pvsystem.temperature.sapm_cell(
                            poa_irradiance['poa_global'], tamb, wspd, tcell_model_parameters['a'],
                            tcell_model_parameters['b'], tcell_model_parameters["deltaT"])

                        # Calculate effective irradiance
                        effective_irradiance = pvsystem.sapm_effective_irradiance(
                            poa_irradiance['poa_direct'], poa_irradiance['poa_diffuse'], air_mass, aoi, pv_module)

                        # Calculate power
                        sapm_out = pvsystem.sapm(effective_irradiance, tcell, pv_module)

                        # Save output to numpy
                        p_mp[analog_index, lead_time_index, day_index, station_index] = sapm_out["p_mp"]

            if not silent:
                pbar.next()

    if not silent:
        pbar.finish()

    return p_mp
