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

import os
import math
import yaml

import numpy as np
import pandas as pd

from netCDF4 import Dataset
from os import listdir, path
from functools import partial
from progress.bar import IncrementalBar
from tqdm.contrib.concurrent import process_map
from pvlib import pvsystem, irradiance, iotools, location, atmosphere, temperature


def read_yaml(file):
    with open(os.path.expanduser(file)) as f:
        content = yaml.load(f, Loader=yaml.FullLoader)
    return content


def get_nc_dim_size(nc_file, dim_name, parallel_nc):
    nc = Dataset(nc_file, 'r', parallel=parallel_nc)
    return nc.dimensions[dim_name].size


def simulate_sun_positions_by_station(station_index, days, lead_times, latitudes, longitudes, solar_position_method):
    """
    This is the worker function of calculating sun positions at a specified station index. This function should
    be called from the parallel version of this function, `simulate_sun_positions`.

    Direct use of this function is discouraged. Please use the parallel version of this fucntion,
    `simulate_sun_positions`.

    :param station_index: A station index for simulation
    :param days: See `simulate_sun_positions`
    :param lead_times: See `simulate_sun_positions`
    :param latitudes: See `simulate_sun_positions`
    :param longitudes: See `simulate_sun_positions`
    :param solar_position_method: See `simulate_sun_positions`
    :return:
    """

    # Initialization
    num_lead_times, num_days = len(lead_times), len(days)

    dni_extra = np.zeros((num_lead_times, num_days))
    air_mass = np.zeros((num_lead_times, num_days))
    zenith = np.zeros((num_lead_times, num_days))
    apparent_zenith = np.zeros((num_lead_times, num_days))
    azimuth = np.zeros((num_lead_times, num_days))

    # Determine the current location
    current_location = location.Location(latitude=latitudes[station_index], longitude=longitudes[station_index])

    for day_index in range(num_days):
        for lead_time_index in range(num_lead_times):

            # Determine the current time
            current_posix = days[day_index] + lead_times[lead_time_index]
            current_time = pd.Timestamp(current_posix, tz="UTC", unit='s')

            # Calculate sun position
            solar_position = current_location.get_solarposition(
                current_time, method=solar_position_method, numthreads=1)

            # Calculate extraterrestrial DNI
            dni_extra[lead_time_index, day_index] = irradiance.get_extra_radiation(current_time)

            # Calculate air mass
            air_mass[lead_time_index, day_index] = atmosphere.get_relative_airmass(solar_position["apparent_zenith"])

            # Store other keys
            zenith[lead_time_index, day_index] = solar_position["zenith"]
            apparent_zenith[lead_time_index, day_index] = solar_position["apparent_zenith"]
            azimuth[lead_time_index, day_index] = solar_position["azimuth"]

    return [dni_extra, air_mass, zenith, apparent_zenith, azimuth]


def simulate_sun_positions(days, lead_times, latitudes, longitudes,
                           solar_position_method="nrel_numpy",
                           disable_progress_bar=False, cores=1):
    """
    Simulations from each input sun position and each input time.
    """
    assert (len(latitudes) == len(longitudes)), "Numbers of latitudes and longitudes are not consistent"

    # Create a wrapper function for parallelization
    wrapper = partial(simulate_sun_positions_by_station,
                      days=days, lead_times=lead_times, latitudes=latitudes, longitudes=longitudes,
                      solar_position_method=solar_position_method)

    # parallel processing
    results = process_map(wrapper, range(len(latitudes)), max_workers=cores, disable=disable_progress_bar)

    # Initialize output variables
    sky_dict = {
        "dni_extra": np.stack([result[0] for result in results], axis=2),
        "air_mass": np.stack([result[1] for result in results], axis=2),
        "zenith": np.stack([result[2] for result in results], axis=2),
        "apparent_zenith": np.stack([result[3] for result in results], axis=2),
        "azimuth": np.stack([result[4] for result in results], axis=2)
    }

    return sky_dict


def simulate_power_by_station(station_index, ghi, tamb, wspd, albedo, days, lead_times,
                              air_mass, dni_extra, zenith, apparent_zenith, azimuth,
                              surface_tilt, surface_azimuth, pv_module, tcell_model_parameters):

    """
    Simulates PV energy production
    """

    # Sanity check
    assert 0 < station_index < ghi.shape[3], 'Invalid station index'

    # Determine the dimensions
    num_analogs = ghi.shape[0]
    num_lead_times = ghi.shape[1]
    num_days = ghi.shape[2]

    # Initialization
    p_mp = np.zeros((num_analogs, num_lead_times, num_days))

    for day_index in range(num_days):
        for lead_time_index in range(num_lead_times):

            # Determine the current time
            current_posix = days[day_index] + lead_times[lead_time_index]
            current_time = pd.Timestamp(current_posix, tz="UTC", unit='s')

            for analog_index in range(num_analogs):
                ghi_ = ghi[analog_index, lead_time_index, day_index, station_index]

                if ghi_ == 0:
                    continue

                albedo_ = albedo[analog_index, lead_time_index, day_index, station_index]
                wspd_ = wspd[analog_index, lead_time_index, day_index, station_index]
                tamb_ = tamb[analog_index, lead_time_index, day_index, station_index]
                air_mass_ = air_mass[lead_time_index, day_index, station_index]
                dni_extra_ = dni_extra[lead_time_index, day_index, station_index]
                zenith_ = zenith[lead_time_index, day_index, station_index]
                apparent_zenith_ = apparent_zenith[lead_time_index, day_index, station_index]
                azimuth_ = azimuth[lead_time_index, day_index, station_index]

                # Decompose DNI from GHI
                dni_dict = irradiance.disc(ghi_, zenith_, current_time)

                # Calculate POA sky diffuse
                poa_sky_diffuse = irradiance.haydavies(
                    surface_tilt, surface_azimuth, ghi_, dni_dict["dni"], dni_extra_, apparent_zenith_, azimuth_)

                # Calculate POA ground diffuse
                poa_ground_diffuse = irradiance.get_ground_diffuse(surface_tilt, ghi_, albedo_)

                # Calculate angle of incidence
                aoi = irradiance.aoi(surface_tilt, surface_azimuth, apparent_zenith_, azimuth_)

                # Calculate POA total
                poa_irradiance = irradiance.poa_components(aoi, dni_dict["dni"], poa_sky_diffuse, poa_ground_diffuse)

                # Calculate cell temperature
                tcell = pvsystem.temperature.sapm_cell(
                    poa_irradiance['poa_global'], tamb_, wspd_, tcell_model_parameters['a'],
                    tcell_model_parameters['b'], tcell_model_parameters["deltaT"])

                # Calculate effective irradiance
                effective_irradiance = pvsystem.sapm_effective_irradiance(
                    poa_irradiance['poa_direct'], poa_irradiance['poa_diffuse'], air_mass_, aoi, pv_module)

                # Calculate power
                sapm_out = pvsystem.sapm(effective_irradiance, tcell, pv_module)

                # Save output to numpy
                p_mp[analog_index, lead_time_index, day_index] = sapm_out["p_mp"]

    return p_mp


def simulate_power(power_varname, power_longname, scenarios, nc,
                   ghi, tamb, wspd, alb, days, lead_times,
                   air_mass, dni_extra, zenith, apparent_zenith, azimuth,
                   parallel_nc=False, cores=1, verbose=True, output_stations_index=None,
                   disable_progress_bar=False, timer=None, temperature=temperature):

    # Sanity check
    assert ghi.shape == tamb.shape, "ghi.shape != tamb.shape"
    assert ghi.shape == wspd.shape, "ghi.shape != wspd.shape"
    assert ghi.shape == alb.shape, "ghi.shape != albedo.shape"
    assert ghi.shape[1] == len(lead_times), "ghi.shape[1] != len(lead_times)"
    assert ghi.shape[2] == len(days), "ghi.shape[2] != len(days)"
    assert ghi.shape[1:4] == air_mass.shape, "ghi.shape[1:4] != air_mass.shape"
    assert ghi.shape[1:4] == dni_extra.shape, "ghi.shape[1:4] != dni_extra.shape"
    assert ghi.shape[1:4] == zenith.shape, "ghi.shape[1:4] != zenith.shape"
    assert ghi.shape[1:4] == apparent_zenith.shape, "ghi.shape[1:4] != apparent_zenith.shape"
    assert ghi.shape[1:4] == azimuth.shape, "ghi.shape[1:4] != azimuth.shape"

    num_scenarios = scenarios.total_scenarios()
    num_analogs = ghi.shape[0]
    num_stations = ghi.shape[3]

    if output_stations_index is None:
        output_stations_index = list(range(num_stations))
        assert len(output_stations_index) == nc.dimensions['num_stations'].size, \
            'Output stations index do not match station dimension'

    for scenario_index in range(num_scenarios):

        timer.start('Simulate scenario {:05d}'.format(scenario_index))

        if verbose:
            print("Simulating scenario {}/{}: {}".format(scenario_index + 1, num_scenarios, power_longname))

        # Extract current scenario
        current_scenario = scenarios.get_scenario(scenario_index)

        # Create a group for the current scenario
        nc_output_group = nc.createGroup("PV_simulation_scenario_" + '{:05d}'.format(scenario_index))

        # Write the scenario to the group
        for key, value in current_scenario.items():
            nc_output_group.setncattr(key, value)

        # Check whether I should add the dimension for single-member cases (e.g. forecasts and analysis)
        if num_analogs == 1:
            if 'single_member' not in nc_output_group.dimensions:
                nc_output_group.createDimension('single_member', size=1)

        # Create an array to store power at maximum-power point
        nc_power = nc_output_group.variables.get(power_varname)

        if nc_power is None:
            if num_analogs == 1:
                nc_power = nc_output_group.createVariable(
                    power_varname, "f8", ("single_member", "num_flts", "num_test_times", "num_stations"))
            else:
                nc_power = nc_output_group.createVariable(
                    power_varname, "f8", ("num_analogs", "num_flts", "num_test_times", "num_stations"))

        if parallel_nc:
            nc_power.set_collective(True)

        nc_power.long_name = power_longname

        # Create a wrapper function
        wrapper = partial(
            simulate_power_by_station, ghi=ghi, tamb=tamb, wspd=wspd, albedo=alb,
            days=days, lead_times=lead_times, air_mass=air_mass, dni_extra=dni_extra,
            zenith=zenith, apparent_zenith=apparent_zenith, azimuth=azimuth,
            surface_tilt=current_scenario["surface_tilt"],
            surface_azimuth=current_scenario["surface_azimuth"],
            pv_module=pvsystem.retrieve_sam("SandiaMod")[current_scenario["pv_module"]],
            tcell_model_parameters=temperature.TEMPERATURE_MODEL_PARAMETERS[
                "sapm"][current_scenario["tcell_model_parameters"]])

        # Simulate with the current scenario
        results = process_map(wrapper, range(num_stations), max_workers=cores, disable=disable_progress_bar)

        timer.stop()
        timer.start('Write scenario {:05d}'.format(scenario_index))

        # Write results to the NetCDF file
        for index in range(num_stations):
            nc_power[:, :, :, output_stations_index[index]] = results[index]

        timer.stop()

    return


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
