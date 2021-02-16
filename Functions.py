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

from tqdm import tqdm
from netCDF4 import Dataset
from functools import partial
from tqdm.contrib.concurrent import process_map
from pvlib import pvsystem, irradiance, location, atmosphere, temperature


##############
# Simulation #
##############

def simulate_sun_positions_by_station(station_index, solar_position_method, days, lead_times, latitudes, longitudes):
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
    :return: A list with DNI, air mass, zenith, apparent zenith, and azimuth.
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
                           disable_progress_bar=False, cores=1, verbose=True):
    """
    Simulate sun positions

    :param days: A sequence of days in UNIX time format
    :param lead_times: A sequence of lead times in UNIX time format
    :param latitudes: A sequence of latitudes
    :param longitudes: A sequence of longitudes
    :param solar_position_method: The method to use for calculating solar positions
    :param disable_progress_bar: Whether to hide the progress bar
    :param cores: The number of cores to use
    :param verbose: Whether to be verbose
    :return: A disctionary with simulated results
    """
    assert (len(latitudes) == len(longitudes)), "Numbers of latitudes and longitudes are not consistent"

    # Define a simple wrapper
    wrapper = partial(simulate_sun_positions_by_station, days=days, lead_times=lead_times,
                      latitudes=latitudes, longitudes=longitudes, solar_position_method=solar_position_method)

    # parallel processing
    if verbose:
        print('Calculating sky conditions ...')

    if cores == 1:
        results = [None] * len(latitudes)
        for station_index in tqdm(range(len(latitudes)), disable=disable_progress_bar):
            results[station_index] = wrapper(station_index)
    else:
        results = process_map(wrapper, range(len(latitudes)), max_workers=cores, disable=disable_progress_bar,
                              chunksize=1 if len(latitudes) < 1000 else int(len(latitudes) / 100))

    # Initialize output variables
    sky_dict = {
        "dni_extra": np.stack([result[0] for result in results], axis=2),
        "air_mass": np.stack([result[1] for result in results], axis=2),
        "zenith": np.stack([result[2] for result in results], axis=2),
        "apparent_zenith": np.stack([result[3] for result in results], axis=2),
        "azimuth": np.stack([result[4] for result in results], axis=2)
    }

    return sky_dict


def simulate_power_by_station(
        station_index, surface_tilt, surface_azimuth, pv_module, tcell_model_parameters,
        ghi, tamb, wspd, albedo, days, lead_times, air_mass, dni_extra, zenith, apparent_zenith, azimuth):
    """
    This is the worker function for simulating power at a specified location. This function should be used inside
    of `simulate_power` and direct usage is discouraged.

    :param station_index: A station index
    :param ghi: See `simulate_power`
    :param tamb: See `simulate_power`
    :param wspd: See `simulate_power`
    :param albedo: See `simulate_power`
    :param days: See `simulate_power`
    :param lead_times: See `simulate_power`
    :param air_mass: See `simulate_power`
    :param dni_extra: See `simulate_power`
    :param zenith: See `simulate_power`
    :param apparent_zenith: See `simulate_power`
    :param azimuth: See `simulate_power`
    :param surface_tilt: See `simulate_power`
    :param surface_azimuth: See `simulate_power`
    :param pv_module: A PV module name
    :param tcell_model_parameters: A cell module name
    :return: A list with power, cell temperature, and the effective irradiance
    """

    # Sanity check
    assert 0 <= station_index < ghi.shape[3], 'Invalid station index'

    # Determine the dimensions
    num_analogs = ghi.shape[0]
    num_lead_times = ghi.shape[1]
    num_days = ghi.shape[2]

    # Initialization
    p_mp = np.zeros((num_analogs, num_lead_times, num_days))
    tcell = np.zeros((num_analogs, num_lead_times, num_days))
    effective_irradiance = np.zeros((num_analogs, num_lead_times, num_days))
    pv_module = pvsystem.retrieve_sam("SandiaMod")[pv_module]
    tcell_model_parameters = temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"][tcell_model_parameters]

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

                ##########################################################################################
                #                                                                                        #
                #                     Core procedures of simulating power at one location                #
                #                                                                                        #
                ##########################################################################################

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
                tcell[analog_index, lead_time_index, day_index] = pvsystem.temperature.sapm_cell(
                    poa_irradiance['poa_global'], tamb_, wspd_, tcell_model_parameters['a'],
                    tcell_model_parameters['b'], tcell_model_parameters["deltaT"])

                # Calculate effective irradiance
                effective_irradiance[analog_index, lead_time_index, day_index] = pvsystem.sapm_effective_irradiance(
                    poa_irradiance['poa_direct'], poa_irradiance['poa_diffuse'], air_mass_, aoi, pv_module)

                # Calculate power
                sapm_out = pvsystem.sapm(effective_irradiance[analog_index, lead_time_index, day_index],
                                         tcell[analog_index, lead_time_index, day_index], pv_module)

                # Save output to numpy
                p_mp[analog_index, lead_time_index, day_index] = sapm_out["p_mp"]

    return [p_mp, tcell, effective_irradiance]


def simulate_power(group_name, scenarios, nc,
                   ghi, tamb, wspd, alb, days, lead_times,
                   air_mass, dni_extra, zenith, apparent_zenith, azimuth,
                   parallel_nc=False, cores=1, verbose=True, output_stations_index=None,
                   disable_progress_bar=False, timer=None, skip_existing_scenario=False):
    """
    Simulate power and write to a specific group in the NetCDF file.

    :param group_name: The group name to be created under the scenario group
    :param scenarios: The scenarios to simulate
    :param nc: An opened Dataset from netCDF4 with write access
    :param ghi: Golbal horizontal irradiance
    :param tamb: Ambient temperature
    :param wspd: wind speed
    :param alb: albedo
    :param days: A sequence of test days in UNIX time format
    :param lead_times: A sequence of lead times in UNIX time format
    :param air_mass: Air mass from simulated sun positions
    :param dni_extra: Extraterrestrial direct normal irradiance from simulated sun positions
    :param zenith: Zenith from simulated sun positions
    :param apparent_zenith: Apparent zenith from simulated sun positions
    :param azimuth: Azimuth
    :param parallel_nc: Whether to have parallel access to NetCDF
    :param cores: Number of cores to use
    :param verbose: Whether to be verbose
    :param output_stations_index: Station indices used when writing output
    :param disable_progress_bar: Whether to hide progress bars
    :param timer: A simple timer
    :param skip_existing_scenario: Whether to skip simulation if a scenario already exists
    """

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

        timer.start('Simulate scenario {:05d} for {}'.format(scenario_index, group_name))

        if verbose:
            print("Simulating scenario {}/{} with sub-group name '{}'".format(
                scenario_index + 1, num_scenarios, group_name))

        # Extract current scenario
        current_scenario = scenarios.get_scenario(scenario_index)

        # Create a group for the current scenario
        scenario_name = "PV_simulation_scenario_" + '{:05d}'.format(scenario_index)

        if skip_existing_scenario:
            if scenario_name in nc.groups:
                if group_name in nc.groups[scenario_name].groups:

                    if verbose:
                        print("Skip simulating {} for scenario {}/{}".format(group_name, scenario_index + 1, num_scenarios))

                    timer.stop()
                    continue

        nc_scenario_group = nc.createGroup(scenario_name)

        # Write the scenario to the group
        for key, value in current_scenario.items():
            nc_scenario_group.setncattr(key, value)

        # Check whether I should add the dimension for single-member cases (e.g. forecasts and analysis)
        if num_analogs == 1:
            if 'single_member' not in nc_scenario_group.dimensions:
                nc_scenario_group.createDimension('single_member', size=1)

            output_dims = ("single_member", "num_flts", "num_test_times", "num_stations")

        else:
            output_dims = ("num_analogs", "num_flts", "num_test_times", "num_stations")

        # Create a wrapper function for this iteration
        wrapper = partial(simulate_power_by_station,
                          surface_tilt=current_scenario["surface_tilt"],
                          surface_azimuth=current_scenario["surface_azimuth"],
                          pv_module=current_scenario["pv_module"],
                          tcell_model_parameters=current_scenario["tcell_model_parameters"],
                          ghi=ghi, tamb=tamb, wspd=wspd, albedo=alb, days=days, lead_times=lead_times, air_mass=air_mass,
                          dni_extra=dni_extra, zenith=zenith, apparent_zenith=apparent_zenith, azimuth=azimuth)

        # Simulate with the current scenario
        if cores == 1:
            results = [None] * num_stations
            for station_index in tqdm(range(num_stations), disable=disable_progress_bar):
                results[station_index] = wrapper(station_index)
        else:
            results = process_map(wrapper, range(num_stations), max_workers=cores, disable=disable_progress_bar,
                                  chunksize=1 if num_stations < 1000 else int(num_stations / 100))

        results = {
            "power": np.stack([result[0] for result in results], axis=3),
            "tcell": np.stack([result[1] for result in results], axis=3),
            "effective_irradiance": np.stack([result[2] for result in results], axis=3)
        }

        timer.stop()
        timer.start('Write scenario {:05d}'.format(scenario_index))

        # Write results to the NetCDF file
        if verbose:
            print("Writing scenario {}/{}".format(scenario_index + 1, num_scenarios))

        write_array_dict(nc_scenario_group, group_name, results, output_dims, parallel_nc, output_stations_index)
        timer.stop()


#######################
# MPI parallelization #
#######################

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


############
# File I/O #
############

def read_yaml(file):
    with open(os.path.expanduser(file)) as f:
        content = yaml.load(f, Loader=yaml.FullLoader)
    return content


def get_nc_dim_size(nc_file, dim_name, parallel_nc):
    nc = Dataset(nc_file, 'r', parallel=parallel_nc)
    return nc.dimensions[dim_name].size


def write_array_dict(nc, group_name, d, dimensions, parallel_nc, output_stations_index=None):
    """
    Write a dictionary of arrays into a new group under the given nc device.

    :param nc: An NetCDF4 Dataset or Group
    :param group_name: Group name to be created
    :param d: A dictionary of arrays with the same dimensions
    :param dimensions: The dimension names to be associated in the NetCDF file for arrays
    :param parallel_nc: Whether to set collective access
    :param output_stations_index: The indices used when writing to the output variable. This is
    assumed to be the last dimension.
    """

    # Sanity check
    d_dims = [v.shape for v in d.values()]
    assert all([d_dims[0] == dim for dim in d_dims]), 'All arrays should have the same shape in the dictionary'
    assert len(d_dims[0]) == len(dimensions), 'The specified dimensions do not match array dimensions'

    # Create a group for the current scenario
    nc_group = nc.createGroup(group_name)

    # Write variables
    for k, v in d.items():

        # Get a variable device
        var = nc_group.variables.get(k)

        if var is None:
            var = nc_group.createVariable(k, "f8", dimensions)

        # Set parallel access
        if parallel_nc:
            var.set_collective(True)

        if output_stations_index is None:
            var[:] = v
        else:
            var[..., output_stations_index] = v


def recursive_summary_dict(d, prefix='--'):
    msg = ''

    for key, value in d.items():
        if isinstance(value, dict):
            msg += '{} {}:\n'.format(prefix, key)
            msg += recursive_summary_dict(value, prefix + ' --')
        else:
            msg += '{} {}: shape {}\n'.format(prefix, key, value.shape)

    return msg


def read_array_dict(nc, group_name, parallel_nc, stations_index=None):
    """
    Read a dictionary of arrays

    :param nc: An NetCDF4 Dataset or Group
    :param group_name: Group name to be created
    :param parallel_nc: Whether to set collective access
    :param stations_index: The indices used when reading partial arrays. This is assumed to be the last dimension.
    :return: A dictionary of arrays
    """

    # Initialization
    d = {}

    # Get access to the group
    nc_group = nc.groups[group_name]

    # Read data
    for k, v in nc_group.variables.items():

        if parallel_nc:
            v.set_collective(True)

        if stations_index is None:
            d[k] = v[:]
        else:
            d[k] = v[..., stations_index]

    return d
