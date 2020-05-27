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
import math
import argparse
import datetime

# Scientific python add-ons
import yaml
import numpy as np
from netCDF4 import Dataset
from pvlib import temperature

# Self hosted modules
from Scenarios import Scenarios
from Functions import simulate_sun_positions, simulate_power_batch, get_start_index, get_end_index

# Performance add-on
from mpi4py import MPI


def run_pv_simulations_with_analogs(
        nc_file, variable_dict, scenarios, progress=True, solar_position_method="nrel_numpy",
        downscale=None, profile_memory=False, simple_clock=False):
    """
    Simulates the power output ensemble given a weather analog file and several scenarios.

    The workflow for PV simulation is referenced from the following pvlib tutorial:

    https://nbviewer.ipython.org/github/pvlib/pvlib-python/blob/master/docs/tutorials/forecast_to_power.ipynb

    :param nc_file: A weather analog file
    :param variable_dict: A variable dictionary
    :param scenarios: Scenarios to simulation
    :param progress: Whether to show progress information
    :param solar_position_method: The method to use for calculating solar position
    :param downscale: The denominator to downscale the computation for testing purposes
    :return: None
    """

    # Sanity checks
    if not isinstance(scenarios, Scenarios):
        raise Exception("Input scenarios should an object of the class Scenarios")

    if not os.path.isfile(nc_file):
        raise Exception("File {} not exists".format(nc_file))

    # Get the current MPI process information
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_procs = comm.Get_size()

    if progress and rank == 0:
        print("Running PV simulation with AnEn ...")

    if rank == 0 and simple_clock:
        timestamps = [time()]
        log_names = []

    # Open the NetCDF file
    if num_procs == 1:
        parallel_netcdf = False
        if progress:
            print("No parallel I/O is enabled when only 1 process is created.")
    else:
        parallel_netcdf = True

    if progress and rank == 0:
        print("Input file is {}".format(nc_file))
    nc = Dataset(nc_file, "a", parallel=parallel_netcdf)

    # Determine the dimensions of the problem to simulate
    num_stations = nc.dimensions['num_stations'].size
    num_days = nc.dimensions['num_test_times'].size
    num_lead_times = nc.dimensions['num_flts'].size
    num_analogs = nc.dimensions['num_analogs'].size
    num_scenarios = scenarios.total_scenarios()

    # Downscale the computation if it is set to an integer number
    if isinstance(downscale, int):
        num_stations = math.ceil(num_stations/downscale)

    if progress and rank == 0:
        summary_info = "Summary information for the current task:\n" + \
                         "-- {} scenarios\n".format(num_scenarios) + \
                         "-- {} stations\n".format(num_stations) + \
                         "-- {} test days\n".format(num_days) + \
                         "-- {} lead times\n".format(num_lead_times) + \
                         "-- {} analog members\n".format(num_analogs) + \
                         "-- {} processes".format(num_procs)
        print(summary_info)

    # Determine the chunk of stations assigned to this current process
    station_index_start = get_start_index(num_stations, num_procs, rank)
    station_index_end = get_end_index(num_stations, num_procs, rank)


    """
    Read analog variables
    """

    if progress: print("Rank #{} reading analogs from stations [{}, {})".format(
        rank, station_index_start, station_index_end))

    # These are high dimensional variables
    ghi_anen = nc.variables[variable_dict["ghi"]]
    albedo_anen = nc.variables[variable_dict["alb"]]
    wspd_anen = nc.variables[variable_dict["wspd"]]
    tamb_anen = nc.variables[variable_dict["tamb"]]

    # Set collective mode for better I/O
    ghi_anen.set_collective(True)
    albedo_anen.set_collective(True)
    wspd_anen.set_collective(True)
    tamb_anen.set_collective(True)

    # Actually read the subset of values
    ghi_anen = ghi_anen[0:num_analogs, 0:num_lead_times, 0:num_days, station_index_start:station_index_end]
    albedo_anen = albedo_anen[0:num_analogs, 0:num_lead_times, 0:num_days, station_index_start:station_index_end]
    wspd_anen = wspd_anen[0:num_analogs, 0:num_lead_times, 0:num_days, station_index_start:station_index_end]
    tamb_anen = tamb_anen[0:num_analogs, 0:num_lead_times, 0:num_days, station_index_start:station_index_end]

    # Weather output pre-process
    albedo_anen /= 100
    tamb_anen -= 273.15


    """
    Read forecast variables
    """

    if progress: print("Rank #{} reading forecasts from stations [{}, {})".format(
        rank, station_index_start, station_index_end))

    # These are high dimensional variables
    nc_data = nc.groups["Forecasts"].variables["Data"]

    # Set collective mode for better I/O
    nc_data.set_collective(True)

    # Figure out which variables to read
    ghi_index, = np.where(nc.groups["Forecasts"].variables["ParameterNames"][:] == variable_dict["ghi_fcsts"])
    albedo_index, = np.where(nc.groups["Forecasts"].variables["ParameterNames"][:] == variable_dict["alb_fcsts"])
    wspd_index, = np.where(nc.groups["Forecasts"].variables["ParameterNames"][:] == variable_dict["wspd_fcsts"])
    tamb_index, = np.where(nc.groups["Forecasts"].variables["ParameterNames"][:] == variable_dict["tamb_fcsts"])

    if not (len(ghi_index) == 1 and len(albedo_index) == 1 and len(wspd_index) == 1 and len(tamb_index) == 1):
        raise Exception("Failed to find some forecast variables from file {}".format(nc_file))

    # Actually read the subset of values
    ghi_fcsts = nc_data[0:num_lead_times, 0:num_days, station_index_start:station_index_end, ghi_index]
    albedo_fcsts = nc_data[0:num_lead_times, 0:num_days, station_index_start:station_index_end, albedo_index]
    wspd_fcsts = nc_data[0:num_lead_times, 0:num_days, station_index_start:station_index_end, wspd_index]
    tamb_fcsts = nc_data[0:num_lead_times, 0:num_days, station_index_start:station_index_end, tamb_index]

    # Weather output pre-process
    albedo_fcsts /= 100
    tamb_fcsts -= 273.15

    # Reorganize the dimensions
    ghi_fcsts = np.transpose(ghi_fcsts, (3, 0, 1, 2))
    albedo_fcsts = np.transpose(albedo_fcsts, (3, 0, 1, 2))
    wspd_fcsts = np.transpose(wspd_fcsts, (3, 0, 1, 2))
    tamb_fcsts = np.transpose(tamb_fcsts, (3, 0, 1, 2))


    """
    Read meta variables
    """

    if progress: print("Rank #{} reading meta variables from stations [{}, {})".format(
        rank, station_index_start, station_index_end))

    # These are single dimensional vectors
    nc_lat = nc.variables[variable_dict["lat"]][station_index_start:station_index_end]
    nc_lon = nc.variables[variable_dict["lon"]][station_index_start:station_index_end]
    nc_day = nc.variables[variable_dict["date"]][0:num_days]
    nc_flt = nc.variables[variable_dict["flt"]][0:num_lead_times]


    """
    Read observation values
    """

    if progress: print("Rank #{} reading observations from stations [{}, {})".format(
        rank, station_index_start, station_index_end))

    # These are high dimensional variables
    nc_data = nc.groups["Observations"].variables["Data"]
    obs_times = nc.groups["Observations"].variables["Times"]

    # Set collective mode for better I/O
    nc_data.set_collective(True)
    obs_times.set_collective(True)

    # Figure out which variables to read
    ghi_index, = np.where(nc.groups["Observations"].variables["ParameterNames"][:] == variable_dict["ghi_obs"])
    albedo_index, = np.where(nc.groups["Observations"].variables["ParameterNames"][:] == variable_dict["alb_obs"])
    wspd_index, = np.where(nc.groups["Observations"].variables["ParameterNames"][:] == variable_dict["wspd_obs"])
    tamb_index, = np.where(nc.groups["Observations"].variables["ParameterNames"][:] == variable_dict["tamb_obs"])

    if not (len(ghi_index) == 1 and len(albedo_index) == 1 and len(wspd_index) == 1 and len(tamb_index) == 1):
        raise Exception("Failed to find some observation variables from file {}".format(nc_file))

    # Reshape data
    nc_data = nc_data[:, station_index_start:station_index_end,
              (ghi_index[0], albedo_index[0], wspd_index[0], tamb_index[0])]
    obs_times = obs_times[:]

    nc_data_reshape = np.zeros((num_lead_times, num_days, nc_data.shape[1], nc_data.shape[2]))
    for lead_time_index in range(nc_data_reshape.shape[0]):
        for day_index in range(nc_data_reshape.shape[1]):
            obs_time_index, = np.where(obs_times == nc_flt[lead_time_index] + nc_day[day_index])
            if len(obs_time_index) == 1:
                nc_data_reshape[lead_time_index, day_index, :, :] = nc_data[obs_time_index, :, :]
            else:
                nc_data_reshape[lead_time_index, day_index, :, :] = np.nan

    nc_data_reshape = np.transpose(nc_data_reshape, (3, 0, 1, 2))

    # Actually read the subset of values
    ghi_obs = nc_data_reshape[0, :, :, :]
    albedo_obs = nc_data_reshape[1, :, :, :]
    wspd_obs = nc_data_reshape[2, :, :, :]
    tamb_obs = nc_data_reshape[3, :, :, :]

    # Weather output pre-process
    albedo_obs /= 100
    tamb_obs -= 273.15

    if progress:
        print("Rank #{} finished reading data".format(rank))

    if rank == 0 and simple_clock:
        timestamps.append(time())
        log_names.append('Reading file')

    # Pre-calculate air mass and extraterrestrial irradiance from solar positions
    sky_dict = simulate_sun_positions(nc_day, nc_flt, nc_lat, nc_lon, solar_position_method, rank)

    if rank == 0 and simple_clock:
        timestamps.append(time())
        log_names.append('Sun position calculation')

    # Batch run simulations
    simulate_power_batch(
        p_mp_varname="analogs", p_mp_longname="Maximum power simulation from analogs",
        num_scenarios=num_scenarios, num_analogs=num_analogs, num_lead_times=num_lead_times, num_days=num_days,
        nc=nc, nc_ghi=ghi_anen, nc_tamb=tamb_anen, nc_wspd=wspd_anen,
        nc_albedo=albedo_anen, nc_day=nc_day, nc_flt=nc_flt,
        sky_dict=sky_dict, temperature=temperature, scenarios=scenarios,
        simple_clock=simple_clock, timestamps=timestamps, log_names=log_names,
        station_index_start=station_index_start, station_index_end=station_index_end,
        rank=rank, progress=progress)

    simulate_power_batch(
        p_mp_varname="forecasts", p_mp_longname="Maximum power simulation from WRF NAM forecasts",
        num_scenarios=num_scenarios, num_analogs=1, num_lead_times=num_lead_times, num_days=num_days,
        nc=nc, nc_ghi=ghi_fcsts, nc_tamb=tamb_fcsts, nc_wspd=wspd_fcsts,
        nc_albedo=albedo_fcsts, nc_day=nc_day, nc_flt=nc_flt,
        sky_dict=sky_dict, temperature=temperature, scenarios=scenarios,
        simple_clock=simple_clock, timestamps=timestamps, log_names=log_names,
        station_index_start=station_index_start, station_index_end=station_index_end,
        rank=rank, progress=progress)

    simulate_power_batch(
        p_mp_varname="analysis", p_mp_longname="Maximum power simulation from WRF NAM analysis",
        num_scenarios=num_scenarios, num_analogs=1, num_lead_times=num_lead_times, num_days=num_days,
        nc=nc, nc_ghi=ghi_obs, nc_tamb=tamb_obs, nc_wspd=wspd_obs,
        nc_albedo=albedo_obs, nc_day=nc_day, nc_flt=nc_flt,
        sky_dict=sky_dict, temperature=temperature, scenarios=scenarios,
        simple_clock=simple_clock, timestamps=timestamps, log_names=log_names,
        station_index_start=station_index_start, station_index_end=station_index_end,
        rank=rank, progress=progress)

    if progress and rank == 0:
        print("Power simulation is complete!")

        if profile_memory:
            print("Rank #{} heap usage:".format(rank))
            heap_usage = hpy().heap()
            print(heap_usage)
            print()
            print(heap_usage.more)

    nc.close()

    if rank == 0 and simple_clock:
        if len(log_names) == len(timestamps) - 1:
            # This is expected
            print("****************** Simple Clock Summary ******************")
            print("Total wall time: {}s".format(timestamps[-1] - timestamps[0]))
            for i in range(len(log_names)):
                print("{}: {}s".format(log_names[i], timestamps[i + 1] - timestamps[i]))
            print("*************** End of Simple Clock Summary **************")

        else:
            # This is not expected
            print("Something unexpected happened. The lengths are not correct. I'm dumping the results.")
            print(log_names)
            print(timestamps)

    return


if __name__ == '__main__':

    # Get MPI info
    rank = MPI.COMM_WORLD.rank

    welcome_msg = "evergreen -- our pursuit of a more sustainable future"
    nc_message = "An NetCDF file with weather analogs generated by PAnEn (https://weiming-hu.github.io/AnalogsEnsemble/"
    map_message = "A variable map that translates variable names in the NetCDF file to internally required names"
    scenario_message = "A dictionary defines key values for scenarios"

    # Define arguments
    parser = argparse.ArgumentParser(description=welcome_msg)
    parser.add_argument('--nc', help=nc_message, required=False, default="analogs_PA.nc")
    parser.add_argument('--map', help=map_message, required=False, default="variable-map.yaml")
    parser.add_argument('--scenario', help=scenario_message, required=False, default="scenarios.yaml")
    parser.add_argument('--silent', help="No progress information", action='store_true', default=False)
    parser.add_argument('--solar', help="Method for solar position calculation", default="nrel_numpy")
    parser.add_argument('--profile', help="Turn on profiling", action='store_true', default=False)
    parser.add_argument('--profiler', default='pyinstrument', help="Either pyinstrument or yappi")
    parser.add_argument('--downscale', default=None, type=int, help="Subset the computation for testing")

    # Parse arguments
    args = parser.parse_args()

    if not args.silent and rank == 0:
        print(welcome_msg)

    # This is my input NetCDF file generated from Analog Ensemble
    nc_file = os.path.expanduser(args.nc)

    """
    Variable Mapping
    
    This is a dictionary of variables names in the PV simulator (keys) and the Analog Ensemble weather file (values)
    
    Please modify through the yaml file.
    
    Do NOT modify the dictionary keys as they will be indexed internally.
    
    Only change the values to the existing variable names in your weather analog files. 
    """
    with open(os.path.expanduser(args.map)) as f:
        variable_dict = yaml.load(f, Loader=yaml.FullLoader)


    """
    Simulation Scenarios
    
    This class defines the scenarios to be simulated. Scenarios are defined by the combination of a single value
    from each of the keys.
    
    To modify scenarios, only change, through the yaml file, the values for existing keys.
    Currently, it does NOT support adding or removing keys.
    
    Values must be a sequence-list type, e.g. a list or a tuple.
    """
    # Read the YAML file
    with open(os.path.expanduser(args.scenario)) as f:
        scenarios = yaml.load(f, Loader=yaml.FullLoader)

    # Convert dictionary to Scenarios
    scenarios = Scenarios(scenarios)

    # Start a profiler
    if args.profile:

        if args.profiler == "yappi":
            import yappi
            yappi.start()

        elif args.profiler == "pyinstrument":
            from pyinstrument import Profiler
            profiler = Profiler()
            profiler.start()

        elif args.profiler == "memory":
            from guppy import hpy

        elif args.profiler == "line_profiler":
            try:
                # Use the decorator from kernprof to profile the main simulation function
                run_pv_simulations_with_analogs = profile(run_pv_simulations_with_analogs)
            except:
                raise Exception("Failed with function decorating. Did you properly use kernprof ?")

        elif args.profiler == "simple":
            from time import time

        else:
            raise Exception("Unsupported profiler: {}".format(args.profiler))

    # Run the simulator
    run_pv_simulations_with_analogs(
        nc_file=nc_file, variable_dict=variable_dict, scenarios=scenarios, progress=not args.silent,
        solar_position_method=args.solar, downscale=args.downscale,
        profile_memory=(args.profile and args.profiler == "memory"),
        simple_clock=(args.profile and args.profiler == "simple"))

    if args.profile:
        if args.profiler == "yappi":
            current_time = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

            func_stats = yappi.get_func_stats()

            func_stats.save("yappi_{}_rank-{}.log".format(current_time, rank), type="callgrind")

            if rank == 0:
                print("You have chosen yappi as the profiler tool.")
                print("yappi_* files are created. You can use callgrind tools for visualization")

        elif args.profiler == "pyinstrument":
            profiler.stop()

            if rank == 0:
                print(profiler.output_text(unicode=True, color=True))

