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
import argparse
import datetime

# Scientific python add-ons
import yaml
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from pvlib import location, atmosphere, temperature, pvsystem, irradiance

# Visualization add-ons
from progress.bar import IncrementalBar

# Self hosted modules
from Functions import simulate_single_instance, get_start_index, get_end_index
from Scenarios import Scenarios

# Performance add-ons
from mpi4py import MPI


def run_pv_simulations_with_analogs(nc_file, variable_dict, scenarios, progress=True,
                                    max_num_stations=None, early_stopping=False,
                                    early_stopping_count=5):
    """
    Simulates the power output ensemble given a weather analog file and several scenarios.

    The workflow for PV simulation is referenced from the following pvlib tutorial:

    https://nbviewer.ipython.org/github/pvlib/pvlib-python/blob/master/docs/tutorials/forecast_to_power.ipynb

    :param nc_file: A weather analog file
    :param variable_dict: A variable dictionary
    :param scenarios: Scenarios to simulation
    :param progress: Whether to show progress information
    :param early_stopping: Stop the simulation earlier when profiling is engaged
    :param early_stopping_count: A small number to signify the early stopping of the simulation.
    :return:
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
        print("{} processes have been started for the simulation".format(num_procs))

    # Open the NetCDF file
    if num_procs == 1:
        parallel_netcdf = False
        if progress:
            print("No parallel I/O is enabled when only 1 process is created.")
    else:
        parallel_netcdf = True

    nc = Dataset(nc_file, "a", parallel=parallel_netcdf, comm=comm, info=MPI.Info())

    # Determine the dimensions of the problem to simulate
    num_stations = nc.dimensions['num_stations'].size
    num_days = nc.dimensions['num_test_times'].size
    num_lead_times = nc.dimensions['num_flts'].size
    num_analogs = nc.dimensions['num_analogs'].size
    num_scenarios = scenarios.total_scenarios()

    # If the max number of stations is defined, constrain the number of total stations to simulate
    if max_num_stations is not None and max_num_stations < num_stations:
        num_stations = max_num_stations
        if rank == 0:
            print("The total number of stations to simulate is effectively limited to {}".format(num_stations))

    if progress and rank == 0:
        dimension_info = "Total simulation dimensions:\n" + \
                         "-- {} scenarios\n".format(num_scenarios) + \
                         "-- {} stations\n".format(num_stations) + \
                         "-- {} test days\n".format(num_days) + \
                         "-- {} lead times\n".format(num_lead_times) + \
                         "-- {} analog members\n".format(num_analogs) + \
                         "-- {} processes".format(num_procs)
        print(dimension_info)

    # Determine the chunk of stations assigned to this current process
    station_index_start = get_start_index(num_stations, num_procs, rank)
    station_index_end = get_end_index(num_stations, num_procs, rank)
    num_sub_stations = station_index_end - station_index_start

    if num_sub_stations == 1:
        msg = "Rank #{} only processes 1 station but at least 2 stations. Reduce the number of processes.".format(rank)
        raise Exception(msg)

    # Extract variables for the current stations
    if progress:
        print("Rank #{} reading data from stations [{}, {})".format(rank, station_index_start, station_index_end))

    # These are high dimensional arrays
    nc_ghi = nc.variables[variable_dict["ghi"]][:, :, :, station_index_start:station_index_end]
    nc_albedo = nc.variables[variable_dict["alb"]][:, :, :, station_index_start:station_index_end]
    nc_wspd = nc.variables[variable_dict["wspd"]][:, :, :, station_index_start:station_index_end]
    nc_tamb = nc.variables[variable_dict["tamb"]][:, :, :, station_index_start:station_index_end]
    
    # These are single dimensional vectors
    nc_lat = nc.variables[variable_dict["lat"]][:]
    nc_lon = nc.variables[variable_dict["lon"]][:]
    nc_day = nc.variables[variable_dict["date"]][:]
    nc_flt = nc.variables[variable_dict["flt"]][:]

    if progress:
        print("Rank #{} finished reading data".format(rank))

    # Initialize the array dimensions
    array_dimensions = ("num_analogs", "num_flts", "num_test_times", "num_stations")

    # Initialize a progress bar
    if early_stopping and rank == 0:
        bar_length = early_stopping_count

        msg = "Early stopping is engaged. Progress will be terminated after {} simulated instances".format(
            early_stopping_count * num_analogs * num_lead_times)
        msg += " including {} * {} lead times * {} analog members".format(
            early_stopping_count, num_lead_times, num_analogs)
        print(msg)

    else:
        bar_length = num_scenarios * num_stations * num_days

    if progress:
        if num_procs == 1:
            pbar = IncrementalBar("PV simulation rank #{}".format(rank), max=bar_length)
            pbar.suffix = '%(percent).1f%% - %(eta)ds'
        else:
            progress_threshold = bar_length / 100
            progress_value = 0
            progress_count = 0

    # Batch run simulations
    for scenario_index in range(num_scenarios):

        # Extract current scenario
        current_scenario = scenarios.get_scenario(scenario_index)

        # Create a group for the current scenario
        nc_output_group = nc.createGroup("PV_simulation_scenario_" + '{:05d}'.format(scenario_index))

        # Write the scenario to the group
        for key, value in current_scenario.items():
            nc_output_group.setncattr(key, value)

        # Create an array to store power at maximum-power point
        p_mp = np.zeros((num_analogs, num_lead_times, num_days, num_sub_stations))
        nc_p_mp = nc_output_group.variables.get("p_mp")

        if nc_p_mp is None:
            nc_p_mp = nc_output_group.createVariable("p_mp", "f8", array_dimensions)

        nc_p_mp.long_name = "power at maximum-power point"

        # Copy values from the current scenarios
        surface_tilt = current_scenario["surface_tilt"]
        surface_azimuth = current_scenario["surface_azimuth"]
        tcell_model_parameters = temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"][
            current_scenario["tcell_model_parameters"]]
        pv_module = pvsystem.retrieve_sam("SandiaMod")[current_scenario["pv_module"]]

        for station_index in range(num_sub_stations):

            # Determine the current location
            current_location = location.Location(latitude=nc_lat[station_index], longitude=nc_lon[station_index])

            for day_index in range(num_days):
                for lead_time_index in range(num_lead_times):

                    # Determine the current time
                    current_posix = nc_day[day_index] + nc_flt[lead_time_index]
                    current_time = pd.Timestamp(current_posix, tz="UTC", unit='s')

                    # Calculate sun position
                    #
                    # TODO: What is the relation between solar position and air pressure?
                    #
                    solar_position = current_location.get_solarposition(current_time, method="nrel_numba")

                    # Calculate extraterrestrial DNI
                    #
                    # TODO: There are different models to estimate the extraterrestrial radiation.
                    #
                    dni_extra = irradiance.get_extra_radiation(current_time)

                    # Calculate air mass
                    air_mass = atmosphere.get_relative_airmass(solar_position["apparent_zenith"])

                    for analog_index in range(num_analogs):

                        # Extra weather forecasts
                        ghi = nc_ghi[analog_index, lead_time_index, day_index, station_index]
                        albedo = nc_albedo[analog_index, lead_time_index, day_index, station_index] / 100
                        wspd = nc_wspd[analog_index, lead_time_index, day_index, station_index]
                        tamb = nc_tamb[analog_index, lead_time_index, day_index, station_index] - 273.15

                        # Simulate a single instance power output
                        sapm_out = simulate_single_instance(
                            ghi, dni_extra, tamb, wspd, albedo, current_time, surface_tilt, surface_azimuth,
                            pv_module, air_mass, tcell_model_parameters, solar_position)

                        # Assign results
                        p_mp[analog_index, lead_time_index, day_index, station_index] = sapm_out["p_mp"][0]

                # Update the progress bar
                if progress:
                    if num_procs == 1:
                        pbar.next()
                    else:
                        progress_count += 1
                        if progress_count > progress_threshold:
                            progress_count = 0
                            progress_value += 1
                            print("Rank #{} finished {}% ...".format(rank, progress_value))

                # Subtract 1 from early stopping counter
                early_stopping_count -= 1

                if early_stopping and early_stopping_count == 0:
                    nc_p_mp[:, :, :, station_index_start:station_index_end] = p_mp
                    nc.close()

                    if num_procs == 1:
                        pbar.finish()

                    if progress and rank == 0:
                        print("Simulation terminated due to the profiling tool engaged")

                    return

        # Write the simulation results with the current scenario to the NetCDF file
        nc_p_mp[:, :, :, station_index_start:station_index_end] = p_mp

    nc.close()

    if num_procs == 1:
        pbar.finish()

    if progress and rank == 0:
        print("PV simulation is complete!")
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
    parser.add_argument('--nc', help=nc_message, required=False, default="~/data/analogs_PA.nc")
    parser.add_argument('--map', help=map_message, required=False, default="variable-map.yaml")
    parser.add_argument('--scenario', help=scenario_message, required=False, default="scenarios.yaml")
    parser.add_argument('--silent', help="No progress information", action='store_true', default=False)
    parser.add_argument('--profile', help="Turn on profiling", action='store_true', default=False)
    parser.add_argument('--profiler', default='pyinstrument', help="Either pyinstrument or yappi")
    parser.add_argument('--stations', default=None, type=int, 
                        help="Limit the number of stations to simulate [useful in testing].")

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

        elif args.profiler == "line_profiler":
            try:
                # Use the decorator from kernprof to profile the main simulation function
                run_pv_simulations_with_analogs = profile(run_pv_simulations_with_analogs)
            except:
                raise Exception("Failed with function decorating. Did you properly use kernprof ?")

        else:
            raise Exception("Unsupported profiler: {}".format(args.profiler))

    # Run the simulator
    run_pv_simulations_with_analogs(
        nc_file=nc_file, variable_dict=variable_dict, scenarios=scenarios,
        progress=not args.silent, max_num_stations=args.stations, early_stopping=args.profile)

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
