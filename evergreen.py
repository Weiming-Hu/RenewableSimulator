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
from netCDF4 import Dataset
from pvlib import temperature, pvsystem

# Self hosted modules
from Scenarios import Scenarios
from Functions import simulate_sun_positions, simulate_power, get_start_index, get_end_index

# Performance add-on
from mpi4py import MPI


def run_pv_simulations_with_analogs(
        nc_file, variable_dict, scenarios, progress=True, solar_position_method="nrel_numpy",
        downscale=None, profile_memory=False):
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
        #num_days = math.ceil(num_days/downscale)
        #num_analogs = math.ceil(num_analogs/downscale)
        #num_scenarios = math.ceil(num_scenarios/downscale)

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

    # Extract variables for the current stations
    if progress:
        print("Rank #{} reading data from stations [{}, {})".format(rank, station_index_start, station_index_end))

    # These are high dimensional variables
    nc_ghi = nc.variables[variable_dict["ghi"]]
    nc_albedo = nc.variables[variable_dict["alb"]]
    nc_wspd = nc.variables[variable_dict["wspd"]]
    nc_tamb = nc.variables[variable_dict["tamb"]]

    # Actually read the subset of values
    nc_ghi = nc_ghi[0:num_analogs, 0:num_lead_times, 0:num_days, station_index_start:station_index_end]
    nc_albedo = nc_albedo[0:num_analogs, 0:num_lead_times, 0:num_days, station_index_start:station_index_end]
    nc_wspd = nc_wspd[0:num_analogs, 0:num_lead_times, 0:num_days, station_index_start:station_index_end]
    nc_tamb = nc_tamb[0:num_analogs, 0:num_lead_times, 0:num_days, station_index_start:station_index_end]

    # These are single dimensional vectors
    nc_lat = nc.variables[variable_dict["lat"]][station_index_start:station_index_end]
    nc_lon = nc.variables[variable_dict["lon"]][station_index_start:station_index_end]
    nc_day = nc.variables[variable_dict["date"]][0:num_days]
    nc_flt = nc.variables[variable_dict["flt"]][0:num_lead_times]

    # Weather output pre-process
    nc_albedo /= 100
    nc_tamb -= 273.15

    if progress:
        print("Rank #{} finished reading data".format(rank))

    # Pre-calculate air mass and extraterrestrial irradiance from solar positions
    sky_dict = simulate_sun_positions(nc_day, nc_flt, nc_lat, nc_lon, solar_position_method, rank)

    # Batch run simulations
    for scenario_index in range(num_scenarios):

        if progress:
            print("Rank #{} simulating scenario {}/{}".format(rank, scenario_index, num_scenarios))

        # Extract current scenario
        current_scenario = scenarios.get_scenario(scenario_index)

        # Create a group for the current scenario
        nc_output_group = nc.createGroup("PV_simulation_scenario_" + '{:05d}'.format(scenario_index))

        # Write the scenario to the group
        for key, value in current_scenario.items():
            nc_output_group.setncattr(key, value)

        # Create an array to store power at maximum-power point
        nc_p_mp = nc_output_group.variables.get("p_mp")

        if nc_p_mp is None:
            nc_p_mp = nc_output_group.createVariable(
                "p_mp", "f8", ("num_analogs", "num_flts", "num_test_times", "num_stations"))

        nc_p_mp.long_name = "power at maximum-power point"

        # Copy values from the current scenarios
        surface_tilt = current_scenario["surface_tilt"]
        surface_azimuth = current_scenario["surface_azimuth"]
        pv_module = pvsystem.retrieve_sam("SandiaMod")[current_scenario["pv_module"]]
        tcell_model_parameters = \
            temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"][current_scenario["tcell_model_parameters"]]

        # Simulate with the current scenario
        p_mp = simulate_power(nc_ghi, nc_tamb, nc_wspd, nc_albedo, nc_day, nc_flt, sky_dict,
                              surface_tilt, surface_azimuth, pv_module, tcell_model_parameters, rank)

        # Write the simulation results with the current scenario to the NetCDF file
        nc_p_mp[0:num_analogs, 0:num_lead_times, 0:num_days, station_index_start:station_index_end] = p_mp

    if progress and rank == 0:
        print("Power simulation is complete!")

        if profile_memory:
            print("Rank #{} heap usage:".format(rank))
            heap_usage = hpy().heap()
            print(heap_usage)
            print()
            print(heap_usage.more)

    nc.close()

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

        else:
            raise Exception("Unsupported profiler: {}".format(args.profiler))

    # Run the simulator
    run_pv_simulations_with_analogs(
        nc_file=nc_file, variable_dict=variable_dict, scenarios=scenarios, progress=not args.silent,
        solar_position_method=args.solar, downscale=args.downscale,
        profile_memory=(args.profiler == "memory"))

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

