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

import line_profiler

@profile
def run_pv_simulations_with_analogs(nc_file, variable_dict, scenarios, progress=True,
                                    early_stopping=False, max_num_stations=None):
    """
    Simulates the power output ensemble given a weather analog file and several scenarios.

    The workflow for PV simulation is referenced from the following pvlib tutorial:

    https://nbviewer.ipython.org/github/pvlib/pvlib-python/blob/master/docs/tutorials/forecast_to_power.ipynb

    :param nc_file: A weather analog file
    :param variable_dict: A variable dictionary
    :param scenarios: Scenarios to simulation
    :param progress: Whether to show progress information
    :param early_stopping: Stop the simulation earlier when profiling is engaged
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
            print("Serial I/O is used when only 1 process is created.")
    else:
        parallel_netcdf = True

    nc = Dataset(nc_file, "a", parallel=parallel_netcdf, comm=comm, info=MPI.Info())

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
    num_scenarios = scenarios.total_scenarios()

    if max_num_stations is not None:
        if max_num_stations < num_stations:
            num_stations = max_num_stations

            if rank == 0:
                print("max_num_stations is effectively set to {}".format(max_num_stations))

    if progress and rank == 0:
        dimension_info = "Total simulation dimensions:\n" + \
                         "-- {} scenarios\n".format(num_scenarios) + \
                         "-- {} stations\n".format(num_stations) + \
                         "-- {} test days\n".format(num_days) + \
                         "-- {} lead times\n".format(num_lead_times) + \
                         "-- {} analog memebrs\n".format(num_analogs) + \
                         "-- {} processes".format(num_procs)
        print(dimension_info)

    # Determine the chunk of stations allocated to this current process
    station_index_start = get_start_index(num_stations, num_procs, rank)
    station_index_end = get_end_index(num_stations, num_procs, rank)
    num_sub_stations = station_index_end - station_index_start
    if num_sub_stations == 1:
        msg = "Rank #{} only processes 1 station but at least 2 stations".format(rank)
        raise Exception(msg)

    # Extract variables for the current stations
    if progress:
        print("Rank #{} reading data for stations [{}, {})".format(rank, station_index_start, station_index_end))

    nc_ghi = nc_vars["ghi"][:, :, :, station_index_start:station_index_end]
    nc_albedo = nc_vars["alb"][:, :, :, station_index_start:station_index_end]
    nc_wspd = nc_vars["wspd"][:, :, :, station_index_start:station_index_end]
    nc_tamb = nc_vars["tamb"][:, :, :, station_index_start:station_index_end]

    if progress:
        print("Rank #{} finished reading data".format(rank))

    # Initialize the array dimensions
    array_dimensions = ("num_analogs", "num_flts", "num_test_times", "num_stations")

    # If profiling is used, I explicitly terminate the program earlier after certain simulations.
    early_stopping_count = 5

    # Initialize progress bar
    if early_stopping and rank == 0:
        bar_length = early_stopping_count

        msg = "Early stopping is engaged. Progress will be terminated after {} simulated instances".format(
            early_stopping_count * num_analogs * num_lead_times)
        msg += " including {} * {} lead times * {} analog members".format(
            early_stopping_count, num_lead_times, num_analogs)

        print(msg)

    else:
        bar_length = num_scenarios * num_stations * num_days

    pbar = IncrementalBar("PV simulation rank #{}".format(rank), max=bar_length)
    pbar.suffix = '%(percent).1f%% - %(eta)ds'

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
        nc_p_mp = nc_output_group.variables.get("p_mp")
        p_mp = np.zeros((num_analogs, num_lead_times, num_days, num_sub_stations))

        if nc_p_mp is None:
            nc_p_mp = nc_output_group.createVariable("p_mp", "f8", array_dimensions)

        nc_p_mp.long_name = "power at maximum-power point"

        # Copy values from the current scenarios
        surface_tilt = current_scenario["surface_tilt"]
        surface_azimuth = current_scenario["surface_azimuth"]
        tcell_model_parameters = temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"][
            current_scenario["tcell_model_parameters"]]
        pv_module = pvsystem.retrieve_sam("SandiaMod")[current_scenario["pv_module"]]

        # TODO: Uncomment this if you are using liujordan separation model
        # forecast_model = forecast.NAM()

        for station_index in range(num_sub_stations):

            # Determine the current location
            latitude = nc_vars["lat"][station_index].data
            longitude = nc_vars["lon"][station_index].data

            current_location = location.Location(latitude=latitude, longitude=longitude)

            for day_index in range(num_days):

                for lead_time_index in range(num_lead_times):

                    # Determine the current time
                    current_posix = nc_vars["date"][day_index] + nc_vars["flt"][lead_time_index]
                    current_time = pd.Timestamp(current_posix, tz="UTC", unit='s')

                    # Calculate sun position
                    #
                    # TODO: What is the relation between solar position and air pressure?
                    #
                    solar_position = current_location.get_solarposition(current_time)

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
                        nc_p_mp[analog_index, lead_time_index, day_index, station_index] = sapm_out["p_mp"]

                # Update the progress bar
                if progress and rank == 0:
                    pbar.next()

                # Subtract 1 from early stopping counter
                early_stopping_count -= 1

                if early_stopping and early_stopping_count == 0:
                    pbar.finish()
                    nc.close()

                    if progress and rank == 0:
                        print("Simulation terminated due to the profiling tool engaged")

                    return

    pbar.finish()
    nc.close()

    if progress and rank == 0:
        print("PV simulation is complete!")
    return


if __name__ == '__main__':

    # Get MPI info
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    welcome_msg = "evergreen -- our pursuit of a more sustainable future"
    nc_message = "an NetCDF file with weather analogs generated by PAnEn (https://weiming-hu.github.io/AnalogsEnsemble/"

    # Define arguments
    parser = argparse.ArgumentParser(description=welcome_msg)
    parser.add_argument('--nc', help=nc_message, required=False, default="~/data/analogs_PA.nc")
    parser.add_argument('--silent', help="No progress information", action='store_true', default=False)
    parser.add_argument('--profile', help="Turn on profiling", action='store_true', default=False)
    parser.add_argument('--profiler', default='pyinstrument', help="Either pyinstrument or yappi")
    parser.add_argument('--stop', help="Early stop [useful in testing]", action='store_true', default=False)
    parser.add_argument('--stations', default=None, type=int, 
                        help="Limit the number of stations to simulate [useful in testing].")

    # Parse arguments
    args = parser.parse_args()

    # This is my input NetCDF file generated from Analog Ensemble
    nc_file =  os.path.expanduser(args.nc)

    """
    Variable Mapping
    
    This is a dictionary of variables names in the PV simulator (keys) and the Analog Ensemble weather file (values)
    
    Please modify this dictionary if you have different names defined in the weather analog files.
    
    Do NOT modify the dictionary keys as they will be indexed internally.
    
    Only change the values to the existing variable names in your weather analog files. 
    """
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

    """
    Simulation Scenarios
    
    This class defines the scenarios to be simulated. Scenarios are defined by the combination of a single value
    from each of the keys.
    
    To modify scenarios, only change the values for existing keys. Currently, it does NOT support adding or
    removing keys.
    
    Values must be a sequence-list type, e.g. a list or a tuple.
    """
    # Initialization. Do not remove.
    scenarios = Scenarios()

    # Change values of keys as you want
    scenarios["surface_tilt"] = [0]
    scenarios["surface_azimuth"] = [180]
    scenarios["tcell_model_parameters"] = ["open_rack_glass_glass"]
    scenarios["pv_module"] = ["Silevo_Triex_U300_Black__2014_"]

    # Start a profiler
    if args.profile:
        args.stop = True

        if args.profiler == "yappi":
            import yappi
            yappi.start()

        elif args.profiler == "pyinstrument":
            from pyinstrument import Profiler
            profiler = Profiler()
            profiler.start()

        else:
            raise Exception("Unsupported profiler: {}".format(args.profiler))

    if not args.profile and not args.stop:
        # If no prifiling tools engaged and not early stopping,
        # this is in operation mode, and then use Cython optimization
        #
        import pyximport; pyximport.install()

    # Run the simulator
    run_pv_simulations_with_analogs(nc_file, variable_dict, scenarios,
                                    progress=not args.silent,
                                    early_stopping=args.stop,
                                    max_num_stations=args.stations)

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
