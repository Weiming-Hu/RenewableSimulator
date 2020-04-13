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

# Scientific python add-ons
import yaml
import numpy as np
import pandas as pd
from pvlib import location, atmosphere, temperature, pvsystem, irradiance

# Visualization add-ons
from progress.bar import IncrementalBar

# Self hosted modules
from Scenarios import Scenarios
from Functions import read_hourly_SURFRAD, simulate_single_instance


def run_pv_simulation_with_surfrad(output_file_prefix, scenarios, year_folder, progress=True):

    # Read hourly SURFRAD data
    yearly_data, meta = read_hourly_SURFRAD(year_folder, progress)

    # Determine location
    current_location = location.Location(latitude=meta["latitude"], longitude=meta["longitude"])

    # Determine the total number of scenarios
    num_scenarios = scenarios.total_scenarios()

    # Initialize a progress bar
    pbar = IncrementalBar("PV simulation", max=yearly_data.shape[0] * num_scenarios)
    pbar.suffix = '%(percent).1f%% - %(eta)ds'

    for scenario_index in range(num_scenarios):

        current_scenario = scenarios.get_scenario(scenario_index)

        # get the current scenario settings
        surface_tilt = current_scenario["surface_tilt"]
        surface_azimuth = current_scenario["surface_azimuth"]
        tcell_model_parameters = temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"][
            current_scenario["tcell_model_parameters"]]
        pv_module = pvsystem.retrieve_sam("SandiaMod")[current_scenario["pv_module"]]

        # Initialize empty list for simulated maximum power
        p_mp = ['MaximumPowerOutput']
        times = ['Time']

        for row_index in range(yearly_data.shape[0]):

            # Extract current row and values
            current_row = yearly_data.iloc[row_index]
            ghi = current_row['ghi']
            tamb = current_row['temp_air']
            wspd = current_row['wind_speed']

            albedo = 0
            if ghi != 0:
                albedo = current_row['uw_solar'] / current_row['ghi']

            # Get current time
            datetime_str = '/'.join([str(int(x)) for x in current_row[['year', 'month', 'day', 'hour', 'minute']]])
            current_time = pd.to_datetime(datetime_str, format="%Y/%m/%d/%H/%M")
            times.append(current_time.strftime(format = "%Y-%m-%d-%H-%M"))

            # Calculate sun position
            solar_position = current_location.get_solarposition(current_time)

            # Calculate extraterrestrial DNI
            dni_extra = irradiance.get_extra_radiation(current_time)

            # Calculate air mass
            air_mass = atmosphere.get_relative_airmass(solar_position["apparent_zenith"])

            # Simulate a single instance power output
            sapm_out = simulate_single_instance(
                ghi, dni_extra, tamb, wspd, albedo, current_time, surface_tilt, surface_azimuth,
                pv_module, air_mass, tcell_model_parameters, solar_position)

            # Store results
            p_mp.append(sapm_out["p_mp"][0])

            if progress:
                pbar.next()

        np.savetxt(output_file_prefix + "_scenario-{:05d}.csv".format(scenario_index),
                   [p for p in zip(times, p_mp)], delimiter=',', fmt='%s')

    if progress:
        pbar.finish()

    return meta


if __name__ == '__main__':

    scenario_message = "A dictionary defines key values for scenarios"
    welcome_msg = "Simulate PV power supply with data from SURFRAD as ground truth. " + \
        "To successfully run this program, you need to provide the folder root to your data." + \
        "It is assumed that SURFRAD data files (*.dat) are organized in <root folder>/<location>/<year>, " + \
        "e.g. /Volumes/WD3TB/SURFRAD/extracted/Boulder_CO/2018"

    # Define arguments
    parser = argparse.ArgumentParser(description=welcome_msg)
    parser.add_argument('--root', help="Root data folder", required=False, default="/Volumes/WD3TB/SURFRAD/extracted")
    parser.add_argument('--output', help="Output data folder", required=False, default="SURFRAD")
    parser.add_argument('--scenario', help=scenario_message, required=False, default="scenarios.yaml")
    parser.add_argument('--silent', help="No progress information", action='store_true', default=False)

    # Parse arguments
    args = parser.parse_args()

    # Check output folder
    if os.path.exists(args.output):
        if os.path.isdir(args.output):
            print("Output folder {} exists. Contents will be overwritten!".format(args.output))
        else:
            raise Exception('{} is not a folder'.format(args.output))
    else:
        os.makedirs(args.output)

    # Read YAML files
    with open(args.scenario) as f:
        scenarios = yaml.load(f, Loader=yaml.FullLoader)

    # Create an object of the class Scenarios
    scenarios = Scenarios(scenarios)

    # Expand user folder if needed
    root_folder = os.path.expanduser(args.root)

    # Walk through each location and each year folder to simulation
    simulation_folders = []
    for place in os.listdir(root_folder):
        for year in os.listdir(os.path.join(root_folder, place)):
            simulation_folders.append(os.path.join(root_folder, place, year))

    if len(simulation_folders) == 0:
        msg = "No <location>/<year> folders found. Did you set the correct root folder ({}) ?".format(root_folder)
        raise Exception(msg)

    if not args.silent:
        print("The following simulations will be carried out:\n" + '\n'.join(simulation_folders))

    coords = pd.DataFrame(index=["latitude", "longitude"])

    for simulation_folder in simulation_folders:
        if not args.silent:
            print("Simulating {} ...".format(simulation_folder))

        # Use the location and the year as the prefix
        output_file_prefix = os.path.join(args.output, '-'.join(simulation_folder.split('/')[-2:]))

        # Catch the meta information
        meta = run_pv_simulation_with_surfrad(output_file_prefix, scenarios, simulation_folder, not args.silent)

        # Record the location
        coords[meta["name"]] = [meta["latitude"], meta["longitude"]]

    # Write coordinates to CSV
    coords.to_csv(os.path.join(args.output, 'coordinates.csv'))

    if not args.silent:
        print("Simulation with SURFRAD data is complete!")

