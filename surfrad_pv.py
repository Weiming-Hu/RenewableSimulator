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
from datetime import datetime

# Scientific python add-ons
import yaml
import numpy as np
import pandas as pd
from pvlib import temperature, pvsystem

# Self hosted modules
from Scenarios import Scenarios
from Functions import read_hourly_surfrad, simulate_power, simulate_sun_positions


def run_pv_simulation_with_surfrad(output_file, scenarios, year_folder, progress=True):

    # Read hourly SURFRAD data
    yearly_data, meta = read_hourly_surfrad(year_folder, progress)

    # Convert value table to numpy array to use the same the workflow as evergreen
    ghi_arr = np.array(yearly_data.ghi).reshape((1, 1, yearly_data.shape[0], 1))
    albedo_arr = np.array(yearly_data.uw_solar/yearly_data.ghi).reshape((1, 1, yearly_data.shape[0], 1))
    wspd_arr = np.array(yearly_data.wind_speed).reshape((1, 1, yearly_data.shape[0], 1))
    tamb_arr = np.array(yearly_data.temp_air).reshape((1, 1, yearly_data.shape[0], 1))

    # Confine the calculated albedo
    albedo_arr[:, :, yearly_data.ghi <= 0, :] = 0
    albedo_arr[albedo_arr > 1] = 0

    # Append the calculated albedo to the original data frame
    yearly_data['calculated_albedo'] = albedo_arr[0, 0, :, 0]

    # Calculate lead times from observation times
    posix_flt = [0]
    posix_start = datetime(1970, 1, 1)
    datetime_strs = ['{}/{}/{} {}:{}'.format(
        int(row['year']), int(row['month']), int(row['day']), int(row['hour']), int(row['minute']))
        for _, row in yearly_data.iterrows()]
    posix_days = [(datetime.strptime(datetime_str, '%Y/%m/%d %H:%M') - posix_start).total_seconds()
                  for datetime_str in datetime_strs]

    # Pre-calculate air mass and extraterrestrial irradiance from solar positions
    sky_dict = simulate_sun_positions(posix_days, posix_flt, [meta['latitude']], [meta['longitude']], 'nrel_numba', 0)
    yearly_data['dni_extra'] = sky_dict["dni_extra"][0, :, 0]
    yearly_data['air_mass'] = sky_dict["air_mass"][0, :, 0]
    yearly_data['zenith'] = sky_dict["zenith"][0, :, 0]
    yearly_data['azimuth'] = sky_dict["azimuth"][0, :, 0]

    # Determine the total number of scenarios
    num_scenarios = scenarios.total_scenarios()

    for scenario_index in range(num_scenarios):

        if progress:
            print("Simulating scenario {}/{}".format(scenario_index, num_scenarios))

        current_scenario = scenarios.get_scenario(scenario_index)

        # get the current scenario settings
        surface_tilt = current_scenario["surface_tilt"]
        surface_azimuth = current_scenario["surface_azimuth"]
        pv_module = pvsystem.retrieve_sam("SandiaMod")[current_scenario["pv_module"]]
        tcell_model_parameters = \
            temperature.TEMPERATURE_MODEL_PARAMETERS["sapm"][current_scenario["tcell_model_parameters"]]

        p_mp = simulate_power(ghi_arr, tamb_arr, wspd_arr, albedo_arr, posix_days, posix_flt, sky_dict,
                              surface_tilt, surface_azimuth, pv_module, tcell_model_parameters, 0)

        # Append the results to a data frame
        column_name = "power_scenario_{:05d}".format(scenario_index)
        yearly_data[column_name] = p_mp[0, 0, :, 0]

    # Write the current data frame to a CSV file
    if progress:
        print("Writing to {}".format(output_file))
    yearly_data.to_csv(output_file)

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
    
    args.output = os.path.expanduser(args.output)
    args.root = os.path.expanduser(args.root)

    # Check output folder
    if os.path.exists(args.output):
        if os.path.isdir(args.output):
            print("Output folder {} exists. Contents will be overwritten!".format(args.output))
        else:
            raise Exception('{} is not a folder'.format(args.output))
    else:
        os.makedirs(args.output)

    # Read YAML files
    with open(os.path.expanduser(args.scenario)) as f:
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
        output_file = os.path.join(args.output, '-'.join(simulation_folder.split('/')[-2:]) + '.csv')

        # Catch the meta information
        meta = run_pv_simulation_with_surfrad(output_file, scenarios, simulation_folder, not args.silent)

        # Record the location
        coords[meta["name"]] = [meta["latitude"], meta["longitude"]]

    # Write coordinates to CSV
    coords.to_csv(os.path.join(args.output, 'coordinates.csv'))

    # Write scenarios to YAML
    scenarios.write_all(os.path.join(args.output, 'scenarios.yaml'))

    if not args.silent:
        print("Simulation with SURFRAD data is complete!")

