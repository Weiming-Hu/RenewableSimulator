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
import gc

import numpy as np

from Timer import Timer
from netCDF4 import Dataset
from Scenarios import Scenarios
from Functions import simulate_sun_positions, simulate_power, read_yaml


###################
# Class Simulator #
###################

class Simulator:
    def __init__(self, scenarios, stations_index=None, verbose=True):

        if isinstance(scenarios, str):
            self.scenarios = Scenarios(read_yaml(scenarios))
        else:
            assert isinstance(scenarios, Scenarios), 'Only accept Scenarios or yaml file path'
            self.scenarios = scenarios

        self.stations_index = stations_index
        self.verbose = verbose
        self.timer = Timer()

        if self.stations_index is not None:
            assert isinstance(stations_index, list), 'Station indices should be a list'
            assert len(stations_index) == len(set(stations_index)), 'Duplicates found in station indices'

            # Make sure the station indices are sorted
            self.stations_index = stations_index
            self.stations_index.sort()

    def simulate(self):
        raise NotImplementedError

    def summary(self):
        raise NotImplementedError


###############################
# Class SimulatorSolarAnalogs #
###############################

class SimulatorSolarAnalogs(Simulator):

    def __init__(self, nc_file, variable_dict, scenarios, solar_position_method='nrel_numpy',
                 parallel_nc=False, stations_index=None,
                 cores=1, verbose=True, disable_progress_bar=False):

        super().__init__(scenarios, stations_index, verbose)

        if self.verbose:
            print('Initializing PV simulation with Analog Ensemble ...')

        # Sanity checks
        assert os.path.isfile(nc_file), '{} does not exist'.format(nc_file)

        # Initialization
        self.cores = cores
        self.parallel_nc = parallel_nc
        self.variable_dict = variable_dict
        self.nc_file = os.path.expanduser(nc_file)
        self.solar_position_method = solar_position_method
        self.disable_progress_bar = disable_progress_bar

        # Process variable dict
        if isinstance(self.variable_dict, str):
            self.variable_dict = read_yaml(self.variable_dict)

        self.simulation_data = {
            'longitudes': None,
            'latitudes': None,
            'test_times': None,
            'lead_times': None,
            'analogs': {'ghi': None, 'alb': None, 'wspd': None, 'tamb': None},
            'fcsts': {'ghi': None, 'alb': None, 'wspd': None, 'tamb': None},
            'obs': {'ghi': None, 'alb': None, 'wspd': None, 'tamb': None}
        }

        # Read data from NetCDF generated from Analog Ensemble
        self._read_simulation_data()

        # Pre-calculate air mass and extraterrestrial irradiance from solar positions
        self.timer.start('Calculate sun positions')

        if self.verbose:
            print('Simulating sun positions ...')

        sky = simulate_sun_positions(
            days=self.simulation_data['test_times'],
            lead_times=self.simulation_data['lead_times'],
            latitudes=self.simulation_data['latitudes'],
            longitudes=self.simulation_data['longitudes'],
            solar_position_method=self.solar_position_method,
            disable_progress_bar=self.disable_progress_bar,
            cores=self.cores)

        # Merge results
        self.simulation_data = {**self.simulation_data, **sky}
        self.timer.stop()

    def simulate(self):

        # Sanity check
        required = ['analogs', 'fcsts', 'obs', 'test_times', 'lead_times',
                    'air_mass', 'dni_extra', 'zenith', 'apparent_zenith', 'azimuth']
        assert all([key in self.simulation_data.keys() for key in required]), 'Not properly initialized'

        # Open Connection
        self.timer.start('Open a read connection')
        nc = Dataset(self.nc_file, 'a', parallel=self.parallel_nc)
        self.timer.stop()

        for key, name in {'analogs': 'Power simulation with analogs',
                          'fcsts': 'Power simulation with forecasts',
                          'obs': 'Power simulation with observations'}.items():

            simulate_power(key, name, self.scenarios, nc,
                           self.simulation_data[key]['ghi'], self.simulation_data[key]['tamb'],
                           self.simulation_data[key]['wspd'], self.simulation_data[key]['alb'],
                           self.simulation_data['test_times'], self.simulation_data['lead_times'],
                           self.simulation_data['air_mass'], self.simulation_data['dni_extra'],
                           self.simulation_data['zenith'], self.simulation_data['apparent_zenith'],
                           self.simulation_data['azimuth'], self.parallel_nc, self.cores, self.verbose,
                           self.stations_index, self.disable_progress_bar, self.timer)

        nc.close()

        if self.verbose:
            print("Power simulation is complete!")

    def summary(self):
        summary_info = '*************** Summary ***************' + \
                       'Overview:\n' + \
                       '-- {} scenarios\n'.format(self.scenarios.total_scenarios()) + \
                       '-- {} stations\n'.format(len(self.simulation_data['longitudes'])) + \
                       '-- {} test times\n'.format(len(self.simulation_data['test_times'])) + \
                       '-- {} lead times\n'.format(len(self.simulation_data['lead_times'])) + \
                       '-- {} analog members\n'.format(self.simulation_data['analogs']['ghi'].shape[0]) + \
                       '-- {} cores\n'.format(self.cores)

        def recursive_summary(d, prefix='--'):
            msg = ''

            for key, value in d.items():
                if isinstance(value, dict):
                    msg += '{} {}:\n'.format(prefix, key)
                    msg += recursive_summary(value, prefix + ' --')
                else:
                    msg += '{} {}: shape {}\n'.format(prefix, key, value.shape)

            return msg

        summary_info += '\nSimulation data summary:\n' + recursive_summary(self.simulation_data)
        summary_info += '*********** End of Summary ************'

        return summary_info

    def _read_simulation_data(self):

        # Open connection
        self.timer.start('Open a read connection')

        if self.verbose:
            print('Opening a {} connection to {} ...'.format('parallel' if self.parallel_nc else 'sequential', self.nc_file))

        nc = Dataset(self.nc_file, 'r', parallel=self.parallel_nc)

        # Determine station indices
        if self.stations_index is None:
            num_stations = nc.dimensions['num_stations'].size
            self.stations_index = list(range(num_stations))

        self.timer.stop()

        ####################
        # Read analog data #
        ####################
        self.timer.start('Read analogs')

        if self.verbose:
            print('Reading analogs ...')

        for value_key in self.simulation_data['analogs'].keys():

            # Get the variable object
            variable = nc.variables[self.variable_dict[value_key]]

            # Allow parallel access if parallel access
            if self.parallel_nc:
                variable.set_collective(True)

            # Dimensions are [members, lead times, test times, stations]
            self.simulation_data['analogs'][value_key] = variable[:, :, :, self.stations_index]

        self.timer.stop()

        ######################################
        # Read forecast and observation data #
        ######################################
        self.timer.start('Read forecasts and observations')

        if self.verbose:
            print('Reading forecasts and observations ...')

        for type_key in ['fcsts', 'obs']:

            # Determine the group and variables
            nc_group = nc.groups['Forecasts' if type_key == 'fcsts' else 'Observations']
            nc_data = nc_group.variables['Data']
            nc_parameters = nc_group.variables['ParameterNames']

            if self.parallel_nc:
                nc_data.set_collective(True)
                nc_parameters.set_collective(True)

            nc_parameters = nc_parameters[:]

            for value_key in self.simulation_data[type_key].keys():
                value_index = np.where(nc_parameters == self.variable_dict['{}_{}'.format(value_key, type_key)])[0]
                assert len(value_index) == 1, 'Failed to find index for {} in for key {}'.format(value_key, type_key)

                if type_key == 'fcsts':
                    # Dimensions are [lead times, test times, stations, 1 (variable)]
                    value = nc_data[:, :, self.stations_index, value_index]

                    # Transpose dimensions to be [1 (member), lead times, test times, stations]
                    self.simulation_data[type_key][value_key] = np.transpose(value, (3, 0, 1, 2))

                else:
                    # Dimensions are [test times, stations, 1 (variable)]
                    value = nc_data[:, self.stations_index, value_index]

                    # Transpose dimensions to be [1 (member), test times, stations]
                    value = np.transpose(value, (2, 0, 1))

                    # Expand a dimension to be [1 (member), 1 (lead time) , test times, stations]
                    self.simulation_data[type_key][value_key] = np.expand_dims(value, axis=1)

        self.timer.stop()

        ##################
        # Read Meta data #
        ##################
        self.timer.start('Read meta')

        for key in ['longitudes', 'latitudes']:
            variable = nc.variables[self.variable_dict[key]]

            if self.parallel_nc:
                variable.set_collective(True)

            self.simulation_data[key] = variable[self.stations_index]

        for key in ['test_times', 'lead_times']:
            variable = nc.variables[self.variable_dict[key]]

            if self.parallel_nc:
                variable.set_collective(True)

            self.simulation_data[key] = variable[:]

        self.timer.stop()

        ###################
        # Unit conversion #
        ###################

        self.timer.start('Preprocess simulation data')

        # Change units
        if self.verbose:
            print('Changing units ...')

        for type_key in ['analogs', 'fcsts', 'obs']:

            # Albedo from percentage to decimal
            self.simulation_data[type_key]['alb'] /= 100.0

            # Temperature from Kelvin to Celsius
            self.simulation_data[type_key]['tamb'] -= 273.15

        # Initialize dimensions
        if self.verbose:
            print('Align observations ...')

        num_stations = len(self.simulation_data['longitudes'])
        num_lead_times = len(self.simulation_data['lead_times'])
        num_test_times = len(self.simulation_data['test_times'])

        # Initialize arrays
        obs_dict = {
            'ghi': np.full((1, num_lead_times, num_test_times, num_stations), np.nan),
            'alb': np.full((1, num_lead_times, num_test_times, num_stations), np.nan),
            'wspd': np.full((1, num_lead_times, num_test_times, num_stations), np.nan),
            'tamb': np.full((1, num_lead_times, num_test_times, num_stations), np.nan)
        }

        # Read observation times
        obs_times = nc.groups["Observations"].variables["Times"]

        if self.parallel_nc:
            obs_times.set_collective(True)

        obs_times = obs_times[:]

        # Reshape observations to expand dimensions of lead times and test times
        for lead_time_index in range(num_lead_times):
            for day_index in range(num_test_times):

                fcst_time = self.simulation_data['lead_times'][lead_time_index] + \
                            self.simulation_data['test_times'][day_index]
                obs_time_index, = np.where(obs_times == fcst_time)

                if len(obs_time_index) == 1:
                    for key, value in obs_dict.items():
                        value[0, lead_time_index, day_index] = self.simulation_data['obs'][key][0, 0, obs_time_index]

        nc.close()

        self.simulation_data['obs'] = obs_dict
        self.timer.stop()

        gc.collect()
