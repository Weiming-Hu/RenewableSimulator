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

import numpy as np

from time import time
from netCDF4 import Dataset
from Scenarios import Scenarios
from Functions import simulate_sun_positions, read_yaml


###################
# Class Simulator #
###################

class Simulator:
    def __init__(self, scenarios, simple_clock=False, stations_index=None, verbose=True):

        if isinstance(scenarios, str):
            self.scenarios = Scenarios(read_yaml(scenarios))
        else:
            assert isinstance(scenarios, Scenarios), 'Only accept Scenarios or yaml file path'
            self.scenarios = scenarios

        self.stations_index = stations_index
        self.verbose = verbose
        self.simple_clock = simple_clock
        self.simple_clock = {'timestamps': [time()], 'log_names': []}

        if self.stations_index is not None:
            assert isinstance(stations_index, list), 'Station indices should be a list'
            assert len(stations_index) == len(set(stations_index)), 'Duplicates found in station indices'

            # Make sure the station indices are sorted
            self.stations_index = stations_index
            self.stations_index.sort()

    def simulate(self):
        raise NotImplementedError

    def write(self):
        raise NotImplementedError

    def summary(self):
        raise NotImplementedError

    def _log_event(self, log_name):
        if self.simple_clock:
            self.simple_clock['timestamps'].append(time())
            self.simple_clock['log_names'].append(log_name)


###############################
# Class SimulatorSolarAnalogs #
###############################

class SimulatorSolarAnalogs(Simulator):

    def __init__(self, nc_file, variable_dict, scenarios, solar_position_method='nrel_numpy',
                 simple_clock=False, parallel_nc=False, stations_index=None,
                 cores=1, verbose=True, disable_progress_bar=False):

        super().__init__(scenarios, simple_clock, stations_index, verbose)

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
            'lon': None,
            'lat': None,
            'test_times': None,
            'lead_times': None,
            'analogs': {'ghi': None, 'alb': None, 'wspd': None, 'tamb': None},
            'fcsts': {'ghi': None, 'alb': None, 'wspd': None, 'tamb': None},
            'obs': {'ghi': None, 'alb': None, 'wspd': None, 'tamb': None}
        }

        # Read data from NetCDF generated from Analog Ensemble
        self._read_simulation_data()

        # Pre-calculate air mass and extraterrestrial irradiance from solar positions
        self.simulation_data['sky'] = simulate_sun_positions(
            days=self.simulation_data['test_times'],
            lead_times=self.simulation_data['lead_times'],
            lats=self.simulation_data['lats'],
            lons=self.simulation_data['lons'],
            solar_position_method=self.solar_position_method,
            silent=not self.verbose)

        self._log_event('Calculate sun positions')

    def simulate(self):
        pass

    def write(self):
        pass

    def summary(self):
        summary_info = 'Simulation summary:\n' + \
                       '-- {} scenarios\n'.format(self.scenarios.total_scenarios()) + \
                       '-- {} stations\n'.format(len(self.simulation_data['lon'])) + \
                       '-- {} test times\n'.format(len(self.simulation_data['test_times'])) + \
                       '-- {} lead times\n'.format(len(self.simulation_data['lead_times'])) + \
                       '-- {} analog members\n'.format(self.simulation_data['analogs']['ghi'].shape[0]) + \
                       '-- {} processes to be created'.format(self.cores)

        return summary_info

    def _read_simulation_data(self):

        # Open connection
        if self.verbose:
            print('Open {} connection to {}'.format('parallel' if self.parallel_nc else 'sequential', self.nc_file))

        nc = Dataset(self.nc_file, 'r', parallel=self.parallel_nc)

        # Determine station indices
        if self.stations_index is None:
            num_stations = nc.dimensions['num_stations'].size
            self.stations_index = list(range(num_stations))

        ####################
        # Read analog data #
        ####################

        for value_key in self.simulation_data['analogs'].keys():

            # Get the variable object
            variable = nc.variables[self.variable_dict[value_key]]

            # Allow parallel access if parallel access
            if self.parallel_nc:
                variable.set_collective(True)

            # Dimensions are [members, lead times, test times, stations]
            self.simulation_data['analogs'][value_key] = variable[:, :, :, self.stations_index]

        self._log_event('Read analogs')

        ######################################
        # Read forecast and observation data #
        ######################################

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
                    self.simulation_data[type_key][value_key] = np.transpose(value, (2, 0, 1))

        self._log_event('Read forecasts and observations')

        ##################
        # Read Meta data #
        ##################

        self.simulation_data['lon'] = nc.variables[self.variable_dict['lon']][self.stations_index]
        self.simulation_data['lat'] = nc.variables[self.variable_dict['lat']][self.stations_index]
        self.simulation_data['test_times'] = nc.variables[self.variable_dict['test_times']][:]
        self.simulation_data['lead_times'] = nc.variables[self.variable_dict['lead_times']][:]

        ###################
        # Unit conversion #
        ###################

        for type_key in ['analogs', 'fcsts', 'obs']:

            # Albedo from percentage to decimal
            self.simulation_data[type_key]['alb'] /= 100.0

            # Temperature from Kelvin to Celsius
            self.simulation_data[type_key]['tamb'] -= 273.15

        nc.close()

        self._log_event('Preprocess simulation data')
