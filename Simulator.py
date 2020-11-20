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
import gc

from glob import glob
from Functions import *
from Timer import Timer
from pvlib import iotools
from netCDF4 import Dataset
from Scenarios import Scenarios


###################
# Class Simulator #
###################

class Simulator:
    def __init__(self, scenarios, verbose=True):

        if isinstance(scenarios, str):
            self.scenarios = Scenarios(read_yaml(os.path.expanduser(scenarios)))
        else:
            assert isinstance(scenarios, Scenarios), 'Only accept Scenarios or yaml file path'
            self.scenarios = scenarios

        self.verbose = verbose
        self.timer = Timer()

    def simulate(self):
        raise NotImplementedError

    def summary(self):
        raise NotImplementedError


########################
# Class SimulatorSolar #
########################

class SimulatorSolar(Simulator):
    def __init__(self, nc_file, scenarios, solar_position_method='nrel_numpy',
                 parallel_nc=False, stations_index=None, read_sky_conditions=True,
                 cores=1, verbose=True, disable_progress_bar=False):

        super().__init__(scenarios, verbose)

        self.cores = cores
        self.parallel_nc = parallel_nc
        self.stations_index = stations_index
        self.nc_file = os.path.expanduser(nc_file)
        self.read_sky_conditions = read_sky_conditions
        self.disable_progress_bar = disable_progress_bar
        self.solar_position_method = solar_position_method

        if self.stations_index is not None:
            assert isinstance(stations_index, list), 'Station indices should be a list'
            assert len(stations_index) == len(set(stations_index)), 'Duplicates found in station indices'

            # Make sure the station indices are sorted
            self.stations_index = stations_index
            self.stations_index.sort()

        # A predefined dictionary format
        self.simulation_data = {'longitudes': None, 'latitudes': None, 'test_times': None, 'lead_times': None}

        # This function populate the dictionary with simulation data
        self._read_simulation_data()
        gc.collect()

        # Prepare air mass and extraterrestrial irradiance from solar positions
        self._prepare_sky_conditions()
        gc.collect()

    def simulate(self):
        raise NotImplementedError

    def summary(self):
        summary_info = '*************** Summary ***************\n' + \
                       'Overview:\n' + \
                       '-- {} scenarios\n'.format(self.scenarios.total_scenarios()) + \
                       '-- {} stations\n'.format(len(self.simulation_data['longitudes'])) + \
                       '-- {} test times\n'.format(len(self.simulation_data['test_times'])) + \
                       '-- {} lead times\n'.format(len(self.simulation_data['lead_times'])) + \
                       '-- {} cores\n'.format(self.cores)

        summary_info += '\nSimulation data summary:\n' + recursive_summary_dict(self.simulation_data)
        summary_info += '*********** End of Summary ************'

        return summary_info

    def _read_simulation_data(self):
        raise NotImplementedError

    def _prepare_sky_conditions(self):

        nc = Dataset(self.nc_file, 'a', parallel=self.parallel_nc)

        if 'SkyConditions' in nc.groups and self.read_sky_conditions:
            self.timer.start('Read sky conditions')

            if self.verbose:
                print('Reading sky conditions ...')

            sky = read_array_dict(nc, 'SkyConditions', self.parallel_nc, self.stations_index)

            # Sanity check
            required = ('dni_extra', 'air_mass', 'zenith', 'apparent_zenith', 'azimuth')
            assert all([k in sky.keys() for k in required]), 'Sky condition missing. Require {}'.format(required)

            expected_dims = (nc.dimensions['num_flts'].size,
                             nc.dimensions['num_test_times'].size,
                             len(self.stations_index))

            d_dims = [v.shape for v in sky.values()]
            assert all([d_dims[0] == dim for dim in d_dims]), 'All arrays should have the same shape in the dictionary'
            assert d_dims[0] == expected_dims, 'Expect dimensions {} got {}'.format(expected_dims, d_dims[0])

            self.timer.stop()

        else:
            self.timer.start('Calculate sky conditions')

            sky = simulate_sun_positions(
                days=self.simulation_data['test_times'],
                lead_times=self.simulation_data['lead_times'],
                latitudes=self.simulation_data['latitudes'],
                longitudes=self.simulation_data['longitudes'],
                solar_position_method=self.solar_position_method,
                disable_progress_bar=self.disable_progress_bar,
                cores=self.cores, verbose=self.verbose)

            self.timer.stop()
            self.timer.start('Write sky conditions')

            if self.verbose:
                print('Writing sky conditions ...')

            write_array_dict(nc, 'SkyConditions', sky, ('num_flts', 'num_test_times', 'num_stations'),
                             self.parallel_nc, self.stations_index)

            self.timer.stop()

        # Merge results
        assert all([k not in sky.keys() for k in self.simulation_data.keys()]), 'Duplicate names found during merging'
        self.simulation_data = {**self.simulation_data, **sky}


###############################
# Class SimulatorSolarSurfrad #
###############################

class SimulatorSolarSurfrad(SimulatorSolar):

    def __init__(self, data_folder, nc_file, scenarios, file_ext='.dat', solar_position_method='nrel_numpy',
                 cores=1, verbose=True, disable_progress_bar=False):

        if verbose:
            print('Initializing PV simulation with SURFRAD ...')

        self.data_folder = os.path.expanduser(data_folder)
        assert os.path.isdir(self.data_folder), '{} does not exist'.format(self.data_folder)

        # Extract available data files
        match = os.path.join(self.data_folder, "**/*{}".format(file_ext))
        self.data_files = [f for f in glob(match, recursive=True)]
        assert len(self.data_files) > 0, 'No files have been found with the match criteria {}'.format(match)

        super().__init__(nc_file, scenarios, solar_position_method, False, None,
                         False, cores, verbose, disable_progress_bar)

    def simulate(self):

        # Sanity check
        required = ['surfrad', 'test_times', 'lead_times', 'air_mass',
                    'dni_extra', 'zenith', 'apparent_zenith', 'azimuth']
        assert all([key in self.simulation_data.keys() for key in required]), 'Not properly initialized'

        # Open Connection
        nc = Dataset(self.nc_file, 'a', parallel=self.parallel_nc)

        simulate_power('surfrad', self.scenarios, nc,
                       self.simulation_data['surfrad']['ghi'], self.simulation_data['surfrad']['tamb'],
                       self.simulation_data['surfrad']['wspd'], self.simulation_data['surfrad']['alb'],
                       self.simulation_data['test_times'], self.simulation_data['lead_times'],
                       self.simulation_data['air_mass'], self.simulation_data['dni_extra'],
                       self.simulation_data['zenith'], self.simulation_data['apparent_zenith'],
                       self.simulation_data['azimuth'], self.parallel_nc, self.cores, self.verbose,
                       self.stations_index, self.disable_progress_bar, self.timer)

        nc.close()
        gc.collect()

        if self.verbose:
            print("Power simulation is complete!")

    def _read_simulation_data(self):

        #############
        # Read data #
        #############

        self.timer.start('Read files')

        if self.verbose:
            print('Reading {} SURFRAD files ...'.format(len(self.data_files)))

        data = process_map(iotools.read_surfrad, self.data_files, max_workers=self.cores,
                           disable=self.disable_progress_bar,
                           chunksize=1 if len(self.data_files) < 1000 else int(len(self.data_files) / 100))

        self.timer.stop()

        ########################
        # Lists to data frames #
        ########################

        self.timer.start('Lists to data frames')

        if self.verbose:
            print('Formatting SURFRAD data ...')

        # Add stations to each separate data frame
        for value in data:
            value[0]['_station_'] = value[1]['name']

        data = {
            'meta': pd.DataFrame([list(v[1].values()) for v in data], columns=data[0][1].keys()),
            'value': pd.concat([v[0] for v in data])
        }

        # Remove duplicated rows
        data['meta'] = data['meta'].drop_duplicates()

        self.timer.stop()

        #########################
        # Data frames to arrays #
        #########################

        self.timer.start('Data frames to arrays')

        # Calculate unix time
        posix_start = pd.Timestamp('1970/1/1', tz='UTC')
        times_delta = data['value'].index - posix_start
        unix_times = times_delta.total_seconds().to_numpy()
        data['value'].index = unix_times
        unix_times = np.unique(unix_times)

        self.simulation_data['stations'] = data['meta'].name.to_numpy()
        self.simulation_data['longitudes'] = data['meta'].longitude.to_numpy()
        self.simulation_data['latitudes'] = data['meta'].latitude.to_numpy()
        self.simulation_data['test_times'] = unix_times
        self.simulation_data['lead_times'] = np.array([0])

        with np.errstate(invalid='ignore'):

            # Save the shared data to disk
            if self.verbose:
                print('Preparing data for parallel processing ...')

            tmp_file = save_dict({'data': data, 'simulation_data': self.simulation_data})

            if self.verbose:
                print('Temporary file saved to {}'.format(tmp_file))

            # Define a simple wrapper
            def wrapper(station_index, tmp, length):
                locals().update(read_dict(tmp))
                return SimulatorSolarSurfrad._align_data(station_index, data, simulation_data, length)

            wrapper = partial(wrapper, tmp=tmp_file, length=len(self.simulation_data['test_times']))

            surfrad = process_map(wrapper, range(len(self.simulation_data['stations'])), max_workers=self.cores,
                                  disable=self.disable_progress_bar,
                                  chunksize=1 if len(self.data_files) < 1000 else int(len(self.data_files) / 100))

            # Remove the temporary file
            os.remove(tmp_file)

            self.simulation_data['surfrad'] = {
                'ghi': np.expand_dims(np.expand_dims(np.stack([v[0] for v in surfrad], axis=1), 0), 0),
                'uw_solar': np.expand_dims(np.expand_dims(np.stack([v[1] for v in surfrad], axis=1), 0), 0),
                'alb': np.expand_dims(np.expand_dims(np.stack([v[2] for v in surfrad], axis=1), 0), 0),
                'wspd': np.expand_dims(np.expand_dims(np.stack([v[3] for v in surfrad], axis=1), 0), 0),
                'tamb': np.expand_dims(np.expand_dims(np.stack([v[4] for v in surfrad], axis=1), 0), 0),
            }

        self.timer.stop()

        ##################################
        # Save converted simulation data #
        ##################################

        self.timer.start('Save converted simulation data')

        if self.verbose:
            print('Writing simulation data to NetCDF ...')

        nc = Dataset(self.nc_file, 'w')

        # Create dimensions
        nc.createDimension('num_stations', len(self.simulation_data['stations']))
        nc.createDimension('num_test_times', len(self.simulation_data['test_times']))
        nc.createDimension('num_flts', 1)
        nc.createDimension('num_analogs', 1)

        # Create variables and copy values
        array_dims = ('num_analogs', 'num_flts', 'num_test_times', 'num_stations')

        var = nc.createVariable('stations', np.str, ('num_stations'))
        var[:] = self.simulation_data['stations']

        var = nc.createVariable('longitudes', np.float, ('num_stations'))
        var[:] = self.simulation_data['longitudes']

        var = nc.createVariable('latitudes', np.float, ('num_stations'))
        var[:] = self.simulation_data['latitudes']

        var = nc.createVariable('test_times', np.int, ('num_test_times'))
        var[:] = self.simulation_data['test_times']

        var = nc.createVariable('lead_times', np.int, ('num_flts'))
        var[:] = self.simulation_data['lead_times']

        write_array_dict(nc, 'surfrad', self.simulation_data['surfrad'], dimensions=array_dims, parallel_nc=False)

        nc.close()
        self.timer.stop()

    @staticmethod
    def _align_data(station_index, data, simulation_data, length):

        # Extract data for this station
        values = data['value'][data['value']._station_ == simulation_data['stations'][station_index]]

        # Figure out what times are available
        indices = [np.where(simulation_data['test_times'] == index)[0][0] for index in values.index]

        # Preprocessing
        ghi = values.ghi.to_numpy()
        uw_solar = values.uw_solar.to_numpy()

        ghi[ghi < 0] = 0
        mask = (ghi > 0) & (ghi > uw_solar)

        alb = np.zeros(len(ghi))
        alb[mask] = uw_solar[mask] / ghi[mask]

        ghi_arr = np.full(length, np.nan)
        uw_solar_arr = np.full(length, np.nan)
        alb_arr = np.full(length, np.nan)
        wspd_arr = np.full(length, np.nan)
        tamb_arr = np.full(length, np.nan)

        ghi_arr[indices] = ghi
        uw_solar_arr[indices] = uw_solar
        alb_arr[indices] = alb
        wspd_arr[indices] = values.wind_speed
        tamb_arr[indices] = values.temp_air

        return [ghi_arr, uw_solar_arr, alb_arr, wspd_arr, tamb_arr]


###############################
# Class SimulatorSolarAnalogs #
###############################

class SimulatorSolarAnalogs(SimulatorSolar):

    def __init__(self, nc_file, variable_dict, scenarios, solar_position_method='nrel_numpy',
                 parallel_nc=False, stations_index=None, read_sky_conditions=True,
                 cores=1, verbose=True, disable_progress_bar=False):

        if verbose:
            print('Initializing PV simulation with Analog Ensemble ...')

        # Sanity checks
        assert os.path.isfile(nc_file), '{} does not exist'.format(nc_file)

        # Initialization
        self.variable_dict = variable_dict

        # Process variable dict
        if isinstance(self.variable_dict, str):
            self.variable_dict = read_yaml(self.variable_dict)

        super().__init__(nc_file, scenarios, solar_position_method, parallel_nc, stations_index,
                         read_sky_conditions, cores, verbose, disable_progress_bar)

    def simulate(self):

        # Sanity check
        required = ['analogs', 'fcsts', 'obs', 'test_times', 'lead_times',
                    'air_mass', 'dni_extra', 'zenith', 'apparent_zenith', 'azimuth']
        assert all([key in self.simulation_data.keys() for key in required]), 'Not properly initialized'

        # Open Connection
        self.timer.start('Open a read connection')
        nc = Dataset(self.nc_file, 'a', parallel=self.parallel_nc)
        self.timer.stop()

        for name in ['analogs', 'fcsts', 'obs']:

            simulate_power(name, self.scenarios, nc,
                           self.simulation_data[name]['ghi'], self.simulation_data[name]['tamb'],
                           self.simulation_data[name]['wspd'], self.simulation_data[name]['alb'],
                           self.simulation_data['test_times'], self.simulation_data['lead_times'],
                           self.simulation_data['air_mass'], self.simulation_data['dni_extra'],
                           self.simulation_data['zenith'], self.simulation_data['apparent_zenith'],
                           self.simulation_data['azimuth'], self.parallel_nc, self.cores, self.verbose,
                           self.stations_index, self.disable_progress_bar, self.timer)

        nc.close()
        gc.collect()

        if self.verbose:
            print("Power simulation is complete!")

    def _read_simulation_data(self):

        # Additional data that will be populated during reading
        self.simulation_data['analogs'] = {'ghi': None, 'alb': None, 'wspd': None, 'tamb': None}
        self.simulation_data['fcsts'] = {'ghi': None, 'alb': None, 'wspd': None, 'tamb': None}
        self.simulation_data['obs'] = {'ghi': None, 'alb': None, 'wspd': None, 'tamb': None}

        # Open connection
        self.timer.start('Open a read connection')

        if self.verbose:
            print('Opening a {} connection to {} ...'.format(
                'parallel' if self.parallel_nc else 'sequential', self.nc_file))

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
            print('Aligning observations ...')

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
