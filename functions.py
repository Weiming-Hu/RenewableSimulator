# "`-''-/").___..--''"`-._
#  (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
#  (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#    _ ..`--'_..-_/  /--'_.' ,'
#  (il),-''  (li),'  ((!.-'
#
# Author: Weiming Hu <weiming@psu.edu>
#         Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#         Department of Geography and Institute for CyberScience
#         The Pennsylvania State University
#

# This file contains the pre-/post- processing functions for SAM model simulation

try:
    from netCDF4 import Dataset
except ImportError:
    print "Please install the netCDF4 module."
    quit(1)

from sscapi import PySSC


# Create mapping information from numeric to string
ssc_data_types = ("Invalid", "String", "Number", "Array", "Matrix", "Table")
ssc_var_types = ("Invalid", "Input", "Output", "Both")



def pvwattsv5(input):

    # Create SSC container and data repository
    ssc = PySSC()
    ssc_data = ssc.data_create()

    # Set data
    resource_data = ssc.data_create()
    ssc.data_set_number(resource_data, 'lat', input['lat'])
    ssc.data_set_number(resource_data, 'lon', input['lon'])
    ssc.data_set_number(resource_data, 'tz',  input['tz'])
    ssc.data_set_number(resource_data, 'elev',  input['elev'])

    ssc.data_set_array(resource_data, 'year',  input['year'])
    ssc.data_set_array(resource_data, 'month',  input['month'])
    ssc.data_set_array(resource_data, 'day',  input['day'])
    ssc.data_set_array(resource_data, 'hour', input['hour'])

    ssc.data_set_array(resource_data, 'dn', input['beam'])
    ssc.data_set_array(resource_data, 'df', input['diffuse'])
    ssc.data_set_array(resource_data, 'wspd', input['wspd'])
    ssc.data_set_array(resource_data, 'tdry', input['tdry'])
    ssc.data_set_array(resource_data, 'albedo', input['albedo'])

    ssc.data_set_table(ssc_data, 'solar_resource_data', resource_data )
    ssc.data_free(resource_data)

    # Create simulation module
    mod = ssc.module_create("pvwattsv5")
    ssc.module_exec_set_print(0)

    ssc.data_set_number(ssc_data, 'system_capacity', 4)
    ssc.data_set_number(ssc_data, 'module_type', 0)
    ssc.data_set_number(ssc_data, 'array_type', 0)
    ssc.data_set_number(ssc_data, 'losses', 14)
    ssc.data_set_number(ssc_data, 'tilt', 15)
    ssc.data_set_number(ssc_data, 'azimuth', 180)
    ssc.data_set_number(ssc_data, 'adjust:constant', 0)

    # Run simulation
    if ssc.module_exec(mod, ssc_data) == 0:
        print 'PVWatts V5 simulation error'
        idx = 1
        msg = ssc.module_log(mod, 0)
        while (msg != None):
            print '\t: ' + msg
            msg = ssc.module_log(mod, idx)
            idx = idx + 1

    # Get the output variable name list
    out_vars_index = _getVariablesIndexByType(ssc, mod)

    # Extract output variables
    ret_dict = {}
    for index in out_vars_index:
        p_info = ssc.module_var_info(mod, index)
        name = ssc.info_name(p_info)
        data_type = ssc.info_data_type(p_info)

        if data_type == ssc_data_types.index("Invalid"):
            print "Skip output variable {}.".format(name)
        elif data_type == ssc_data_types.index("String"):
            ret_dict[name] = ssc.data_get_string(ssc_data, name)
        elif data_type == ssc_data_types.index("Number"):
            ret_dict[name] = ssc.data_get_number(ssc_data, name)
        elif data_type == ssc_data_types.index("Array"):
            ret_dict[name] = ssc.data_get_array(ssc_data, name)
        else:
            print "Skip output variable {} because the type is not supported.".format(name)

    ssc.module_free(mod)
    ssc.data_free(ssc_data)
    return ret_dict



def readWeatherTemplate(number_names, array_names, weather_file = "weatherTemplate.csv"):
    # This function reads the specified variables from a weather file. The type of variables
    # is either a number or an array and they are all read and put into a dictionary, and returned.
    # Unrecognized names will be ignored and a warning message will be printed.
    #

    # Create data repository
    ssc = PySSC()
    ssc_data = ssc.data_create()

    # Build connection to the weather file for reading
    ssc.data_set_string(ssc_data, 'file_name', weather_file)
    ssc.module_exec_simple_no_thread('wfreader', ssc_data)
    if ssc.data_first(ssc_data) == 'file_name':
        print "Something was wrong when reading the weather file {}.".format(weather_file)
        return None

    # Read specified variables from the template file into a dictionary
    data_dict = {}
    for name in array_names:
        if ssc.data_query(ssc_data, name) == ssc_data_types.index("Array"):
            data_dict[name] = ssc.data_get_array(ssc_data, name)
        else:
            print "Cannot find the array {}.".format(name)
    for name in number_names:
        if ssc.data_query(ssc_data, name) == ssc_data_types.index("Number"):
            data_dict[name] = ssc.data_get_number(ssc_data, name)
        else:
            print "Cannot find the number {}.".format(name)

    ssc.data_free(ssc_data)
    return data_dict



def printModuleVariables(mod_name):
    # Print out variable information for a module

    # Create container and module
    ssc = PySSC()
    mod = ssc.module_create(mod_name)

    # Check whether module has been successfully created
    if mod is None:
        print "Unrecognized module name ({}).".format(mod_name)
        return None

    # Print information
    i = 0
    print "VariableType\tGroup\tDataType\tName\tUnit\tLabel\tMeta"
    print "------------------------------------------------------"
    while True:
        p_info = ssc.module_var_info(mod, i)

        if p_info is None:
            print "* End of variables *"
            break

        print "{},\t{},\t{},\t{},\t{},\t{},\t{}".format(
            ssc_var_types[ssc.info_var_type(p_info)],
            ssc.info_group(p_info),
            ssc_data_types[ssc.info_data_type(p_info)],
            ssc.info_name(p_info),
            ssc.info_units(p_info),
            ssc.info_label(p_info),
            ssc.info_meta(p_info))
        i += 1

    ssc.module_free(mod)



def printAvailableModules():
    # Print all supported models from the SAM model.

    # Create a container
    ssc = PySSC()

    # Walk through all available modules
    i = 0
    while True:
        entry = ssc.module_entry(i)
        if entry is None:
            print "* End of modules *"
            break
        print "{}: {}".format(ssc.entry_name(entry), ssc.entry_description(entry))
        i += 1



def getNcdfArray(nc_file, var_names, i_stations = [], i_times = [], i_flts = [], i_members = [], i_col = 0):
    # This function reads a certain part from a NetCDF variable and returns a dict of array.
    # In cases of error, it returns None.
    #

    # Open an NetCDF file
    nc = Dataset(nc_file)

    # Check the NetCDF file
    for var_name in var_names:
        if var_name not in nc.variables:
            print "The NetCDF file {} does not have {}!".format(nc_file, var_name)
            return None

    for name in ("num_test_stations", "num_test_times", "num_flts"):
        if name not in nc.dimensions:
            print "The NetCDF file {} does not have the dimension {}!".format(nc_file, name)
            return None

    if i_col >= 3:
        print "i_col ({}) should be smaller than 3!".format(i_col)
        return None

    # Read data
    ret_dict = {}

    for var_name in var_names:
        arr = None

        if var_name == "Analogs":
            arr = nc.variables[var_name][i_col, i_members, i_flts, i_times, i_stations]

            # Swap dimensions
            arr = arr.swapaxes(0, 3)
            arr = arr.swapaxes(1, 2)

        elif var_name in ("Xs", "Ys"):
            arr = nc.variables[var_name][i_stations]

        elif var_name in ("Times"):
            arr = nc.variables[var_name][i_times]

        elif var_name in ("FLTs"):
            arr = nc.variables[var_name][i_flts]

        ret_dict[var_name] = arr

    # Close the NetCDF file
    nc.close()

    return ret_dict



def _getVariablesIndexByType(ssc, module, type = "Output"):
    # Get the variable names by the variable type

    ret = []
    i = 0
    while True:
        p_info = ssc.module_var_info(module, i)

        if p_info is None:
            break

        if ssc.info_var_type(p_info) == ssc_var_types.index(type):
            ret.append(i)
        i += 1

    return ret
