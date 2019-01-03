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

import functions as funcs
from sscapi import PySSC

def run_pvwattsv5():

    # Collect input data
    template_file = 'weatherTemplate.csv'
    ssc = PySSC()
    scc_data = ssc.data_create()
    ssc.data_set_string(scc_data, 'file_name', template_file)
    ssc.module_exec_simple_no_thread('wfreader', scc_data)

    input = {}
    input['lat'] = ssc.data_get_number(scc_data, 'lat')
    input['lon'] = ssc.data_get_number(scc_data, 'lon')
    input['tz'] = ssc.data_get_number(scc_data, 'tz')
    input['elev'] = ssc.data_get_number(scc_data, 'elev')
    input['year'] = ssc.data_get_array(scc_data, 'year')
    input['month'] = ssc.data_get_array(scc_data, 'month')
    input['day'] = ssc.data_get_array(scc_data, 'day')
    input['hour'] = ssc.data_get_array(scc_data, 'hour')
    input['beam'] = ssc.data_get_array(scc_data, 'beam')
    input['diffuse'] = ssc.data_get_array(scc_data, 'diffuse')
    input['wspd'] = ssc.data_get_array(scc_data, 'wspd')
    input['tdry'] = ssc.data_get_array(scc_data, 'tdry')
    input['albedo'] = ssc.data_get_array(scc_data, 'albedo')

    ssc.data_free(scc_data)

    # Run simulation
    output = funcs.pvwattsv5(input)

    return output


if __name__ == '__main__':

    data_template = funcs.readWeatherTemplate(['lat', 'lon'], ['wspd'])

    funcs.printModuleVariables("pvwattsv5")
    funcs.printModuleVariables("pvwattsv5_1ts")

    funcs.printAvailableModules()

    nc_file = "/home/graduate/wuh20/github/RenewableSimulator/analogs.nc"
    i_stations = [1, 2, 3, 4]
    i_times = [0]
    i_flts = [0, 1, 2, 3, 4, 5]
    i_members = [0, 1, 2]
    var_names = ["Analogs", "Times", "Xs", "Ys"]
    data = funcs.getDataArray(nc_file, var_names, i_stations, i_times, i_flts, i_members)
