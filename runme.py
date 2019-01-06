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



def run_pvwattsv5_1ts():

    # Configure input data
    other_number_input = {}
    other_number_input["year"] = 2018
    other_number_input["month"] = 7
    other_number_input["day"] = 1
    other_number_input["hour"] = 19
    other_number_input["minute"] = 0
    other_number_input["lat"] = 40.799733
    other_number_input["lon"] = -77.8433542
    other_number_input["tz"] = -5
    other_number_input["beam"] = 866.981975
    other_number_input["diffuse"] = 0
    other_number_input["tamb"] = 37.470625
    other_number_input["wspd"] = 3.16
    other_number_input["alb"] = 0.2
    # other_number_input["time_step"] = 
    other_number_input["system_capacity"] = 4
    other_number_input["module_type"] = 0
    other_number_input["dc_ac_ratio"] = 0.9
    # other_number_input["inv_eff"] = 
    other_number_input["losses"] = 14
    other_number_input["array_type"] = 0
    other_number_input["tilt"] = 15
    other_number_input["azimuth"] = 239.082665
    other_number_input["gcr"] = 0.5
    other_number_input["tcell"] = 37.470625
    other_number_input["poa"] = 0

    output = funcs.run_model('pvwattsv5_1ts', other_number_input = other_number_input)
    return output



def run_pvwattsv5():

    template_file = 'weatherTemplate.csv'
    model_name = "pvwattsv5"

    # Collect input data
    ssc = PySSC()
    scc_data = ssc.data_create()
    ssc.data_set_string(scc_data, 'file_name', template_file)
    ssc.module_exec_simple_no_thread('wfreader', scc_data)

    resource_number_input = {}
    resource_number_input['lat'] = ssc.data_get_number(scc_data, 'lat')
    resource_number_input['lon'] = ssc.data_get_number(scc_data, 'lon')
    resource_number_input['tz'] = ssc.data_get_number(scc_data, 'tz')
    resource_number_input['elev'] = ssc.data_get_number(scc_data, 'elev')
    
    other_number_input = {}
    other_number_input['tilt'] = 15
    other_number_input['losses'] = 14
    other_number_input['azimuth'] = 180
    other_number_input['array_type'] = 0
    other_number_input['module_type'] = 0
    other_number_input['system_capacity'] = 4
    other_number_input['adjust:constant'] = 0

    resource_array_input = {}
    resource_array_input['year'] = ssc.data_get_array(scc_data, 'year')
    resource_array_input['month'] = ssc.data_get_array(scc_data, 'month')
    resource_array_input['day'] = ssc.data_get_array(scc_data, 'day')
    resource_array_input['hour'] = ssc.data_get_array(scc_data, 'hour')
    resource_array_input['dn'] = ssc.data_get_array(scc_data, 'beam')
    resource_array_input['df'] = ssc.data_get_array(scc_data, 'diffuse')
    resource_array_input['wspd'] = ssc.data_get_array(scc_data, 'wspd')
    resource_array_input['tdry'] = ssc.data_get_array(scc_data, 'tdry')
    resource_array_input['albedo'] = ssc.data_get_array(scc_data, 'albedo')

    ssc.data_free(scc_data)

    # Run simulation
    output = funcs.run_model(model_name, resource_number_input, resource_array_input, other_number_input)

    return output



if __name__ == '__main__':

    # output = run_pvwattsv5()
    # print "[pvwattsv5] Annual ennergy: {}".format(output['annual_energy'])

    output = run_pvwattsv5_1ts()
    print "[pvwattsv5_1ts] AC power: {}; DC power: {}".format(output["ac"], output["dc"])

    # data_template = funcs.readWeatherTemplate(['lat', 'lon'], ['wspd'])

    # funcs.printModuleVariables("pvwattsv5")
    # funcs.printModuleVariables("pvwattsv5_1ts")

    # funcs.printAvailableModules()

    # nc_file = "/home/graduate/wuh20/github/RenewableSimulator/analogs.nc"
    # i_stations = [1, 2, 3, 4]
    # i_times = [0]
    # i_flts = [0, 1, 2, 3, 4, 5]
    # i_members = [0, 1, 2]
    # var_names = ["Analogs", "Times", "Xs", "Ys"]
    # data = funcs.getDataArray(nc_file, var_names, i_stations, i_times, i_flts, i_members)
