/* 
 *
 * File:   evergreen.cpp
 * Author: Weiming Hu <weiming@psu.edu>
 *
 * Created on February, 2020, 11:10 AM
 */

/** @file */

#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <iomanip>
#include <map>

#include "Ncdf.h"
#include "sscapi.h"
#include "AnEnReadNcdf.h"
#include "Array4DPointer.h"

using namespace std;
using namespace netCDF;

static const int _NUM_HOURS = 24;
static const int _NUM_DAYS = 365;
static const int _NUM_TIMES = _NUM_DAYS * _NUM_HOURS;

using Scenarios = map<string, vector<ssc_number_t> >;
using Scenario = map<string, ssc_number_t>;

void
getScenario(Scenario & scenario, const Scenarios & scenarios, size_t index) {

    // Clear map
    scenario.clear();

    // Initialize values
    size_t accum_len = 1;
    size_t current_index;

    for (auto e : scenarios) {

        // Calculate the current index for this key based on the global index
        if (accum_len == 1) {
            accum_len = e.second.size();
            current_index = fmod(index, accum_len);
        } else {
            current_index = fmod(index / accum_len, e.second.size());
            accum_len *= e.second.size();
        }

        // Insert pair
        scenario.insert(make_pair(e.first, e.second[current_index]));
    }

    return;
}

size_t
countScenarios(const Scenarios & map) {
    size_t total = 1;
    for (auto e : map) total *= e.second.size();
    return total;
}

void
run_pvwattsv5(const string & file_path, string output_var = "ac",
        Verbose verbose = Verbose::Progress) {

    /*
     * Set up configuration schemes
     */
    
    if (verbose >= Verbose::Progress) cout << "Setting up simulation ..." << endl;

    // Define our batches of scenarios to simulate.
    // Scenarios will be created using permutation of all keys.
    //
    Scenarios config_scenarios;

    config_scenarios["losses"] = {0, 0.2};
    config_scenarios["azimuth"] = {0};
    config_scenarios["system_capacity"] = {0.8};
    config_scenarios["tilt"] = {10};

    // Define our fixed configuration
    Scenario config_fixed;

    config_fixed["array_type"] = 0;
    config_fixed["adjust:constant"] = 0;
    config_fixed["verbose"] = 0;

    // TODO: How to set location and time zones
    config_fixed["lat"] = 40;
    config_fixed["lon"] = -70;
    config_fixed["tz"] = -6;


    /*
     * Read analog input
     */
    Array4DPointer dn_arr, df_arr, t_arr, wspd_arr;

    // TODO: Change variable name
    AnEnReadNcdf anen_read(verbose);
    anen_read.readAnalogs(file_path, dn_arr, "2m_relative-humidity");
    anen_read.readAnalogs(file_path, df_arr, "2m_relative-humidity");
    anen_read.readAnalogs(file_path, t_arr, "2m_relative-humidity");
    anen_read.readAnalogs(file_path, wspd_arr, "2m_relative-humidity");

    if (dn_arr.shape()[1] != _NUM_DAYS) throw runtime_error("Analogs must have 365 days for test times");
    if (dn_arr.shape()[2] != _NUM_HOURS) throw runtime_error("Analogs must have 24 hours for forecast lead times");

    /*
     * Run simulation for all
     * - stations
     * - ensemble members
     * - scenarios
     */
    size_t num_stations = dn_arr.shape()[0];
    size_t num_analogs = dn_arr.shape()[3];
    size_t num_scenarios = countScenarios(config_scenarios);

    // Calculate how wide should I pad the string to make pretty variable names
    size_t max_width = to_string(num_scenarios).size();

    // Create indices for dimensions that are not subset
    vector<size_t> all_test_times(_NUM_DAYS), all_flts(_NUM_HOURS);
    iota(all_test_times.begin(), all_test_times.end(), 0);
    iota(all_flts.begin(), all_flts.end(), 0);

    // Initialize an empty scenario
    Scenario config_scenario;

    // Initialize the dimensions for writing analogs
    array<string, 4> dims = {Config::_DIM_STATIONS, Config::_DIM_TEST_TIMES, Config::_DIM_FLTS, Config::_DIM_ANALOGS};

    // Open a NetCDF file to write to
    NcFile nc(file_path, NcFile::FileMode::write, NcFile::FileFormat::nc4);

    // Allocate memory to store simulated ensemble results
    Array4DPointer ac_array(num_stations, _NUM_DAYS, _NUM_HOURS, num_analogs);
    ssc_number_t dn[_NUM_TIMES], df[_NUM_TIMES], tdry[_NUM_TIMES], wspd[_NUM_TIMES];

    if (verbose >= Verbose::Progress) cout << "There are in total " << num_scenarios << " scenarios to simulate." << endl;
    
    for (size_t scenario_i = 0; scenario_i < num_scenarios; ++scenario_i) {

        if (verbose >= Verbose::Progress) cout << "Simulating scenario " << scenario_i << "/" << num_scenarios << endl;

        /*
         * Initialize data container and the module
         */
        ssc_data_t data_container = ssc_data_create();
        ssc_module_t module = ssc_module_create("pvwattsv5");

        // Return error if the module was not created successfully
        if (module == NULL) {
            ssc_data_free(data_container);
            throw runtime_error("could not create pvwattsv5 module");
        }

        // Set fixed configuration
        for (auto element : config_fixed) {
            ssc_data_set_number(data_container, element.first.c_str(), element.second);
        }

        // Get the current scenario
        getScenario(config_scenario, config_scenarios, scenario_i);

        // Set scenario configuration
        for (auto element : config_scenario) {
            ssc_data_set_number(data_container, element.first.c_str(), element.second);
        }

        /*
         * Simulations for stations and analog members are carried out individually 
         */
        for (size_t station_i = 0; station_i < num_stations; ++station_i) {
            for (size_t analog_i = 0; analog_i < num_analogs; ++analog_i) {

                /*
                 * Subset arrays to the current analog member and convert from double to float
                 */
                size_t pos = 0;

                for (auto test_time_i : all_test_times) {
                    for (auto flt_i : all_flts) {

                        dn[pos] = dn_arr.getValue(station_i, test_time_i, flt_i, analog_i);
                        df[pos] = df_arr.getValue(station_i, test_time_i, flt_i, analog_i);
                        tdry[pos] = t_arr.getValue(station_i, test_time_i, flt_i, analog_i);
                        wspd[pos] = wspd_arr.getValue(station_i, test_time_i, flt_i, analog_i);

                        pos++;
                    }
                }


                /*
                 * Create weather data table
                 */
                ssc_data_t weather_table = ssc_data_create();

                // Assign data arrays to the table
                ssc_data_set_array(weather_table, "dn", dn, _NUM_TIMES);
                ssc_data_set_array(weather_table, "df", df, _NUM_TIMES);
                ssc_data_set_array(weather_table, "tdry", tdry, _NUM_TIMES);
                ssc_data_set_array(weather_table, "wspd", wspd, _NUM_TIMES);

                // Set weather table
                ssc_data_set_table(data_container, "solar_resource_data", weather_table);

                // Free the table after it is assigned
                ssc_data_free(weather_table);


                /*
                 * Run simulation
                 */
                if (ssc_module_exec(module, data_container) == 0) {
                    ssc_module_free(module);
                    ssc_data_free(data_container);
                    throw runtime_error("pvwattsv5 simulation failed");
                }

                
                /*
                 * Copy results to Array4D
                 */
                int len = 0;
                ssc_number_t *results = ssc_data_get_array(data_container, output_var.c_str(), &len);

                if (len == _NUM_TIMES) {
                    size_t pos = 0;

                    for (auto test_time_i : all_test_times) {
                        for (auto flt_i : all_flts) {
                            ac_array.setValue(results[pos], station_i, test_time_i, flt_i, analog_i);
                            pos++;
                        }
                    }

                } else {
                    throw runtime_error("pvwattsv5 failed");
                }

            } // End loop for analog members
        } // End loop for stations


        /*
         * Write the simulated ensemble results for the current scenario
         */
        stringstream var_name_ss;
        
        string padded_scenario_i = to_string(scenario_i);
        padded_scenario_i.insert(padded_scenario_i.begin(), max_width - padded_scenario_i.size(), '0');
        
        var_name_ss << output_var << "_" << padded_scenario_i;
        Ncdf::writeArray4D(nc, ac_array, var_name_ss.str(), dims);

        // Write scenarios as variable attributes
        if (verbose >= Verbose::Progress) cout << "Writing results to " << file_path << endl;

        auto var = nc.getVar(var_name_ss.str());
        for (const auto & element : config_scenario) Ncdf::writeAttribute(var, element.first, element.second, NcType::nc_FLOAT);
        for (const auto & element : config_fixed) Ncdf::writeAttribute(var, element.first, element.second, NcType::nc_FLOAT);


        /*
         * Housekeeping
         */
        ssc_module_free(module);
        ssc_data_free(data_container);
        if (verbose >= Verbose::Progress) cout << "Simulations complete!" << endl;

    } // End of loop for scenarios

    return;
}

int main() {
    cout << "evergreen -- our pursuit for a more sustainable future" << endl;

    string file_path("/home/graduate/wuh20/storage/data/NAM/anen_results.nc");
    run_pvwattsv5(file_path);

    return 0;
}
