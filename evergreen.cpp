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
#include <cmath>


#include "Ncdf.h"

#include "Functions.h"
#include "AnEnReadNcdf.h"
#include "Array4DPointer.h"

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"

using namespace std;
using namespace netCDF;
using namespace boost::program_options;

static const int _NUM_HOURS = 24;
static const int _NUM_DAYS = 365;
static const int _NUM_TIMES = _NUM_DAYS * _NUM_HOURS;
//
//
//int
//simpleTimeZone(double lon) {
//    return round(lon * 24 / 360);
//}
//
///**
// * Extract a certain scenario from the predefined scenarios. The scenario is
// * a combination of different configurations in the predefined scenarios. The
// * index is the scenario index. You should use countScenarios first to see
// * the total number of scenarios available to extract.
// * 
// * @param scenario A Scenario. It is essentially a string-number STL map
// * @param scenarios A Scenarios. It is essentially a string-vector STL map
// * @param index The scenario ID to extract. You can use countScenarios
// * to get the total number of scenarios available.
// */
//void
//getScenario(Scenario & scenario, const Scenarios & scenarios, size_t index) {
//
//    // Clear map
//    scenario.clear();
//
//    // Initialize values
//    size_t accum_len = 1;
//    size_t current_index;
//
//    for (auto e : scenarios) {
//
//        // Calculate the current index for this key based on the global index
//        if (accum_len == 1) {
//            accum_len = e.second.size();
//            current_index = fmod(index, accum_len);
//        } else {
//            current_index = fmod(index / accum_len, e.second.size());
//            accum_len *= e.second.size();
//        }
//
//        // Insert pair
//        scenario.insert(make_pair(e.first, e.second[current_index]));
//    }
//
//    return;
//}
//
//size_t
//countScenarios(const Scenarios & map) {
//    size_t total = 1;
//    for (auto e : map) total *= e.second.size();
//    return total;
//}
//
//class AnEnContainer {
//    
////    Times times;
////    Stations stations;
////    Array4DPointer dn_arr, df_arr, t_arr, wspd_arr;
////    AnEnReadNcdf anen_read(verbose);
////    
////    // Variable names are hard coded
////    NcFile nc(file_path, NcFile::FileMode::read);
////    anen_read.read(nc, times, "test_times");
////    anen_read.read(nc, stations);
////    
////    anen_read.readAnalogs(file_path, dn_arr, "DownwardShortwaveRadiation");
////    anen_read.readAnalogs(file_path, df_arr, "UpwardShortwaveRadiation");
////    anen_read.readAnalogs(file_path, t_arr, "temperature_2m");
////    anen_read.readAnalogs(file_path, wspd_arr, "wspd_1000hPa");
////    
////    // Sanity checks
////    if (stations.size() != dn_arr.shape()[0]) throw runtime_error("Stations do not match the array shape");
////    if (times.size() != dn_arr.shape()[1]) throw runtime_error("Times do not match the array shape");
////    
////    // The time span of analogs should not be longer than a year
////    double first_time = times.getTime(0).timestamp;
////    double last_time = times.getTime(times.size() - 1).timestamp;
////    if (last_time - first_time > 365 * )
////    
////    /*
////     * Input alignment
////     * 
////     * NREL SAM simulator requires hourly data for the entire year. We need to align
////     * input analogs with the year-round hourly time series.
////     */
////    if (times.size() == _NUM_TIMES) {
////        // This is the ideal case. I don't need to align the values
////    } else {
////        // There are probably missing days in the input analogs.
////        // I need to create a 
////    }
//};
//
//
//void
//run_pvwattsv5(const string & file_path, Verbose verbose) {
//    
//    /*
//     * Set up scenarios to simulate
//     */
//    if (verbose >= Verbose::Progress) cout << "Setting up pvwattsv5 simulation ..." << endl;
//
//    // Define our batches of scenarios to simulate.
//    // Scenarios will be created using permutation of all keys.
//    //
//    Scenarios config_scenarios;
//
//    config_scenarios["losses"] = {0, 0.2};
//    config_scenarios["azimuth"] = {0};
//    config_scenarios["system_capacity"] = {0.8};
//    config_scenarios["tilt"] = {10};
//
//    // Define our fixed configuration
//    Scenario config_fixed;
//
//    config_fixed["array_type"] = 0;
//    config_fixed["adjust:constant"] = 0;
//    config_fixed["verbose"] = 0;
//    
//    
//    /*
//     * Read analog input
//     */
//    AnEnContainer anen_input(file_path);
//    
//    
//    /*
//     * Run batch simulations
//     */
//    size_t num_scenarios = 10;
//    for (size_t scenario_i = 0; scenario_i < num_scenarios; ++scenario_i) {
//        
//        
//        /*
//         * Initialize data container and the module
//         */
//        ssc_data_t data_container = ssc_data_create();
//        ssc_module_t module = ssc_module_create("pvwattsv5");
//
//        // Return error if the module was not created successfully
//        if (module == NULL) {
//            ssc_data_free(data_container);
//            throw runtime_error("could not create pvwattsv5 module");
//        }
//        
//        config_fixed.set(data_container)
//        config_scenarios.set(data_container, scenario_i);
//        
//        
//        /*
//         * Simulations for stations and analog members are carried out individually 
//         */
//        for (size_t station_i = 0; station_i < num_stations; ++station_i) {
//            
//            // Set time offset
//            anen_input.getTimezone(station_i);
//            
//            for (size_t analog_i = 0; analog_i < num_analogs; ++analog_i) {
//
//                /*
//                 * Subset arrays to the current analog member and convert from double to float
//                 */
//                anen_input.subset(ptr_map, station_i, flt_i, analog_i);
//                
//                AnEnContainer::set(ptr_map, data_container);
//
//                
//                /*
//                 * Run simulation
//                 */
//                // TODO: disable model standard output
//                if (ssc_module_exec(module, data_container) == 0) {
//                    ssc_module_free(module);
//                    ssc_data_free(data_container);
//                    throw runtime_error("pvwattsv5 simulation failed");
//                }
//
//
//                /*
//                 * Copy results to Array4D
//                 */
//                AnEnContainer::copyValues();
//
//            } // End loop for analog members
//        } // End loop for stations
//
//        writeAnalogs();
//        writeScenarios();
//    }
//    
//    writeAnEnContainer();
//    
//}
//
//
//void
//run_pvwattsv5_old(const string & file_path, string output_var, Verbose verbose) {
//
//    /*
//     * Set up configuration schemes
//     */
//    if (verbose >= Verbose::Progress) cout << "Setting up pvwattsv5 simulation ..." << endl;
//
//    // Define our batches of scenarios to simulate.
//    // Scenarios will be created using permutation of all keys.
//    //
//    Scenarios config_scenarios;
//
//    config_scenarios["losses"] = {0, 0.2};
//    config_scenarios["azimuth"] = {0};
//    config_scenarios["system_capacity"] = {0.8};
//    config_scenarios["tilt"] = {10};
//
//    // Define our fixed configuration
//    Scenario config_fixed;
//
//    config_fixed["array_type"] = 0;
//    config_fixed["adjust:constant"] = 0;
//    config_fixed["verbose"] = 0;
//
//    
//    /*
//     * Read analog input
//     */
//    Array4DPointer dn_arr, df_arr, t_arr, wspd_arr;
//
//    // TODO: Change variable name
//    AnEnReadNcdf anen_read(verbose);
//    anen_read.readAnalogs(file_path, dn_arr, "temperature_2m");
//    anen_read.readAnalogs(file_path, df_arr, "wspd_1000hPa");
//    anen_read.readAnalogs(file_path, t_arr, "temperature_2m");
//    anen_read.readAnalogs(file_path, wspd_arr, "wspd_1000hPa");
//    
//    // TODO: Change this
//    for (size_t i = 0; i < dn_arr.num_elements(); ++i) {
//        dn_arr.getValuesPtr()[i] = 5;
//        df_arr.getValuesPtr()[i] = 3;
//        t_arr.getValuesPtr()[i] = 3;
//        wspd_arr.getValuesPtr()[i] = 3;
//    }
//
////    if (dn_arr.shape()[1] != _NUM_DAYS) throw runtime_error("Analogs must have 365 days for test times");
////    if (dn_arr.shape()[2] != _NUM_HOURS) throw runtime_error("Analogs must have 24 hours for forecast lead times");
//
//    /*
//     * Run simulation for all
//     * - stations
//     * - ensemble members
//     * - scenarios
//     */
//    size_t num_stations = dn_arr.shape()[0];
//    size_t num_analogs = dn_arr.shape()[3];
//    size_t num_scenarios = countScenarios(config_scenarios);
//
//    // Calculate how wide should I pad the string to make pretty variable names
//    size_t max_width = to_string(num_scenarios).size();
//
//    // Create indices for dimensions that are not subset
//    vector<size_t> all_test_times(_NUM_DAYS), all_flts(_NUM_HOURS);
//    iota(all_test_times.begin(), all_test_times.end(), 0);
//    iota(all_flts.begin(), all_flts.end(), 0);
//
//    // Initialize an empty scenario
//    Scenario config_scenario;
//
//    // Initialize the dimensions for writing analogs
//    array<string, 4> dims = {Config::_DIM_STATIONS, Config::_DIM_TEST_TIMES, Config::_DIM_FLTS, Config::_DIM_ANALOGS};
//
//    // Open a NetCDF file to write to
//    NcFile nc(file_path, NcFile::FileMode::write, NcFile::FileFormat::nc4);
//
//    // Allocate memory to store simulated ensemble results
//    Array4DPointer ac_array(num_stations, _NUM_DAYS, _NUM_HOURS, num_analogs);
//    ssc_number_t dn[_NUM_TIMES], df[_NUM_TIMES], tdry[_NUM_TIMES], wspd[_NUM_TIMES];
//
//    if (verbose >= Verbose::Progress) cout << "There are in total " << num_scenarios << " scenarios to simulate." << endl;
//
//    for (size_t scenario_i = 0; scenario_i < num_scenarios; ++scenario_i) {
//
//        if (verbose >= Verbose::Progress) cout << "Simulating scenario " << scenario_i << "/" << num_scenarios << endl;
//
//        /*
//         * Initialize data container and the module
//         */
//        ssc_data_t data_container = ssc_data_create();
//        ssc_module_t module = ssc_module_create("pvwattsv5");
//
//        // Return error if the module was not created successfully
//        if (module == NULL) {
//            ssc_data_free(data_container);
//            throw runtime_error("could not create pvwattsv5 module");
//        }
//
//        // Set fixed configuration
//        for (auto element : config_fixed) {
//            ssc_data_set_number(data_container, element.first.c_str(), element.second);
//        }
//
//        // Get the current scenario
//        getScenario(config_scenario, config_scenarios, scenario_i);
//
//        // Set scenario configuration
//        for (auto element : config_scenario) {
//            ssc_data_set_number(data_container, element.first.c_str(), element.second);
//        }
//
//        /*
//         * Simulations for stations and analog members are carried out individually 
//         */
//        for (size_t station_i = 0; station_i < num_stations; ++station_i) {
//            for (size_t analog_i = 0; analog_i < num_analogs; ++analog_i) {
//
//                /*
//                 * Subset arrays to the current analog member and convert from double to float
//                 */
//                size_t pos = 0;
//
//                for (auto test_time_i : all_test_times) {
//                    for (auto flt_i : all_flts) {
//
//                        dn[pos] = dn_arr.getValue(station_i, test_time_i, flt_i, analog_i);
//                        df[pos] = df_arr.getValue(station_i, test_time_i, flt_i, analog_i);
//                        tdry[pos] = t_arr.getValue(station_i, test_time_i, flt_i, analog_i);
//                        wspd[pos] = wspd_arr.getValue(station_i, test_time_i, flt_i, analog_i);
//
//                        pos++;
//                    }
//                }
//
//
//                /*
//                 * Create weather data table
//                 */
//                ssc_data_t weather_table = ssc_data_create();
//
//                // Assign data arrays to the table
//                ssc_data_set_array(weather_table, "dn", dn, _NUM_TIMES);
//                ssc_data_set_array(weather_table, "df", df, _NUM_TIMES);
//                ssc_data_set_array(weather_table, "tdry", tdry, _NUM_TIMES);
//                ssc_data_set_array(weather_table, "wspd", wspd, _NUM_TIMES);
//
//                // Set weather table
//                ssc_data_set_table(data_container, "solar_resource_data", weather_table);
//
//                // Free the table after it is assigned
//                ssc_data_free(weather_table);
//
//
//                /*
//                 * Run simulation
//                 */
//                // TODO: disable model standard output
//                if (ssc_module_exec(module, data_container) == 0) {
//                    ssc_module_free(module);
//                    ssc_data_free(data_container);
//                    throw runtime_error("pvwattsv5 simulation failed");
//                }
//
//
//                /*
//                 * Copy results to Array4D
//                 */
//                int len = 0;
//                ssc_number_t *results = ssc_data_get_array(data_container, output_var.c_str(), &len);
//
//                if (len == _NUM_TIMES) {
//                    size_t pos = 0;
//
//                    for (auto test_time_i : all_test_times) {
//                        for (auto flt_i : all_flts) {
//                            ac_array.setValue(results[pos], station_i, test_time_i, flt_i, analog_i);
//                            pos++;
//                        }
//                    }
//
//                } else {
//                    throw runtime_error("pvwattsv5 failed");
//                }
//
//            } // End loop for analog members
//        } // End loop for stations
//
//
//        /*
//         * Write the simulated ensemble results for the current scenario
//         */
//        stringstream var_name_ss;
//
//        string padded_scenario_i = to_string(scenario_i);
//        padded_scenario_i.insert(padded_scenario_i.begin(), max_width - padded_scenario_i.size(), '0');
//
//        var_name_ss << output_var << "_" << padded_scenario_i;
//        Ncdf::writeArray4D(nc, ac_array, var_name_ss.str(), dims);
//
//        // Write scenarios as variable attributes
//        if (verbose >= Verbose::Progress) cout << "Writing results to " << file_path << endl;
//
//        auto var = nc.getVar(var_name_ss.str());
//        for (const auto & element : config_scenario) Ncdf::writeAttribute(var, element.first, element.second, NcType::nc_FLOAT);
//        for (const auto & element : config_fixed) Ncdf::writeAttribute(var, element.first, element.second, NcType::nc_FLOAT);
//
//
//        /*
//         * Housekeeping
//         */
//        ssc_module_free(module);
//        ssc_data_free(data_container);
//
//    } // End of loop for scenarios
//
//    if (verbose >= Verbose::Progress) cout << "Simulations complete!" << endl;
//    return;
//}
//
//int main(int argc, char** argv) {
//
//    // Define variables that should be set from command line arguments
//    string file_path, out_var;
//    Verbose verbose;
//    int verbose_int;
//
//    // Define available options
//    options_description desc("Available options");
//    desc.add_options()
//            ("help,h", "Print help information for options")
//            ("anen", value<string>(&file_path)->required(), "The NetCDF file for analogs")
//            ("output-var", value<string>(&out_var)->default_value("ac"), "The simulation variable to output from SSC")
//            ("verbose,v", value<int>(&verbose_int)->default_value(2), "The verbose level (0 - 4)");
//
//    // Parse the command line arguments
//    variables_map vm;
//    store(parse_command_line(argc, argv, desc), vm);
//
//    if (vm.count("help") || argc == 1) {
//        cout << desc << endl;
//        return 0;
//    }
//
//    notify(vm);
//
//    // Convert verbose
//    verbose = Functions::itov(verbose_int);
//    if (verbose >= Verbose::Progress) cout << "evergreen -- our pursuit for a more sustainable future" << endl;
//
//    run_pvwattsv5(file_path, out_var, verbose);
//
//    return 0;
//}

int main() {
    return 0;
}