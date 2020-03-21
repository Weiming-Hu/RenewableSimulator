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
#include <boost/program_options.hpp>

#include "Scenarios.h"
#include "Functions.h"
#include "AnEnContainer.h"

using namespace std;
using namespace boost::program_options;

void run_pvwatts() {

    // TODO: These are input
    Verbose verbose = Verbose::Progress;
    string file_path = "This is input";

    /*
     * Set up scenarios to simulate
     */
    if (verbose >= Verbose::Progress) cout << "Setting up simulation scenarios ..." << endl;

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


    /*
     * Read analog input
     */
    AnEnContainer anen_input(file_path);
    
    size_t num_stations = anen_input.stations().size();
    size_t num_days = anen_input.times().size();
    size_t num_flts = anen_input.flts().size();
    size_t num_analogs = anen_input.num_analogs();

    
    /*
     * Initialize data container and modules
     */
    ssc_data_t ssc_container = ssc_data_create();
    ssc_module_t irradproc = ssc_module_create("irradproc");
    ssc_module_t pvwatts = ssc_module_create("pvwattsv5_1ts");

    // Return error if any modules were not created successfully
    if (irradproc == NULL || pvwatts == NULL) {
        ssc_data_free(ssc_container);
        throw runtime_error("Failed to create all SSC modules");
    }
    
    // Take care of no-provide attributes
    ssc_number_t *skip_DHI = new ssc_number_t[num_analogs];
    ssc_number_t *inout_tcell = new ssc_number_t[num_analogs];
    ssc_number_t *inout_poa = new ssc_number_t[num_analogs];
    
    ssc_data_set_array(irradproc, "diffuse", skip_DHI, num_analogs);
    ssc_data_set_array(pvwatts, "tcell", inout_tcell, num_analogs);
    ssc_data_set_array(pvwatts, "poa", inout_poa, num_analogs);
    
    
    /*
     * Run batch simulations
     */
    size_t num_scenarios = config_scenarios.totalScenarios();
    
    if (verbose >= Verbose::Progress) cout
            << "There are in total " << num_scenarios << " scenarios, "
            << num_stations << " stations, " << num_days << " days, " 
            << num_flts << " lead times to simulate." << endl;
    
    for (size_t scenario_i = 0; scenario_i < num_scenarios; ++scenario_i) {

        // Set the current scenario
        config_fixed.set(ssc_container);
        config_scenarios.set(ssc_container, scenario_i);

        // Looping through all stations, days, and lead times.
        // Each iteration simulates all members at once.
        //
        for (size_t station_i = 0; station_i < num_stations; ++station_i) {
            
            // Set location based attributes
            
            
            for (size_t day_i = 0; day_i < num_days; ++day_i) {
                for (size_t flt_i = 0; flt_i < num_flts; ++flt_i) {

                    
                    /*
                     * Step 1: Reindl et al. decomposition from GHI to DHI
                     */
                    
                    
                    /*
                     * Step 2: irradproc module from GHI and DHI to DNI
                     */
                    
                    
                    /*
                     * Step 3: PV power output simulation from DHI and DNI
                     */
                    
                } // End loop for lead times
            } // End loop for days
        } // End loop for stations

//        writeAnalogs();
//        writeScenarios();
    }

    ssc_module_free(pvwatts);
    ssc_module_free(irradproc);
    ssc_data_free(ssc_container);
    
    delete [] skip_DHI;
    delete [] inout_poa;
    delete [] inout_tcell;
    
    return;
}



//
//
//void
//run_pvwattsv5(const string & file_path, Verbose verbose) {
//    
//    
//    
//    
//    
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

int main(int argc, char** argv) {

    // Define variables that should be set from command line arguments
    string file_path, out_var;
    Verbose verbose;
    int verbose_int;

    // Define a welcoming message
    stringstream msg;
    msg << "evergreen -- our pursuit of a more sustainable future" << endl
            << "Developed by Weiming Hu [weiming-hu.github.io]" << endl
            << "Issues @ https://github.com/Weiming-Hu/RenewableSimulator/issues";

    // Define available options
    options_description desc("Available options");
    desc.add_options()
            ("help,h", "Print help information for options")
            ("anen", value<string>(&file_path)->required(), "The NetCDF file for analogs")
            //            ("output-var", value<string>(&out_var)->default_value("ac"), "The simulation variable to output from SSC")
            ("verbose,v", value<int>(&verbose_int)->default_value(2), "The verbose level (0 - 4)");

    // Parse the command line arguments
    variables_map vm;
    store(parse_command_line(argc, argv, desc), vm);

    if (vm.count("help") || argc == 1) {
        cout << msg.str() << endl << endl << desc << endl;
        return 0;
    }

    notify(vm);

    // Convert verbose
    verbose = Functions::itov(verbose_int);
    if (verbose >= Verbose::Progress) cout << msg.str() << endl;

    run_pvwatts();

    return 0;
}
