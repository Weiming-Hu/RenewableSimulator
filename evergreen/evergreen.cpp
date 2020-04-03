/* 
 *
 * File:   evergreen.cpp
 * Author: Weiming Hu <weiming@psu.edu>
 *
 * Created on February, 2020, 11:10 AM
 */

/** @file */

#include <cmath>
#include <netcdf>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <boost/program_options.hpp>

#include "Ncdf.h"
#include "Config.h"
#include "Scenarios.h"
#include "Functions.h"
#include "AnEnContainer.h"
#include "FunctionsEvergreen.h"

using namespace std;
using namespace boost::program_options;

// Radians to degrees conversion factor
static const double R2D = 180 / M_PI;

void run_pvwatts(const string & file_path, Verbose verbose) {

    /*
     * Set up scenarios to simulate
     */
    if (verbose >= Verbose::Progress) cout << "Setting up simulation scenarios ..." << endl;

    // Define our batches of scenarios to simulate.
    // Scenarios will be created using permutation of all keys.
    //
    Scenarios pvwatts_scenarios;
    pvwatts_scenarios["system_capacity"] = {1, 2, 3, 4};   // System size (DC nameplate) in kW
    pvwatts_scenarios["module_type"] = {0, 1, 2};          // 0 for standard, 1 for premium, and 2 for thin film
    pvwatts_scenarios["losses"] = {0, 20, 40, 80, 99};     // Percent of total system loss within [-5, 99]
    
    // Define our fixed configuration
    Scenario pvwatts_config;
    pvwatts_config["time_step"] = 1;            // SSC default
    pvwatts_config["dc_ac_ratio"] = 1.1;        // SSC default
    pvwatts_config["inv_eff"] = 96;             // SSC default
    pvwatts_config["array_type"] = 0;           // Fixed
    pvwatts_config["gcr"] = 0.4;                // SSC default ground coverage ratio
    pvwatts_config["tilt"] = 0;                 // Horizontal


    /*
     * Read analog input
     */
    AnEnContainer anen_input(file_path, {
        "DownwardShortwaveRadiation", "Temperature_2m", "WindSpeed_10m", "Albedo"
    }, verbose);

    size_t num_stations = anen_input.stations().size();
    size_t num_days = anen_input.times().size();
    size_t num_flts = anen_input.flts().size();
    size_t num_analogs = anen_input.num_analogs();


    /*
     * Initialize data container and modules
     */
    ssc_data_t pvwatts_data = ssc_data_create();
    ssc_module_t pvwatts_module = ssc_module_create("pvwattsv5_1ts");

    // Return error if any modules were not created successfully
    if (pvwatts_module == NULL) throw runtime_error("Failed to create all SSC modules");

    // Take care of no-provide attributes
    ssc_number_t inout_tcell = NAN, inout_poa = NAN;

    // Set the fixed configurations
    pvwatts_config.set(pvwatts_data);
    

    /*
     * Run batch simulations
     */
    size_t num_scenarios = pvwatts_scenarios.totalScenarios();
    ssc_number_t ac_scalar, dc_scalar, azimuth_rad, elevation_rad, julian_day;
    
    Array4DPointer ac(num_stations, num_days, num_flts, num_analogs);
    Array4DPointer dc(num_stations, num_days, num_flts, num_analogs);
    Array4DPointer dni(num_stations, num_days, num_flts, num_analogs);
    Array4DPointer dhi(num_stations, num_days, num_flts, num_analogs);

    array<string, 4> array_dims {Config::_DIM_STATIONS, Config::_DIM_TEST_TIMES, Config::_DIM_FLTS, Config::_DIM_ANALOGS};
    
    // Open the file for writing
    netCDF::NcFile nc(file_path, netCDF::NcFile::FileMode::write);

    if (verbose >= Verbose::Progress) cout << "There are in total " << num_scenarios << " scenarios, "
            << num_stations << " stations, " << num_days << " days, " << num_flts << " lead times to simulate." << endl;

    for (size_t scenario_i = 0; scenario_i < num_scenarios; ++scenario_i) {

        // Set the current scenario
        if (verbose >= Verbose::Progress) cout << "Simulating scenario " << scenario_i << "/" << num_scenarios << endl;
        pvwatts_scenarios.set(pvwatts_data, scenario_i);

        // Looping through all stations, days, and lead times.
        // Each iteration simulates all members at once.
        //
        for (size_t station_i = 0; station_i < num_stations; ++station_i) {

            // Get longitude and latitude for the current station
            ssc_number_t lon, lat;
            anen_input.getCoordinate(lon, lat, station_i);
            
            // Set coordinates
            ssc_data_set_number(pvwatts_data, "lon", lon);
            ssc_data_set_number(pvwatts_data, "lat", lat);
            
            // Set time zone
            ssc_number_t tz = anen_input.getTimeZone(station_i);
            ssc_data_set_number(pvwatts_data, "tz", tz);

            for (size_t day_i = 0; day_i < num_days; ++day_i) {
                for (size_t flt_i = 0; flt_i < num_flts; ++flt_i) {
                    
                    // Get current time in local time zone and UTC
                    ssc_number_t local_year, local_month, local_day, local_hour, local_minute;
                    ssc_number_t UTC_year, UTC_month, UTC_day, UTC_hour, UTC_minute;
                    
                    anen_input.getLocalDateTime(local_year, local_month, local_day, local_hour, local_minute, day_i, flt_i, tz);
                    anen_input.getUTCDateTime(UTC_year, UTC_month, UTC_day, UTC_hour, UTC_minute, day_i, flt_i);

                    // Set local time to pvwatts module
                    ssc_data_set_number(pvwatts_data, "year", local_year);
                    ssc_data_set_number(pvwatts_data, "month", local_month);
                    ssc_data_set_number(pvwatts_data, "day", local_day);
                    ssc_data_set_number(pvwatts_data, "hour", local_hour);
                    ssc_data_set_number(pvwatts_data, "minute", local_minute);

                    if (verbose >= Verbose::Debug) cout << "Local date time: " << local_year << "/" << local_month << "/" << local_day
                            << " " << local_hour << ":" << local_minute << " " << tz << "; UTC date time: " << UTC_year << "/"
                            << UTC_month << "/" << UTC_day << " " << UTC_hour << ":" << UTC_minute << endl;
                    
                    // Calculate current azimuth and other variables
                    FunctionsEvergreen::sun_pos(azimuth_rad, elevation_rad, julian_day, lon, lat,
                            UTC_year, UTC_month, UTC_day, UTC_hour, 0, 0);
                    
                    // Set azimuth
                    ssc_data_set_number(pvwatts_data, "azimuth", azimuth_rad * R2D);

                    for (size_t analog_i = 0; analog_i < num_analogs; ++analog_i) {
                        
                        // Get data from analogs
                        ssc_number_t alb = anen_input.getValue(station_i, day_i, flt_i, analog_i, "Albedo");
                        ssc_number_t wspd = anen_input.getValue(station_i, day_i, flt_i, analog_i, "WindSpeed_10m");
                        ssc_number_t tamb = anen_input.getValue(station_i, day_i, flt_i, analog_i, "Temperature_2m");
                        ssc_number_t ghi = anen_input.getValue(station_i, day_i, flt_i, analog_i, "DownwardShortwaveRadiation");

                        /*
                         * Hard coded operation !!!!
                         */
                        // temparature unit is in Kelvin. Convert it to Celsius
                        tamb -= 273.15;
                        
                        // Albedo is a percentage. Convert it to a ratio.
                        alb /= 100;

                        /*
                         * Step 1: Reindl et al. decomposition from GHI to DHI and DNI
                         */
                        ssc_number_t dhi_scalar, dni_scalar;
                        FunctionsEvergreen::decompose_ghi(ghi, dhi_scalar, dni_scalar, elevation_rad, julian_day);
                        
                        // Save results
                        dni.setValue(dni_scalar, station_i, day_i, flt_i, analog_i);
                        dhi.setValue(dhi_scalar, station_i, day_i, flt_i, analog_i);
                        
                        if (std::isnan(dhi_scalar) || std::isnan(dni_scalar)) {
                            ac.setValue(NAN, station_i, day_i, flt_i, analog_i);
                            dc.setValue(NAN, station_i, day_i, flt_i, analog_i);
                            continue;
                        }

                        
                        /*
                         * Step 2: PV power output simulation from DHI and DNI
                         */
                        ssc_data_set_number(pvwatts_data, "beam", dni_scalar);
                        ssc_data_set_number(pvwatts_data, "diffuse", dhi_scalar);
                        ssc_data_set_number(pvwatts_data, "alb", alb);
                        ssc_data_set_number(pvwatts_data, "tamb", tamb);
                        ssc_data_set_number(pvwatts_data, "wspd", wspd);

                        // Set no-provide variables. Setting them to NA avoids the impact from last iteration
                        ssc_data_set_number(pvwatts_data, "tcell", inout_tcell);
                        ssc_data_set_number(pvwatts_data, "poa", inout_poa);
                        
                        // Simulate AC and DC
                        if (ssc_module_exec(pvwatts_module, pvwatts_data) == 0) {
                            if (verbose >= Verbose::Warning) cerr << "Error from pvwatts at iteration ["
                                    << station_i << "," << day_i << "," << flt_i << "," << analog_i << "] with the input data: "
                                    << FunctionsEvergreen::toString(pvwatts_data, "pvwattsv5_1ts") << endl;
                            
                            ac.setValue(NAN, station_i, day_i, flt_i, analog_i);
                            dc.setValue(NAN, station_i, day_i, flt_i, analog_i);
                            
                            continue;
                        } else {
                            
                            // Get AC and DC
                            ssc_data_get_number(pvwatts_data, "ac", &ac_scalar);
                            ssc_data_get_number(pvwatts_data, "dc", &dc_scalar);

                            // Assign values
                            ac.setValue(ac_scalar, station_i, day_i, flt_i, analog_i);
                            dc.setValue(dc_scalar, station_i, day_i, flt_i, analog_i);
                        }

                    } // End loop for analog members
                } // End loop for lead times
            } // End loop for days
        } // End loop for stations

        // Create a group for this scenario
        stringstream scenario_name;
        scenario_name << "scenario_" << scenario_i;
        auto nc_group = nc.addGroup(scenario_name.str());
        
        Ncdf::writeArray4D(nc_group, ac, "ac", array_dims);
        Ncdf::writeArray4D(nc_group, dc, "dc", array_dims);
        Ncdf::writeArray4D(nc_group, dni, "dni", array_dims);
        Ncdf::writeArray4D(nc_group, dhi, "dhi", array_dims);
        
        pvwatts_scenarios.write(nc_group, scenario_i);
        pvwatts_config.write(nc_group);
    }
    
    ssc_module_free(pvwatts_module);
    ssc_data_free(pvwatts_data);

    return;
}

int main(int argc, char** argv) {

    // Define variables that should be set from command line arguments
    string file_path, out_var;
    Verbose verbose;
    int verbose_int;

    // Define a welcoming message
    stringstream msg;
    msg << "evergreen -- our pursuit of a more sustainable future" << endl
            << "Developed by Weiming Hu [weiming-hu.github.io]" << endl
            << "Issues @ https://github.com/Weiming-Hu/RenewableSimulator/issues" << endl;

    // Define available options
    options_description desc("Available options");
    desc.add_options()
            ("help,h", "Print help information for options")
            ("anen", value<string>(&file_path)->required(), "The NetCDF file for analogs")
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

    run_pvwatts(file_path, verbose);

    return 0;
}
