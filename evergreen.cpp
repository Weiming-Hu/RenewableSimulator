/* 
 *
 * File:   evergreen.cpp
 * Author: Weiming Hu <weiming@psu.edu>
 *
 * Created on February, 2020, 11:10 AM
 */

/** @file */

#include <string>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iterator>

#include "sscapi.h"

#include "Times.h"

using namespace std;

static const int _NUM_TIMES = 8760;


void printVersion() {
    
    int version = ssc_version();

    const char *build = ssc_build_info();
    cout << "SSC version " << version << ": " << build << endl;
    
    return;
}

void test_pvwattsv5() {
    cout << "Test module pvwattsv5" << endl;
    
    /*
     * Prepare array data
     */
    ssc_number_t dn[_NUM_TIMES], df[_NUM_TIMES], tdry[_NUM_TIMES], wspd[_NUM_TIMES];
    
    fill(dn, dn + _NUM_TIMES, 3.5);
    fill(df, df + _NUM_TIMES, 23);
    fill(tdry, tdry + _NUM_TIMES, 27);
    fill(wspd, wspd + _NUM_TIMES, 4);
    
    /*
     * Create weather data table
     */
    ssc_data_t weather_table = ssc_data_create();
    
    // Assign data arrays to the table
    ssc_data_set_array(weather_table, "dn", dn, _NUM_TIMES);
    ssc_data_set_array(weather_table, "df", df, _NUM_TIMES);
    ssc_data_set_array(weather_table, "tdry", tdry, _NUM_TIMES);
    ssc_data_set_array(weather_table, "wspd", wspd, _NUM_TIMES);
    
    ssc_data_set_number(weather_table, "lat", 40);
    ssc_data_set_number(weather_table, "lon", 70);
    ssc_data_set_number(weather_table, "tz", -8);

    
    /*
     * Create the data container
     */
    ssc_data_t data_container = ssc_data_create();
    
    ssc_data_set_number(data_container, "system_capacity", 1);
    ssc_data_set_number(data_container, "losses", 0);
    ssc_data_set_number(data_container, "array_type", 0);
    ssc_data_set_number(data_container, "tilt", 20);
    ssc_data_set_number(data_container, "azimuth", 180);
    ssc_data_set_number(data_container, "adjust:constant", 0);

    // Weather table is the different one from others
    ssc_data_set_table(data_container, "solar_resource_data", weather_table);
    
    // Free the table after it is assigned
    ssc_data_free(weather_table);
    
    
    /*
     * Run simulation
     */
    ssc_module_t module = ssc_module_create("pvwattsv5");
    
    // Return error if the module was not created successfully
    if (module == NULL) {
        ssc_data_free(data_container);
        throw runtime_error("could not create pvwattsv5 module");
    }
    
    // Run module
    if (ssc_module_exec(module, data_container) == 0) {
        ssc_module_free(module);
        ssc_data_free(data_container);
        throw runtime_error("pvwattsv5 simulation failed");
    }
    
    
    /*
     * Retrieve results
     */
    int len = 0;
    ssc_number_t *ac = ssc_data_get_array(data_container, "ac", &len);
    
    if (ac == NULL) {
        throw runtime_error("variable 'ac' not found");
    } else {
        cout << "ac: ";
        
        int count = 0;
        
        for (int i = 0; i < len; ++i) {
            cout << ac[i] << ',';
            count++;
            
            if (count == 24) {
                cout << endl;
                count = 0;
            }
            
        }
    }
    
    
    /*
     * Housekeeping
     */
    ssc_module_free(module);
    ssc_data_free(data_container);
    
    return;
}


int main() {
    cout << "evergreen -- our pursuit for a more sustainable future" << endl;
    
    test_pvwattsv5();

    Time time(100);
    
    return 0;
}
