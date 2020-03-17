#include <string>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "sscapi.h"


using namespace std;

static const int _NUM_TIMES = 8760;

void test_irradproc() {
    cout << "Test module irradproc" << endl;

    // Create module
    ssc_module_t module = ssc_module_create("irradproc");
    if (module == NULL) throw runtime_error("could not create irradproc");

    // Create an empty data container
    ssc_data_t data_container = ssc_data_create();

    // Assign values to the data container
    ssc_data_set_number(data_container, "irrad_mode", 2);
    ssc_data_set_number(data_container, "albedo_const", 0.2);
    ssc_data_set_number(data_container, "lat", 70);
    ssc_data_set_number(data_container, "lon", -40);
    ssc_data_set_number(data_container, "tz", -4);
    ssc_data_set_number(data_container, "sky_model", 0);
    ssc_data_set_number(data_container, "track_mode", 0);
    ssc_data_set_number(data_container, "azimuth", 180);
    ssc_data_set_number(data_container, "tilt", 0);
    ssc_data_set_number(data_container, "rotlim", 45);
    ssc_data_set_number(data_container, "backtrack", 0);
    ssc_data_set_number(data_container, "gcr", 0.4);

    static const int total = 3;

    ssc_number_t beam[total];

    ssc_number_t ghi[total] = {10, 300, 100};
    ssc_number_t dhi[total] = {10, 100, 12};

    //ssc_number_t ghi[total] = {300, 280};
    //ssc_number_t dhi[total] = {100, 120};

    ssc_number_t alb[total] = {0.2, 0.2, 0.2};
    ssc_number_t year[total] = {2010, 2010, 2010};
    ssc_number_t month[total] = {8, 8, 8};
    ssc_number_t day[total] = {1, 1, 1};
    ssc_number_t hour[total] = {12, 12, 12};
    ssc_number_t minute[total] = {0, 0, 0};

    ssc_data_set_array(data_container, "beam", beam, total);
    ssc_data_set_array(data_container, "diffuse", dhi, total);
    ssc_data_set_array(data_container, "global", ghi, total);
    ssc_data_set_array(data_container, "albedo", alb, total);
    ssc_data_set_array(data_container, "year", year, total);
    ssc_data_set_array(data_container, "month", month, total);
    ssc_data_set_array(data_container, "day", day, total);
    ssc_data_set_array(data_container, "hour", hour, total);
    ssc_data_set_array(data_container, "minute", minute, total);

    // Run simulation
    if (ssc_module_exec(module, data_container) == 0) {
        ssc_module_free(module);
        ssc_data_free(data_container);
        throw runtime_error("irradproc simulation failed");
    }

    // Retrieve results
    int len;
    ssc_number_t *poa_beam = ssc_data_get_array(data_container, "poa_beam", &len);
    ssc_number_t *poa_skydiff = ssc_data_get_array(data_container, "poa_skydiff", &len);
    ssc_number_t *poa_gnddiff = ssc_data_get_array(data_container, "poa_gnddiff", &len);

    // Print results
    cout << "poa_beam: ";
    for (int i = 0; i < len; ++i) cout << poa_beam[i] << ",";
    cout << endl;
    cout << "poa_skydiff: ";
    for (int i = 0; i < len; ++i) cout << poa_skydiff[i] << ",";
    cout << endl;
    cout << "poa_gnddiff: ";
    for (int i = 0; i < len; ++i) cout << poa_gnddiff[i] << ",";
    cout << endl;

    // Housekeeping
    ssc_module_free(module);
    ssc_data_free(data_container);

    return;
}

void test_pvwattsv5_1ts() {
    cout << "Test module pvwattsv5_1ts" << endl;

    // Create module
    ssc_module_t module = ssc_module_create("pvwattsv5_1ts");
    if (module == NULL) throw runtime_error("could not create pvwattsv5_1ts");

    // Create an empty data container
    ssc_data_t data_container = ssc_data_create();

    // Assign values to the data container
    ssc_data_set_number(data_container, "year", 2010);
    ssc_data_set_number(data_container, "month", 8);
    ssc_data_set_number(data_container, "day", 1);
    ssc_data_set_number(data_container, "hour", 12);
    ssc_data_set_number(data_container, "minute", 0);
    ssc_data_set_number(data_container, "lat", 70);
    ssc_data_set_number(data_container, "lon", -40);
    ssc_data_set_number(data_container, "tz", -4);
    ssc_data_set_number(data_container, "tamb", 30);
    ssc_data_set_number(data_container, "wspd", 5);
    ssc_data_set_number(data_container, "alb", 0.3);
    ssc_data_set_number(data_container, "time_step", 1);
    ssc_data_set_number(data_container, "system_capacity", 5000);
    ssc_data_set_number(data_container, "module_type", 0);
    ssc_data_set_number(data_container, "dc_ac_ratio", 1);
    ssc_data_set_number(data_container, "inv_eff", 90);
    ssc_data_set_number(data_container, "losses", 10);
    ssc_data_set_number(data_container, "array_type", 0);
    ssc_data_set_number(data_container, "tilt", 30);
    ssc_data_set_number(data_container, "azimuth", 180);
    ssc_data_set_number(data_container, "gcr", 0.4);

    ssc_data_set_number(data_container, "beam", 100);
    ssc_data_set_number(data_container, "diffuse", 50);

    ssc_number_t poa, tcell;
    ssc_data_set_number(data_container, "tcell", tcell);
    ssc_data_set_number(data_container, "poa", poa);

    // Run simulation
    if (ssc_module_exec(module, data_container) == 0) {
        ssc_module_free(module);
        ssc_data_free(data_container);
        throw runtime_error("pvwattsv5_1ts simulation failed");
    }

    // Retrieve results
    ssc_number_t ac, dc;
    ssc_data_get_number(data_container, "ac", &ac);
    ssc_data_get_number(data_container, "dc", &dc);
    ssc_data_get_number(data_container, "poa", &poa);

    // Print results
    cout << "(pvwattsv5_1ts) ac: " << ac << ", dc: " << dc << ", poa: " << poa << endl;

    // Housekeeping
    ssc_module_free(module);
    ssc_data_free(data_container);

    return;
}

void test_pvwattsv5() {
    cout << "Test module pvwattsv5" << endl;
    
    /*
     * Prepare array data
     */
    ssc_number_t dn[_NUM_TIMES], df[_NUM_TIMES], tdry[_NUM_TIMES], wspd[_NUM_TIMES], gh[_NUM_TIMES];
    
    fill(gh, gh + _NUM_TIMES, 60);
    fill(dn, dn + _NUM_TIMES, 35);
    fill(df, df + _NUM_TIMES, 5);
    fill(tdry, tdry + _NUM_TIMES, 27);
    fill(wspd, wspd + _NUM_TIMES, 4);
    
    dn[28] = NAN;
//    tdry[3] = NAN;
    
    /*
     * Create weather data table
     */
    ssc_data_t weather_table = ssc_data_create();
    
    // Assign data arrays to the table
    //ssc_data_set_array(weather_table, "df", df, _NUM_TIMES);
    ssc_data_set_array(weather_table, "gh", gh, _NUM_TIMES);
    ssc_data_set_array(weather_table, "tdry", tdry, _NUM_TIMES);
    ssc_data_set_array(weather_table, "wspd", wspd, _NUM_TIMES);
    
    ssc_data_set_number(weather_table, "lat", 40);
    ssc_data_set_number(weather_table, "lon", -70);
    ssc_data_set_number(weather_table, "tz", -8);

    
    /*
     * Create the data container
     */
    ssc_data_t data_container = ssc_data_create();
    
    ssc_data_set_number(data_container, "system_capacity", 1);
    ssc_data_set_number(data_container, "losses", 0);
    ssc_data_set_number(data_container, "array_type", 0);
    ssc_data_set_number(data_container, "tilt", 0);
    ssc_data_set_number(data_container, "azimuth", 100);
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

            if (count == 0) {
                cout << "[" << i << "]: ";
            }

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

    //test_pvwattsv5();
    test_pvwattsv5_1ts();
    //test_irradproc();

    return 0;
}
