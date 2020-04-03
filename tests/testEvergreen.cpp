/*
 * File:   testEvergreen.cpp
 * Author: Weiming Hu <weiming@psu.edu>
 *
 * Created on Aug 4, 2018, 4:09:20 PM
 */

#include <cmath>
#include <sscapi.h>

#include "testEvergreen.h"
#include "FunctionsEvergreen.h"

// Radians to degrees conversion factor
static const double R2D = 180 / M_PI;

using namespace std;

CPPUNIT_TEST_SUITE_REGISTRATION(testEvergreen);

testEvergreen::testEvergreen() {
}

testEvergreen::~testEvergreen() {
}

void
testEvergreen::testSoltrack_() {
    
    /*
     * Compare results from soltrack and NOAA sun position calculator
     * 
     * NOAA calculator: https://www.esrl.noaa.gov/gmd/grad/solcalc/azel.html
     * Pay attention that the longitude direction for NOAA calculator is reversed.
     */
    
    ssc_number_t azimuth, elevation, julian_day, lon, lat;
    int year, month, day, hour, minute, second;
    
    // Test 1
    lon = 105; lat = 40;
    year = 2020; month = 3; day = 25;
    hour = 6; minute = 0; second = 0;    
    
    FunctionsEvergreen::sun_pos(azimuth, elevation, julian_day, lon, lat, year, month, day, hour, minute, second);
    cout << "Azimuth: " << azimuth * R2D << " Elevation: " << elevation * R2D << " Cosine of zenith: " << sin(elevation) << endl;
    return;
}
void
testEvergreen::testPvwattsv5_() {

    ssc_number_t ac, dc;

    ssc_data_t pvwatts_data = ssc_data_create();
    ssc_module_t pvwatts_module = ssc_module_create("pvwattsv5_1ts");

    ssc_data_set_number(pvwatts_data, "year", 2018);
    ssc_data_set_number(pvwatts_data, "month", 12);
    ssc_data_set_number(pvwatts_data, "day", 14);
    ssc_data_set_number(pvwatts_data, "hour", 8);
    ssc_data_set_number(pvwatts_data, "minute", 0);

    ssc_data_set_number(pvwatts_data, "lon", -80.9792);
    ssc_data_set_number(pvwatts_data, "lat", 39.0392);
    ssc_data_set_number(pvwatts_data, "tz", -5);

    ssc_data_set_number(pvwatts_data, "beam", 166.99);
    ssc_data_set_number(pvwatts_data, "diffuse", 16.4272);
    ssc_data_set_number(pvwatts_data, "tamb", 4.75);
    ssc_data_set_number(pvwatts_data, "wspd", 1.1404);
    ssc_data_set_number(pvwatts_data, "alb", 0.13);

    ssc_data_set_number(pvwatts_data, "azimuth", 123.716);

    ssc_data_set_number(pvwatts_data, "time_step", 1);
    ssc_data_set_number(pvwatts_data, "system_capacity", 1);
    ssc_data_set_number(pvwatts_data, "module_type", 0);
    ssc_data_set_number(pvwatts_data, "dc_ac_ratio", 1.1);
    ssc_data_set_number(pvwatts_data, "inv_eff", 96);
    ssc_data_set_number(pvwatts_data, "losses", 0);
    ssc_data_set_number(pvwatts_data, "array_type", 0);
    ssc_data_set_number(pvwatts_data, "tilt", 0);
    ssc_data_set_number(pvwatts_data, "gcr", 0.4);
    
    ssc_number_t tcell = NAN, poa = NAN;
    ssc_data_set_number(pvwatts_data, "tcell", tcell);
    ssc_data_set_number(pvwatts_data, "poa", poa);

    CPPUNIT_ASSERT(ssc_module_exec(pvwatts_module, pvwatts_data) != 0);

    ssc_data_get_number(pvwatts_data, "ac", &ac);
    ssc_data_get_number(pvwatts_data, "dc", &dc);
    ssc_data_get_number(pvwatts_data, "tcell", &tcell);
    ssc_data_get_number(pvwatts_data, "poa", &poa);

    cout << "AC: " << ac << " DC: " << dc << " tcell: " << tcell << " poa: " << poa << endl;

    ssc_module_free(pvwatts_module);
    ssc_data_free(pvwatts_data);

    return;
}

