/* 
 * File:   FunctionsEvergreen.h
 * Author: Weiming Hu <weiming@psu.edu>
 *
 * Created on March 23, 2020, 2:52 PM
 */

#ifndef FUNCTIONSEVERGREEN_H
#define FUNCTIONSEVERGREEN_H

#include "sscapi.h"
#include "string"

namespace FunctionsEvergreen {
    
    using namespace std;
    
    void sun_pos(ssc_number_t & azimuth_rad, ssc_number_t & elevation_rad, ssc_number_t & julian_day,
            ssc_number_t lon, ssc_number_t lat, int UTC_year, int UTC_month, int UTC_day, int UTC_hour, int UTC_minute, double UTC_second);
    
    void decompose_ghi(ssc_number_t ghi, ssc_number_t & dhi, ssc_number_t & dni, ssc_number_t elevation, ssc_number_t julian_day);

    string toString(const ssc_data_t &, const string & module_name);
};

#endif /* FUNCTIONSEVERGREEN_H */

