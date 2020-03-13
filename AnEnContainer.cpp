/* 
 * File:   AnEnContainer.cpp
* Author: Weiming Hu <weiming@psu.edu>
 * 
 * Created on March 13, 2020, 1:07 PM
 */

#include "AnEnContainer.h"

#include <stdexcept>
#include <numeric>

using namespace std;

static const int _NUM_HOURS = 24;
static const int _NUM_DAYS = 365;
static const int _NUM_TIMES = _NUM_DAYS * _NUM_HOURS;

AnEnContainer::AnEnContainer(const string & file_path) {
    
    

    Times times;
    Stations stations;
    Array4DPointer dn_arr, df_arr, t_arr, wspd_arr;
    AnEnReadNcdf anen_read(verbose);

    // Variable names are hard coded
    NcFile nc(file_path, NcFile::FileMode::read);
    anen_read.read(nc, times, "test_times");
    anen_read.read(nc, stations);

    anen_read.readAnalogs(file_path, dn_arr, "DownwardShortwaveRadiation");
    anen_read.readAnalogs(file_path, df_arr, "UpwardShortwaveRadiation");
    anen_read.readAnalogs(file_path, t_arr, "temperature_2m");
    anen_read.readAnalogs(file_path, wspd_arr, "wspd_1000hPa");

    // Sanity checks
    if (stations.size() != dn_arr.shape()[0]) throw runtime_error("Stations do not match the array shape");
    if (times.size() != dn_arr.shape()[1]) throw runtime_error("Times do not match the array shape");

    // The time span of analogs should not be longer than a year
    double first_time = times.getTime(0).timestamp;
    double last_time = times.getTime(times.size() - 1).timestamp;
    if (last_time - first_time > 365 *)

        /*
         * Input alignment
         * 
         * NREL SAM simulator requires hourly data for the entire year. We need to align
         * input analogs with the year-round hourly time series.
         */
        if (times.size() == _NUM_TIMES) {
            // This is the ideal case. I don't need to align the values
        } else {
            // There are probably missing days in the input analogs.
            // I need to create a 
        }
    
    
    // TODO: Don't forget to shift longitudes
    throw runtime_error("not implemented yet");
}

ssc_number_t
AnEnContainer::getTimezone(size_t station_index) const {
    if (stations_.size() - 1 < station_index) throw range_error("Invalid station index");
    
    double lon = stations_.getStation(station_index).getX();
    
    if (-180 <= lon && lon <= 180) return simpleOffset_(lon);
    else throw runtime_error("Station longitudes need shifting");
}

void
AnEnContainer::subset(map<string, Array4DPointer>& ptr_map_subset,
        size_t& station_i, size_t& flt_i, size_t& analog_i) const {
    
    // Clear the map
    ptr_map_subset.clear();
    
    // These are the indices to extract from the data container
    vector<size_t> stations_i{station_i};
    vector<size_t> flts_i{flt_i};
    vector<size_t> analogs_i{analog_i};
    vector<size_t> times_i(_NUM_TIMES);
    iota(times_i.begin(), times_i.end(), 0);

    // Initialize the map
    for (const auto & pair : ptr_map_) {
        
        // Subset array
        Array4DPointer arr_subset;
        pair.second.subset(stations_i, times_i, flts_i, analogs_i, arr_subset);

        // Push this array to the map
        ptr_map_subset.insert(make_pair(pair.first, arr_subset));
    }
    
    return;
}

void
AnEnContainer::set(const map<string, Array4DPointer> & ptr_map,
            ssc_data_t & data_container, const string & name) {
    
    // Sanity check
    if (ptr_map.size() == 0) throw runtime_error("Input map is empty");

    // Initialize a weather table to populate
    ssc_data_t weather_table = ssc_data_create();

    // Assign data arrays to the table
    ssc_number_t ptr[_NUM_TIMES];
    for (const auto & pair : ptr_map) {
        
        // Sanity check
        if (pair.second.num_elements() != _NUM_TIMES) throw runtime_error("The length of Array4DPointer is incorrect");
        
        // Convert from double to ssc_number_t
        const double * ptr_arr = pair.second.getValuesPtr();
        for (size_t i = 0; i < _NUM_TIMES; ++i) ptr[i] = (ssc_number_t) ptr_arr[i];
        
        // Copy values to the weather table
        ssc_data_set_array(weather_table, pair.first.c_str(), ptr, _NUM_TIMES);
    }

    // Set weather table
    ssc_data_set_table(data_container, name.c_str(), weather_table);

    // Free the table after it is assigned
    ssc_data_free(weather_table);
    return;
}

ssc_number_t
AnEnContainer::simpleOffset_(double lon) const {
    return round(lon * 24 / 360);
}