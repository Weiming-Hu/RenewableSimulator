/* 
 * File:   AnEnContainer.cpp
* Author: Weiming Hu <weiming@psu.edu>
 * 
 * Created on March 13, 2020, 1:07 PM
 */

#include "AnEnContainer.h"
#include "AnEnReadNcdf.h"
#include "netcdf"

#include <stdexcept>
#include <numeric>

using namespace std;

static const int _NUM_HOURS = 24;
static const int _NUM_DAYS = 365;
static const int _NUM_TIMES = _NUM_DAYS * _NUM_HOURS;

AnEnContainer::AnEnContainer(const string & file_path, vector<string> var_names, Verbose verbose) {
    
    AnEnReadNcdf anen_read(verbose);
    Array4DPointer array;
    
    for (const auto & var_name : var_names) {
        anen_read.readAnalogs(file_path, array, var_name);
        array_map_.insert(make_pair(var_name, array));
    }
    
    netCDF::NcFile nc(file_path, netCDF::NcFile::FileMode::read);
    anen_read.read(nc, stations_);
    anen_read.read(nc, times_, Config::_TEST_TIMES);
    anen_read.read(nc, flts_, Config::_FLTS);

    // Sanity checks
    if (stations_.size() != array_map_.begin()->second.shape()[0]) throw runtime_error("Stations do not match the array shape");
    if (times_.size() != array_map_.begin()->second.shape()[1]) throw runtime_error("Times do not match the array shape");

    // Shift coordinates
    vector<double> xs, ys;
    stations_.getCoordinates(xs, ys);
    
    if (*max_element(xs.begin(), xs.end()) > 180) stations_.shiftLongitudes();
    return;
}


const Stations &
AnEnContainer::stations() const {
    return stations_;
}

const Times &
AnEnContainer::times() const {
    return times_;
}

const Times &
AnEnContainer::flts() const {
    return flts_;
}

size_t
AnEnContainer::num_analogs() const {
    if (array_map_.empty()) throw runtime_error("Array map is empty. No analog members are found.");
    return array_map_.begin()->second.shape()[3];
}

ssc_number_t
AnEnContainer::getTimeZone(size_t station_index) const {
    if (stations_.size() - 1 < station_index) throw range_error("Invalid station index");
    
    double lon = stations_.getStation(station_index).getX();
    
    if (-180 <= lon && lon <= 180) return simpleOffset_(lon);
    else throw runtime_error("Station longitudes need shifting");
}

void
AnEnContainer::getCoordinate(ssc_number_t & lon, ssc_number_t & lat, size_t station_index) const {
    const auto & station = stations_.getStation(station_index);
    lon = station.getX();
    lat = station.getY();
    return;
}

void
AnEnContainer::getLocalDateTime(ssc_number_t & year, ssc_number_t & month, ssc_number_t & day,
        ssc_number_t & hour, ssc_number_t & minute, size_t time_index, size_t flt_index, int offset) const {
    
    // Convert the selected time to a formatted POSIX time
    Time selected_time = times_.getTime(time_index) + flts_.getTime(flt_index);
    boost::posix_time::ptime current_ptime = boost::posix_time::time_from_string(selected_time.toString());
    
    // Convert from UTC to local time
    current_ptime = current_ptime + boost::posix_time::hours(offset);
    
    // Get components
    parse_date_time(year, month, day, hour, minute, current_ptime);
    return;
}

void
AnEnContainer::getUTCDateTime(ssc_number_t & year, ssc_number_t & month, ssc_number_t & day,
        ssc_number_t & hour, ssc_number_t & minute, size_t time_index, size_t flt_index) const {
    
    // Convert the selected time to a formatted POSIX time
    Time selected_time = times_.getTime(time_index) + flts_.getTime(flt_index);
    boost::posix_time::ptime current_ptime = boost::posix_time::time_from_string(selected_time.toString());
    
    // Get components of the POSIX time
    // Note that it is assumed that the weather time is always in UTC, and therefore
    // no conversion is done.
    //
    parse_date_time(year, month, day, hour, minute, current_ptime);
    return;
}

ssc_number_t
AnEnContainer::getValue(size_t station_index, size_t time_index, size_t flt_index, size_t analog_index, string var_name) const {
    const auto pair = array_map_.find(var_name);
    return pair->second.getValue(station_index, time_index, flt_index, analog_index);
}

ssc_number_t
AnEnContainer::simpleOffset_(double lon) const {
    return round(lon * 24 / 360);
}

void
AnEnContainer::parse_date_time(ssc_number_t& year, ssc_number_t& month, ssc_number_t& day,
        ssc_number_t& hour, ssc_number_t& minute, const boost::posix_time::ptime & ptime) const {
    
    const auto & current_date = ptime.date();
    year = current_date.year();
    month = current_date.month();
    day = current_date.day();

    const auto & current_time = ptime.time_of_day();
    hour = current_time.hours();
    minute = current_time.minutes();
    
    return;
}
