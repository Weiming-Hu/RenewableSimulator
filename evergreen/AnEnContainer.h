/* 
 * File:   AnEnContainer.h
* Author: Weiming Hu <weiming@psu.edu>
 *
 * Created on March 13, 2020, 1:07 PM
 */

#ifndef ANENCONTAINER_H
#define ANENCONTAINER_H

#include <map>
#include <string>
#include <vector>

#include "sscapi.h"
#include "Times.h"
#include "Stations.h"
#include "Array4DPointer.h"
#include "boost/date_time/posix_time/posix_time.hpp"

using ArrayMap = std::map<std::string, Array4DPointer>;

class AnEnContainer {
public:
    AnEnContainer() = delete;
    AnEnContainer(const AnEnContainer& orig) = delete;
    AnEnContainer(const std::string & file_path, const std::vector<std::string> var_names, Verbose verbose);
    virtual ~AnEnContainer() = default;
    
    const Stations & stations() const;
    const Times & times() const;
    const Times & flts() const;
    std::size_t num_analogs() const;
    
    /**
     * Get the time zone in hour offset from GMT for a particular station
     * @param station_index The index of the station
     * @return an hour offset from GMT
     */
    ssc_number_t getTimeZone(std::size_t station_index) const;
    
    /**
     * Get the coordinate of a station
     */
    void getCoordinate(ssc_number_t & lon, ssc_number_t & lat, std::size_t station_index) const;
    
    /**
     * Get the local date and time for a particular forecast time and lead time.
     * The forecast times and lead times are assumed to be in UTC.
     */
    void getLocalDateTime(ssc_number_t & year, ssc_number_t & month, ssc_number_t & day,
            ssc_number_t & hour, ssc_number_t & minute, std::size_t time_index, std::size_t flt_index, int offset) const;
    
    /**
     * Get the UTC date and time for a particular forecast time and lead time.
     * The forecast times and lead times are assumed to be in UTC so this function
     * simply returns the particular forecast time and lead time.
     */
    void getUTCDateTime(ssc_number_t & year, ssc_number_t & month, ssc_number_t & day,
            ssc_number_t & hour, ssc_number_t & minute,
            std::size_t time_index, std::size_t flt_index) const;
    
    /**
     * Get a scalar value from a particular variable.
     */
    ssc_number_t getValue(std::size_t station_index, std::size_t time_index,
            std::size_t flt_index, std::size_t analog_index, std::string var_name) const;

private:
    ArrayMap array_map_;
    Stations stations_;
    Times times_;
    Times flts_;
    
    /**
     * This is a very simple function to calculate the time offset from GMT in hours.
     *
     * (1) It assumes time zones are divided vertically;
     * (2) no daylight saving time is considered.
     * 
     * @param lon A longitude.
     * @return The time offset in hours from GMT
     */
    ssc_number_t simpleOffset_(double lon) const;
    
    void parse_date_time(ssc_number_t & year, ssc_number_t & month, ssc_number_t & day,
            ssc_number_t & hour, ssc_number_t & minute, const boost::posix_time::ptime & ptime) const;
    
};

#endif /* ANENCONTAINER_H */

