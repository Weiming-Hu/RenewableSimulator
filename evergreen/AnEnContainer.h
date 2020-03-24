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

#include "sscapi.h"

#include "Times.h"
#include "Stations.h"
#include "Array4DPointer.h"

using ArrayMap = std::map<std::string, Array4DPointer>;

class AnEnContainer {
public:
    AnEnContainer() = delete;
    AnEnContainer(const AnEnContainer& orig) = delete;
    AnEnContainer(const std::string & file_path);
    virtual ~AnEnContainer() = default;
    
    const Stations & stations() const;
    const Times & times() const;
    const Times & flts() const;
    std::size_t num_analogs() const;
    
    ssc_number_t getTimeZone(std::size_t station_index) const;
    void getCoordinate(ssc_number_t & lon, ssc_number_t & lat, std::size_t station_index) const {};
    void getLocalDateTime(ssc_number_t & year, ssc_number_t & month, ssc_number_t & day,
            ssc_number_t & hour, ssc_number_t & minute,
            std::size_t time_index, std::size_t flt_index) const {};
    void getUTCDateTime(ssc_number_t & year, ssc_number_t & month, ssc_number_t & day,
            ssc_number_t & hour, ssc_number_t & minute,
            std::size_t time_index, std::size_t flt_index) const {};
    ssc_number_t getValue(std::size_t station_index, std::size_t time_index,
            std::size_t flt_index, std::size_t analog_index, std::string var_name) const {return 0;};

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
    
};

#endif /* ANENCONTAINER_H */

