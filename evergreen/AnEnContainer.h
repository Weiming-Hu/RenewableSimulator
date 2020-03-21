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
    
    /**
     * Sets the lat/lon/timezone attributes of in the data container using
     * the location from the specified station.
     * @param data_container The SSC data container to set
     * @param station_index The station index
     */
    void setStation(ssc_data_t & data_container, std::size_t station_index) const;
    
    
    
    
    
    ssc_number_t getTimezone(std::size_t station_index) const;

    void subset(ArrayMap & ptr_map_subset, std::size_t & station_i, std::size_t & flt_i, std::size_t & analog_i) const;
    
    static void set(const ArrayMap & ptr_map, ssc_data_t & data_container, const std::string & name = "solar_resource_data");
    
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

