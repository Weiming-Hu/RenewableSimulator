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


class AnEnContainer {
public:
    AnEnContainer() = delete;
    AnEnContainer(const AnEnContainer& orig) = delete;
    AnEnContainer(const std::string & file_path);
    virtual ~AnEnContainer() = default;
    
    ssc_number_t getTimezone(std::size_t station_index);
    
    void subset(std::map<std::string, Array4DPointer> & ptr_map_subset,
            std::size_t & station_i, std::size_t & flt_i, std::size_t & analog_i) const;
    
    static void set(const std::map<std::string, Array4DPointer> & ptr_map,
            ssc_data_t & data_container);
    
private:
    Stations stations_;
    Times times_;
    Times flts_;
    
    std::map<std::string, Array4DPointer> ptr_map_;
    
    /**
     * This is a very simple function to calculate the time offset from GMT in hours.
     *
     * (1) It assumes time zones are divided vertically;
     * (2) no daylight saving time is considered.
     * 
     * @param lon A longitude.
     * @return The time offset in hours from GMT
     */
    ssc_number_t simpleOffset_(double lon);
    
};

#endif /* ANENCONTAINER_H */

