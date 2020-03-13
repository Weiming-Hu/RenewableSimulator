/* 
 * File:   Scenarios.h
 * Author: Weiming Hu <weiming@psu.edu>
 *
 * Created on March 13, 2020, 1:08 PM
 */

#ifndef SCENARIOS_H
#define SCENARIOS_H

#include <map>
#include <vector>
#include <string>

#include "sscapi.h"

class Scenario : public std::map<std::string, ssc_number_t> {    
public:
    Scenario() = default;
    Scenario(const Scenario& orig) = default;
    virtual ~Scenario() = default;
    
    void set(ssc_data_t & data_container) const;
    
private:

};

class Scenarios : public std::map<std::string, std::vector<ssc_number_t> > {
public:
    Scenarios() = default;
    Scenarios(const Scenarios& orig) = default;
    virtual ~Scenarios() = default;

    void totalScenarios() const;    
    void set(ssc_data_t & data_container, std::size_t index) const;

private:

};

#endif /* SCENARIOS_H */

