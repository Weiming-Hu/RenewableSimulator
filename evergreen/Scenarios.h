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
#include <iostream>

#include "sscapi.h"

class Scenario : public std::map<std::string, ssc_number_t> {    
public:
    Scenario() = default;
    Scenario(const Scenario& orig) = default;
    virtual ~Scenario() = default;
    
    /**
     * Set data container using the current scenario
     * @param data_container
     */
    void set(ssc_data_t & data_container) const;
    
    virtual void print(std::ostream &) const;
    friend std::ostream & operator<<(std::ostream &, const Scenario &);
};

class Scenarios : public std::map<std::string, std::vector<ssc_number_t> > {
public:
    Scenarios() = default;
    Scenarios(const Scenarios& orig) = default;
    virtual ~Scenarios() = default;

    /**
     * Count the total number of available scenarios defined
     */
    size_t totalScenarios() const;    
    
    /**
     * Set the data container using the scenario corresponding to the scenario index.
     * Usually, it is preferred to use the function totalScenarios to check how many
     * scenarios are defined within this object.
     *
     * @param data_container
     * @param index
     */
    void set(ssc_data_t & data_container, std::size_t index) const;
    
    /**
     * Extract a certain scenario from the predefined scenarios. The scenario is
     * a combination of different configurations in the predefined scenarios. The
     * index is the scenario index. You should use countScenarios first to see
     * the total number of scenarios available to extract.
     * 
     * @param scenario A Scenario. It is essentially a string-number STL map
     * @param index The scenario ID to extract. You can use totalScenarios
     * to get the total number of scenarios available.
     */
    void get(Scenario & scenario, std::size_t index) const;
    
    virtual void print(std::ostream &) const;
    virtual void printMore(std::ostream &) const;
    friend std::ostream & operator<<(std::ostream &, const Scenarios &);
};

#endif /* SCENARIOS_H */

