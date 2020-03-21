/* 
 * File:   Scenarios.cpp
 * Author: Weiming Hu <weiming@psu.edu>
 * 
 * Created on March 13, 2020, 1:08 PM
 */

#include "Scenarios.h"
#include "Functions.h"

#include <cmath>

using namespace std;

/******************************************************************************
 *                                    Scenario                                *
 ******************************************************************************/

void
Scenario::set(ssc_data_t & data_container) const {

    for (const auto & e : * this) {
        ssc_data_set_number(data_container, e.first.c_str(), e.second);
    }

    return;
}


void
Scenario::print(ostream & os) const {
    os << "{ ";
    for (const auto & e : *this) {
        os << e.first << ":" << e.second << " ";
    }
    os << "}";
    return;
}

ostream &
operator<<(ostream & os, const Scenario & scenario) {
    scenario.print(os);
    return os;
}


/******************************************************************************
 *                                   Scenarios                                *
 ******************************************************************************/


size_t
Scenarios::totalScenarios() const {

    if (size() == 0) return 0;

    size_t total = 1;
    for (const auto & e : *this) total *= e.second.size();
    return total;
}

void
Scenarios::set(ssc_data_t & data_container, size_t index) const {
    Scenario scenario;
    get(scenario, index);
    scenario.set(data_container);
    return;
}

void
Scenarios::get(Scenario & scenario, size_t index) const {

    // Clear map
    scenario.clear();

    // Initialize values
    size_t accum_len = 1;
    size_t current_index;

    for (const auto & e : * this) {

        // Calculate the current index for this key based on the global index
        if (accum_len == 1) {
            accum_len = e.second.size();
            current_index = fmod(index, accum_len);
        } else {
            current_index = fmod(index / accum_len, e.second.size());
            accum_len *= e.second.size();
        }

        // Insert pair
        scenario.insert(make_pair(e.first, e.second.at(current_index)));
    }

    return;
}

void
Scenarios::print(ostream & os) const {
    os << "Scenarios:" << endl;
    for (const auto & e : *this) {
        os << e.first << ":" << Functions::format(e.second, ",") << endl;
    }
    return;
}

void
Scenarios::printMore(ostream & os) const {
    os << "A list of available scenarios:" << endl;
    
    size_t total = totalScenarios();
    for (size_t i = 0 ; i < total; ++i) {
        Scenario scenario;
        get(scenario, i);
        os << "Scenario #" << i << "/" << total << ": " << scenario << endl;
    }
    
    return;
}

ostream &
operator<<(ostream & os, const Scenarios & scenarios) {
    scenarios.print(os);
    return os;
}
