/* 
 * File:   FunctionsEvergreen.cpp
 * Author: Weiming Hu <weiming@psu.edu>
 * 
 * Created on March 23, 2020, 2:52 PM
 */

#include <cmath>
#include "FunctionsEvergreen.h"

extern "C" {
#include "SolTrack.h"
}

// Solar constant in W/m^2
static const double I_sc = 1367;

void
FunctionsEvergreen::sun_pos(ssc_number_t & azimuth, ssc_number_t & elevation, ssc_number_t & julian_day,
        ssc_number_t lon, ssc_number_t lat, int year, int month, int day, int hour, int minute, double second) {

    /*
     * Use SolTrack to calculate the sun position
     */

    static const int use_degrees = 0;
    static const int use_north_equals_zero = 1;
    static const int compute_refr_equatorial = 0;
    static const int compute_distance = 0;

    // Note that date time should be in UTC
    Time time = {year, month, day, hour, minute, second};
    Location location = {lon / R2D, lat / R2D};

    // Calculate sun position
    Position position;
    SolTrack(time, location, &position, use_degrees,
            use_north_equals_zero, compute_refr_equatorial, compute_distance);

    elevation = position.altitude;
    azimuth = position.azimuthRefract;
    julian_day = computeJulianDay(year, month, day, hour, minute, second);
    return;
}

void
FunctionsEvergreen::decompose_ghi(ssc_number_t ghi, ssc_number_t& dhi, ssc_number_t& dni, ssc_number_t elevation, ssc_number_t julian_day) {
    
    
    /*
     * Reindl et al. Decomposition Model 1 to estimate DHI
     * 
     * Reference: https://pvpmc.sandia.gov/modeling-steps/1-weather-design-inputs/irradiance-and-insolation-2/direct-normal-irradiance/piecewise_decomp-models/
     */

    // Calculate the cosine of zenith angle (sine of elevation angle)
    double cos_zenith = sin(elevation);
    
    // Deal with darkness
    if (cos_zenith <= 0) {
        dhi = NAN;
        dni = NAN;
        return;
    }

    // Correct solar constant for Earth's elliptical orbit
    // 
    // Reference: https://www.itacanet.org/the-sun-as-a-source-of-energy/part-4-irradiation-calculations/
    //
    double I = I_sc * (1 + 0.034 * cos(2 * PI * julian_day / 365.25));
    double k_t = ghi / (I * cos_zenith);

    // Piecewise linear regression from Reindl et al. with constraints
    double k_d;
    if (k_t <= 0.3) {
        k_d = 1.02 - 0.248 * k_t;
        if (k_d > 1) k_d = 1;
    } else if (0.3 < k_t < 0.78) {
        k_d = 1.45 - 1.67 * k_t;
    } else {
        k_d = 0.147;
    }

    dhi = k_d * ghi;
    
    
    /*
     * Calculate DNI from GHI and DHI
     * 
     * Reference: https://pvpmc.sandia.gov/modeling-steps/1-weather-design-inputs/irradiance-and-insolation-2/global-horizontal-irradiance/
     */
    dni = (ghi - dhi) / cos_zenith;
    
    // Take care of negativity
    if (dni < 0) dni = 0;
    
    return;
}
