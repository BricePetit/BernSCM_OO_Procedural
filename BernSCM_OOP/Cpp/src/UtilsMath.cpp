/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "UtilsMath.hpp"

/**
 * Default constructor.
 */
UtilsMath::UtilsMath() = default;

/**
 * Default destructor.
 */
UtilsMath::~UtilsMath() = default;

/**
 * We create our specific log because the owner use a certain flag in order to do the following.
 *
 * @param value The value to use in the logarithm.
 * @return      Return the computed value.
 */
double UtilsMath::myLog(const double &value) {
    double res;
    if (value > 0e0) {
        res = log(value);
    } else {
        res = value;
    }
    return res;
}
