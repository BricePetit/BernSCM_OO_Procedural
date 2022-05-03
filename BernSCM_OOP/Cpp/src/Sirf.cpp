/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "Sirf.hpp"

/**
 * Default Sirf.
 */
Sirf::Sirf() = default;

/**
 * Constructor of the SIRF object.
 *
 * @param t_max     Parametrization range.
 * @param weight    Temperature sensitivity of weights.
 * @param t_scale   Temperature sensitivity of tscales.
 */
Sirf::Sirf(double t_max, double weight[], double t_scale[]) : _t_max(t_max), _weight(weight), _t_scale(t_scale) {}

/**
 * Copy constructor.
 *
 * @param sirf  Object to copy.
 */
Sirf::Sirf(const Sirf &sirf) {
    this->_t_max = sirf._t_max;
    this->_weight = new double[WEIGHT_SIZE];
    memcpy(this->_weight, sirf._weight, sizeof(double) * WEIGHT_SIZE);
    this->_t_scale = new double[T_SCALE_SIZE];
    memcpy(this->_t_scale, sirf._t_scale, sizeof(double) * T_SCALE_SIZE);

}

/**
 * Move constructor.
 *
 * @param sirf  Object to move.
 */
Sirf::Sirf(Sirf &&sirf) noexcept {
    // Copy data.
    this->_t_max = sirf._t_max;
    this->_weight = sirf._weight;
    this->_t_scale = sirf._t_scale;

    // Remove data from the object to copy.
    sirf._t_max = 0e0;
    sirf._weight = nullptr;
    sirf._t_scale = nullptr;
}

/**
 * Move operator.
 *
 * @param sirf  Object to move.
 * @return      Return the current object.
 */
Sirf &Sirf::operator=(Sirf &&sirf) noexcept {
    if (this != &sirf) {
        // Free existing resources
        delete[] this->_weight;
        delete[] this->_t_scale;
        // Copy data
        this->_t_max = sirf._t_max;
        this->_weight = sirf._weight;
        this->_t_scale = sirf._t_scale;
        // Remove data from the other object
        sirf._t_max = 0e0;
        sirf._weight = nullptr;
        sirf._t_scale = nullptr;
    }
    return *this;
}

/**
 * Destructor of the SIRF object.
 */
//Sirf::~Sirf() {
//    delete[] this->_weight;
//    delete[] this->_t_scale;
//}
Sirf::~Sirf() = default;
/**
 * Return the t_max (Temperature max).
 *
 * @return Return the t_max value.
 */
double Sirf::getTMax() const {
    return this->_t_max;
}

/**
 * Getter for the i-th element of the weight array.
 *
 * @param i Index.
 * @return  The i-th value of the array of weights.
 */
double Sirf::getWeightElem(int i) const {
    return this->_weight[i];
}

/**
 * Getter for the i-th element of the t_scale array.
 *
 * @param i Index.
 * @return  The i-th value of the array of t_scale.
 */
double Sirf::getTscaleElem(int i) const {
    return this->_t_scale[i];
}
