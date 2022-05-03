/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "Pirf.hpp"

Pirf::Pirf() = default;

/**
 * This class represents IRF parameter set.
 *
 * @param name      Name of the IRF.
 * @param nscale    Number of finite timescales.
 * @param weight    Array of weight or fraction of input for each box.
 * @param t_scale   Array of time scales of each box (yr).
 */
Pirf::Pirf(string_view name, int nscale, double weight[], double t_scale[]) : _name(name), _nscale(nscale),
                                                                              _weight(weight), _t_scale(t_scale) {}

/**
 * Copy constructor with deepcopy.
 *
 * @param p Pirf object to copy.
 */
Pirf::Pirf(const Pirf &p) noexcept {
    this->_name = p._name;
    this->_nscale = p._nscale;
    this->_weight = new double[WEIGHT_SIZE];
    memcpy(this->_weight, p._weight, sizeof(double) * WEIGHT_SIZE);
    this->_t_scale = new double[T_SCALE_SIZE];
    memcpy(this->_t_scale, p._t_scale, sizeof(double) * T_SCALE_SIZE);
}

/**
 * Move constructor
 *
 * @param pirf Object to move in the current object.
 */
Pirf::Pirf(Pirf &&pirf) noexcept {
    // Copy data.
    this->_name = pirf._name;
    this->_nscale = pirf._nscale;
    this->_weight = pirf._weight;
    this->_t_scale = pirf._t_scale;

    // Remove data from the object to copy.
    pirf._name = "";
    pirf._nscale = 0;
    pirf._weight = nullptr;
    pirf._t_scale = nullptr;
}

/**
 * Move operator.
 *
 * @param pirf  Object to move.
 * @return      Return the current object.
 */
Pirf &Pirf::operator=(Pirf &&pirf) noexcept {
    if (this != &pirf) {
        // Free existing resources
        delete[] this->_weight;
        delete[] this->_t_scale;
        // Copy data
        this->_name = pirf._name;
        this->_nscale = pirf._nscale;
        this->_weight = pirf._weight;
        this->_t_scale = pirf._t_scale;
        // Remove data from the other object
        pirf._weight = nullptr;
        pirf._t_scale = nullptr;
        pirf._name = "";
        pirf._nscale = 0;
    }
    return *this;
}

/**
* Destructor of the PIRF object.
*/
Pirf::~Pirf() = default;

/**
 * Getter for the name of the PIRF.
 *
 * @return  Return the name of the PIRF.
 */
string Pirf::getName() const {
    return this->_name;
}

/**
 * Getter for the nscale value.
 *
 * @return  Return nscale.
 */
int Pirf::getNscale() const {
    return this->_nscale;
}

/**
 * Getter for the i-th element of weight.
 *
 * @param index The element i of the weight.
 * @return      Return the i-th weight value.
 */
double Pirf::getWeightElem(int i) const {
    return this->_weight[i];
}

/**
 * Getter for the i-th element of t_scale.
 *
 * @param i The element i of the t_scale.
 * @return  Return the i-th t_scale value.
 */
double Pirf::getTScaleElem(int i) const {
    return this->_t_scale[i];
}

/**
 * Multiply weights in the list starting from the first element up to the bound by the values.
 *
 * @param bound The index used to define the number of value that we want to multiply.
 * @param value Array of values that we want to use to multiply the weight.
 */
void Pirf::multiplyWeights(int bound, const double *values) {
    for (int i = 0; i < bound; i++) {
        this->_weight[i] *= values[i];
    }
}

/**
 * Divide weights in the list starting from the first element up to the bound by the value.
 *
 * @param bound Upper bound of the list.
 * @param value Value by which we will divide weights.
 */
void Pirf::divideWeights(int bound, double value) {
    for (int i = 0; i < bound; i++) {
        this->_weight[i] /= value;
    }
}

/**
 * Multiply t_scale in the list starting from the first element up to the index i by the values.
 *
 * @param bound Upper bound of the list.
 * @param value Array of values by which we will multiply t_scale values.
 */
void Pirf::multiplyTscales(int bound, const double *values) {
    for (int i = 0; i < bound; i++) {
        this->_t_scale[i] *= values[i];
    }
}
