/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "Atmosphere.hpp"

/**
 * Constructor of the Atmosphere object.
 */
Atmosphere::Atmosphere() = default;

/**
 * Move constructor.
 *
 * @param land  Object to move.
 */
Atmosphere::Atmosphere(Atmosphere &&atmosphere) noexcept {
    // Copy data
    this->_ma = atmosphere._ma;
    this->_e_co2 = atmosphere._e_co2;
    this->_co2_atm0 = atmosphere._co2_atm0;

    // Remove data from the other object
    atmosphere._ma = nullptr;
    atmosphere._e_co2 = nullptr;
    atmosphere._co2_atm0 = 0e0;
}

/**
 * Move operator.
 *
 * @param earth Object to move.
 * @return      Return the current object.
 */
Atmosphere &Atmosphere::operator=(Atmosphere &&atmosphere) noexcept {
    if (this != &atmosphere) {
        // Free existing resources
        delete[] this->_ma;
        delete[] this->_e_co2;

        // Copy data
        this->_ma = atmosphere._ma;
        this->_e_co2 = atmosphere._e_co2;
        this->_co2_atm0 = atmosphere._co2_atm0;

        // Remove data from the other object
        atmosphere._ma = nullptr;
        atmosphere._e_co2 = nullptr;
        atmosphere._co2_atm0 = 0e0;
    }
    return *this;
}

/**
 * Destructor of the Atmosphere object.
 */
Atmosphere::~Atmosphere() {
    delete[] this->_ma;
    delete[] this->_e_co2;
}

/**
 * Initialize the mass of carbon in atmosphere according to the time dimension.
 *
 * @param ntime Time dimension.
 */
void Atmosphere::initialize(int ntime) {
    this->_ma = new double[ntime]{0};
    this->_e_co2 = new double[ntime]{0};
}

/**
 * Initialize the preindustrial equilibrium CO2 concentration.
 */
void Atmosphere::initializeCo2Atm0() {
    this->_co2_atm0 = this->_ma[0] / PPMTOGT;
}

/**
 * Setter for the i-th element of ma.
 *
 * @param i     Index.
 * @param value Value.
 */
void Atmosphere::setMaElem(int i, double value) {
    this->_ma[i] = value;
}

/**
 * Setter for the i-th element of e_co2.
 *
 * @param i     Index.
 * @param value Value.
 */
void Atmosphere::setECo2Elem(int i, double value) {
    this->_e_co2[i] = value;
}

/**
 * Getter for the i-th element of ma (mass of carbon in the atmosphere).
 *
 * @param i Index.
 * @return  Return the i-th value of ma.
 */
double Atmosphere::getMaElem(int i) const {
    return this->_ma[i];
}

/**
 * Getter for the i-th element of eCo2 (CO₂ emissions (GtC/yr)).
 *
 * @param i Index.
 * @return  Return the i-th value of e_co2.
 */
double Atmosphere::getECo2Elem(int i) const {
    return this->_e_co2[i];
}

/**
 * Getter for the preindustrial equilibrium CO2 concentration.
 *
 * @return  Return the co2_atm0 attribute.
 */
double Atmosphere::getCo2Atm0() const {
    return this->_co2_atm0;
}

/**
 * Compute the RF of the atmospheric CO2.
 *
 * @param i Index.
 * @return  Return the RF (radiation forcing) of atmospheric CO2 (Wm⁻²).
 */
double Atmosphere::rfCo2(int i) {
    return RECO2 * UtilsMath::myLog((this->_ma[i] / PPMTOGT) / CO2PREIND);
}
