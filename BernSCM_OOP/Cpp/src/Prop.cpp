/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "Prop.hpp"

/**
 *  Constructor of the propagator.
 *
 * @param nscale    Number of finite timescales.
 * @param x         X value.
 */
Prop::Prop(int nscale, double x) : _nscale(nscale), _x(x) {}

/**
 * Move constructor.
 *
 * @param prop  Propagator to move.
 */
Prop::Prop(Prop &&prop) noexcept {
    this->_nscale = prop._nscale;
    this->_x = prop._x;
    for (int i = 0; i < NSCALEMAX + 1; i++) {
        this->_propm[i] = prop._propm[i];
        this->_propf[i] = prop._propf[i];
        this->_propfo[i] = prop._propfo[i];

        // Remove data from the other object
        prop._propm[i] = 0e0;
        prop._propf[i] = 0e0;
        prop._propfo[i] = 0e0;
    }
    this->_nscale = 0;
    this->_x = 0e0;
}

/**
 * Move operator.
 *
 * @param pirf  Object to move.
 * @return      Return the moved object or the current if not different.
 */
Prop &Prop::operator=(Prop &&prop) noexcept {
    if (this != &prop) {
        // Copy data
        this->_nscale = prop._nscale;
        this->_x = prop._x;
        for (int i = 0; i < NSCALEMAX + 1; i++) {
            this->_propm[i] = prop._propm[i];
            this->_propf[i] = prop._propf[i];
            this->_propfo[i] = prop._propfo[i];

            // Remove data from the other object
            prop._propm[i] = 0e0;
            prop._propf[i] = 0e0;
            prop._propfo[i] = 0e0;
        }
        this->_nscale = 0;
        this->_x = 0e0;
    }
    return *this;
}

/**
 * Destructor of the propagator.
 */
Prop::~Prop() = default;

/**
 * Getter for nscale value.
 *
 * @return  Return the nscale value.
 */
int Prop::getNscale() const {
    return this->_nscale;
}

/**
 * Getter for the i-th element of propm.
 *
 * @param i Index
 * @return  Return the i-th value of propm.
 */
double Prop::getPropmElem(int i) const {
    return this->_propm[i];
}

/**
 * Getter for the i-th element of propf.
 *
 * @param i Index.
 * @return  Return the i-th value of propf.
 */
double Prop::getPropfElem(int i) const {
    return this->_propf[i];
}

// double Prop::getPropfFirstElems(bound) {}

/**
 * Getter for the i-th element of propfo.
 *
 * @param i Index.
 * @return  Return the i-th value of propfo.
 */
double Prop::getPropfoElem(int i) const {
    return this->_propfo[i];
}

/**
 * Calculate propagators.
 *
 * @param pirf  IRF coefficients.
 * @param dt    Time step.
 */
void Prop::computePropagators(Pirf &pirf, double dt) {
    if (!LINEAR) {
        for (int i = 0; i < pirf.getNscale(); i++) {
            if (pirf.getTScaleElem(i) <= TEQUIL) {
                this->_propf[i] = pirf.getWeightElem(i) * pirf.getTScaleElem(i);
                this->_propm[i] = 0e0;
            } else {
                this->_propf[i] =
                        pirf.getWeightElem(i) * pirf.getTScaleElem(i) * (1e0 - exp(-dt / pirf.getTScaleElem(i)));
                this->_propm[i] = exp(-dt / pirf.getTScaleElem(i));
            }
        }
        this->_propf[pirf.getNscale()] = dt * pirf.getWeightElem(pirf.getNscale());
        this->_propm[pirf.getNscale()] = 1e0;
    } else if (LINEAR) {
        for (int i = 0; i < pirf.getNscale(); i++) {
            this->_propf[i] = pirf.getWeightElem(i) * pirf.getTScaleElem(i) *
                              (dt - pirf.getTScaleElem(i) * (1e0 - exp(-dt / pirf.getTScaleElem(i)))) / dt;
            this->_propfo[i] =
                    pirf.getWeightElem(i) * pirf.getTScaleElem(i) * (1e0 - exp(-dt / pirf.getTScaleElem(i))) -
                    this->_propf[i];
        }
        this->_propf[pirf.getNscale()] = pirf.getWeightElem(pirf.getNscale()) * dt / 2e0;
        this->_propfo[pirf.getNscale()] = pirf.getWeightElem(pirf.getNscale()) * dt - this->_propf[pirf.getNscale()];

        for (int i = 0; i < pirf.getNscale(); i++) {
            this->_propm[i] = exp(-dt / pirf.getTScaleElem(i));
        }
        this->_propm[pirf.getNscale()] = 1e0;
    }
    this->_nscale = pirf.getNscale();
}
