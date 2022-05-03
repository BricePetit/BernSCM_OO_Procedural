/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "Land.hpp"

#include <utility>

/**
 * Default constructor.
 */
Land::Land() = default;

/**
 * Constructor of the Land object.
 *
 * @param pirf          IRF coefficients.
 * @param sirf          IRF sensitivities. Can be none according to the model.
 * @param npp0          Preindustrial (potential natural) NPP (net primary production) (Gt/yr).
 * @param fert          Fertilization Coefficient, only in the npp_4box. Can be none according to the model.
 * @param model_name    The name of the model.
 * @param doc           The documentation about the model.
 */
Land::Land(Pirf pirf, Sirf sirf, double npp0, double fert, string_view model_name, string doc) : _pirf(move(pirf)),
                                                                                                   _sirf(move(sirf)),
                                                                                                   _npp0(npp0),
                                                                                                   _fert(fert),
                                                                                                   _model_name(
                                                                                                           model_name),
                                                                                                   _doc(move(doc)) {}

/**
 * Constructor of the Land object. This constructor does not have pirf because not needed for the model.
 *
 * @param pirf          IRF coefficients.
 * @param npp0          Preindustrial (potential natural) NPP (net primary production) (Gt/yr).
 * @param fert          Fertilization Coefficient, only in the npp_4box. Can be none according to the model.
 * @param model_name    The name of the model.
 * @param doc           The documentation about the model.
*/
Land::Land(Pirf pirf, double npp0, double fert, string_view model_name, string doc) : _pirf(move(pirf)), _npp0(npp0),
                                                                                       _fert(fert),
                                                                                       _model_name(model_name),
                                                                                       _doc(move(doc)) {}

/**
 * Move constructor.
 *
 * @param land  Object to move.
 */
Land::Land(Land &&land) noexcept {
    this->_pirf = move(land._pirf);
    this->_sirf = move(land._sirf);
    this->_npp0 = land._npp0;
    this->_fert = land._fert;
    this->_model_name = land._model_name;
    this->_doc = move(land._doc);
    this->_prop = move(land._prop);
    this->_mlk = land._mlk;
    this->_ml = land._ml;
    this->_exceed = land._exceed;
    this->_co2_exceed = land._co2_exceed;
    this->_t_exceed = land._t_exceed;

    // Delete value for the move.
    land._pirf = Pirf("", 0, {}, {});
    land._sirf = Sirf();
    land._npp0 = 0e0;
    land._fert = 0e0;
    land._model_name = "";
    land._doc = "";
    land._mlk = nullptr;
    land._ml = nullptr;
    land._exceed = false;
    land._co2_exceed = false;
    land._t_exceed = false;
}

/**
 * Move operator.
 *
 * @param land  Object to move.
 * @return      Return the current object.
 */
Land &Land::operator=(Land &&land) noexcept {
    if (this != &land) {
        // Free existing resources
        delete[] this->_mlk;
        delete[] this->_ml;
        // Copy data
        this->_pirf = move(land._pirf);
        this->_sirf = move(land._sirf);
        this->_npp0 = land._npp0;
        this->_fert = land._fert;
        this->_model_name = land._model_name;
        this->_prop = move(land._prop);
        this->_doc = move(land._doc);
        this->_mlk = land._mlk;
        this->_ml = land._ml;
        this->_exceed = land._exceed;
        this->_co2_exceed = land._co2_exceed;
        this->_t_exceed = land._t_exceed;
        // Remove data from the other object
        land._pirf = Pirf();
        land._sirf = Sirf();
        land._npp0 = 0e0;
        land._fert = 0e0;
        land._model_name = "";
        land._prop = Prop(0, 0e0);
        land._doc = "";
        land._mlk = nullptr;
        land._ml = nullptr;
        land._exceed = false;
        land._co2_exceed = false;
        land._t_exceed = false;
    }
    return *this;
}

/**
 * Destructor of the Land object.
 */
Land::~Land() {
    delete[] this->_mlk;
    delete[] this->_ml;
}

/**
 * Initialize the stock of carbon.
 *
 * @param ntime The time dimension.
 */
void Land::initCarbonStock(int ntime) {
    this->_mlk = new double[this->_pirf.getNscale() + 1]{0};
    this->_ml = new double[ntime]{0};
    for (int i = 0; i < this->_pirf.getNscale(); i++) {
        this->_mlk[i] = (this->_npp0 * this->_pirf.getWeightElem(i) * this->_pirf.getTScaleElem(i));
        this->_ml[0] += this->_mlk[i];
    }
}

/**
 * Increase mlk at the index by value.
 *
 * @param index Index.
 * @param value New values to add at the mlk array.
 */
void Land::increaseMlkElem(int index, double value) {
    this->_mlk[index] += value;
}

/**
 * Setter for the i-th element of ml.
 *
 * @param i     Index.
 * @param value Value.
 */
void Land::setMlElem(int i, double value) {
    this->_ml[i] = value;
}

/**
 * As the co2 exceeded, we turn the variable to True.
 */
void Land::co2Exceed() {
    this->_co2_exceed = true;
}

/**
 * As the temperature exceeded, we turn the variable to True.
 */
void Land::tExceed() {
    this->_t_exceed = true;
}

/**
 * Getter for the npp0 attribute.
 *
 * @return  Return preindustrial net primary production (Gt/yr).
 */
double Land::getNpp0() const {
    return this->_npp0;
}

/**
 * Getter for the documentation of the model.
 *
 * @return  Return the documentation.
 */
string Land::getDoc() const {
    return this->_doc;
}

/**
 * Getter for the propagator.
 *
 * @return  Return the propagator.
 */
Prop &Land::getPropagator() {
    return this->_prop;
}

/**
 * Getter for the nscale of the propagator.
 *
 * @return  Return the nscale of the propagator
 */
int Land::getPropNscale() const {
    return this->_prop.getNscale();
}

/**
 * Getter for the i-th element of propf attribute from the propagator.
 *
 * @param i Index.
 * @return  Return the i-th value of propf.
 */
double Land::getPropfElem(int i) const {
    return this->_prop.getPropfElem(i);
}

/**
 * Getter for the mlk attribute (Land biosphere C pools (GtC)).
 *
 * @return  Return the mlk array.
 */
double *Land::getMlk() {
    return this->_mlk;
}

/**
 * Getter for the i-th element of ml.
 *
 * @param i Index.
 * @return  Return the i-th value of ml.
 */
double Land::getMlElem(int i) const {
    return this->_ml[i];
}

/**
 * Check if the parametrization range of IRF exceed.
 *
 * @return  Check if it exceeded.
 */
bool Land::isExceed() const {
    return this->_exceed;
}

/**
 * Check if the co2 exceed.
 *
 * @return  Check if it exceeded.
 */
bool Land::isCo2Exceed() const {
    return this->_co2_exceed;
}

/**
 * Check if the temperature exceed.
 *
 * @return  Check if it exceeded.
 */
bool Land::isTExceed() const {
    return this->_t_exceed;
}

/**
 * Calculate propagators.
 *
 * @param dt    Time step.
 * @param pirf  IRF coefficients.
 */
void Land::computePropagators(double dt, Pirf *pirf) {
    if (pirf == nullptr) {
        this->_prop.computePropagators(this->_pirf, dt);
    } else {
        this->_prop.computePropagators(*pirf, dt);
    }
}

/**
 * Wrapper that updates propagators for T-dependent IRF coefficients.
 *
 * @param t     Temperature perturbation.
 * @param dt    Time step.
 */
void Land::computePropTempDependent(double t, double dt) {
    Pirf p = Pirf(this->_pirf);
    int i;
    double tmp_weights[p.getNscale() + 1], tmp_t_scales[p.getNscale()], tmp_sum = 0;
    if (t > this->_sirf.getTMax() && !this->_exceed) {
        cerr << "warning: temperature parametrization range of IRF " << p.getName() << " exceeded ("
             << this->_sirf.getTMax() << "K)." << endl;
        this->_exceed = true;
    }
    for (i = 0; i < p.getNscale() + 1; i++) {
        tmp_weights[i] = exp(this->_sirf.getWeightElem(i) * t);
    }
    p.multiplyWeights(p.getNscale() + 1, tmp_weights);
    for (i = 0; i < p.getNscale() + 1; i++) {
        tmp_sum += p.getWeightElem(i);
    }
    p.divideWeights(p.getNscale() + 1, tmp_sum);
    for (i = 0; i < p.getNscale(); i++) {
        tmp_t_scales[i] = exp(-this->_sirf.getTscaleElem(i) * t);
    }
    p.multiplyTscales(p.getNscale(), tmp_t_scales);
    this->computePropagators(dt, &p);
}

/**
 * Do the computation of the npp according to the model.
 *
 * @param ma        Mass of C in atmosphere.
 * @param t         Global ΔSAT (℃)
 * @param co2_atm0  Preindustrial equilibrium CO2 concentration.
 * @param co2_dep   Co2 dependence.
 * @param t_dep     Temperature dependence.
 * @param deriv     Derivative dNPP/dm.
 * @return          Return the npp.
 */
double Land::computeNpp(double &ma, double &t, double const &co2_atm0, bool &co2_dep, bool &t_dep, bool &deriv) {
    double res = 0e0;
    if (this->_model_name == "HRBM") {
        res = this->computeHrbmNpp(ma, t, co2_dep, t_dep, deriv);
    } else if (this->_model_name == "4box") {
        res = this->compute4boxNpp(ma, co2_atm0, co2_dep, deriv);
    }
    return res;
}

/**
 * Compute the net primary production for the model HRBM.
 *
 * @param ma        Mass of C in atmosphere.
 * @param t         Global ΔSAT (℃)
 * @param co2_dep   Co2 dependence.
 * @param t_dep     Temperature dependence.
 * @param deriv     Derivative dNPP/dm.
 * @return          Return the npp.
 */
double Land::computeHrbmNpp(double &ma, double &t, bool &co2_dep, bool &t_dep, bool &deriv) {
    double npp;
    if (ma / PPMTOGT > 1120 && !this->_co2_exceed) {
        cerr << "warning: CO2 parametrization range for NPP in HRBM exceeded (1120ppm)" << endl;
        this->_co2_exceed = true;
    }
    if (t > 5 && !this->_t_exceed) {
        cerr << "warning: temperature parametrization range for NPP in HRBM exceeded (5K)" << endl;
        this->_t_exceed = true;
    }
    if (!deriv) {
        if (co2_dep) {
            npp = -exp(3.672801e0) + exp(-0.430818e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 1) -
                  exp(-6.145559e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 2) +
                  exp(-12.353878e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 3) -
                  exp(-19.010800e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 4) +
                  exp(-26.183752e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 5) -
                  exp(-34.317488e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 6) -
                  exp(-41.553715e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 7) +
                  exp(-48.265138e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 8) -
                  exp(-56.056095e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 9) +
                  exp(-64.818185e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 10);
        } else {
            npp = this->_npp0;
        }
        if (t_dep) {
            npp *= (1 + 0.11780208e+0 * tanh(t / 0.509312421e+02) + 0.24305130e-02 * tanh(t / 0.885326739e+01));
        }
    } else {
        double npp_dev;
        if (co2_dep) {
            npp_dev = (+exp(-0.430818e0) - 2 * exp(-6.145559e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 1) +
                       3 * exp(-12.353878e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 2) -
                       4 * exp(-19.010800e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 3) +
                       5 * exp(-26.183752e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 4) -
                       6 * exp(-34.317488e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 5) -
                       7 * exp(-41.553715e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 6) +
                       8 * exp(-48.265138e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 7) -
                       9 * exp(-56.056095e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 8) +
                       10 * exp(-64.818185e0 + UtilsMath::myLog(min(ma / PPMTOGT, 1274e0)) * 9)) / PPMTOGT;
        } else {
            npp_dev = 0e0;
        }
        if (t_dep) {
            npp_dev *= (1e0 + 0.11780208e+0 * tanh(t / 0.509312421e+02) + 0.24305130e-02 * tanh(t / 0.885326739e+01));
        }
        npp = npp_dev;
    }
    return npp;
}

/**
 * Compute the net primary production for the 4box model.
 *
 * @param ma        Mass of C in atmosphere.
 * @param co2_atm0  Preindustrial equilibrium CO2 concentration.
 * @param co2_dep   Co2 dependence.
 * @param deriv     Derivative dNPP/dm.
 * @return          Return the npp.
 */
double Land::compute4boxNpp(double &ma, double const &co2_atm0, bool &co2_dep, bool &deriv) const {
    double npp, m0, d_npp;
    if (!deriv) {
        if (co2_dep) {
            m0 = co2_atm0 * PPMTOGT;
            d_npp = this->_npp0 * this->_fert * UtilsMath::myLog(ma / m0);
            npp = d_npp + this->_npp0;
        } else {
            npp = this->_npp0;
        }
    } else {
        double npp_dev;
        if (co2_dep) {
            npp_dev = this->_npp0 * this->_fert * (1 / ma);
        } else {
            npp_dev = 0e0;
        }
        npp = npp_dev;
    }
    return npp;
}
