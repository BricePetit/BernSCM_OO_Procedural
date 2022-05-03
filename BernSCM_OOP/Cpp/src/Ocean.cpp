/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "Ocean.hpp"

#include <utility>

/**
 * Default constructor.
 */
Ocean::Ocean() = default;

/**
 * Constructor of the Ocean object.
 *
 * @param pirf          IRF coefficients refitted with 6 finite time scales.
 * @param hmix          Mixed layer depth (m).
 * @param cp            Heat capacity (J/kg/K).
 * @param dens          Water density (kg/m³).
 * @param dens_c        Water density (value used for DIC).
 * @param aoc           Sea Surface Area (m²).
 * @param kg_aoc        Gas exchange coefficient (1/yr).
 * @param t_chem        SST for calculating pCO2~DIC dependence (K); fixed, not=Temp.
 * @param model_name    The name of the model.
 * @param doc           The documentation about the model.
 */
Ocean::Ocean(Pirf pirf, double hmix, double cp, double dens, double dens_c, double aoc, double kg_aoc, double t_chem,
             string_view model_name, string doc) : _pirf(move(pirf)), _hmix(hmix), _cp(cp), _dens(dens),
                                                   _dens_c(dens_c), _aoc(aoc), _kg_aoc(kg_aoc), _t_chem(t_chem),
                                                   _model_name(model_name), _doc(move(doc)) {
    this->_om_t = PETA * SECTOYEAR / (this->_hmix * this->_cp * this->_dens * this->_aoc);
    this->_om_c = PETA / MUMOL / this->_dens_c / this->_hmix / this->_aoc;
}

/**
 * Move constructor.
 *
 * @param ocean Ocean to move.
 */
Ocean::Ocean(Ocean &&ocean) noexcept {
    this->_pirf = move(ocean._pirf);
    this->_hmix = ocean._hmix;
    this->_cp = ocean._cp;
    this->_dens = ocean._dens;
    this->_dens_c = ocean._dens_c;
    this->_aoc = ocean._aoc;
    this->_kg_aoc = ocean._kg_aoc;
    this->_t_chem = ocean._t_chem;
    this->_om_t = ocean._om_t;
    this->_om_c = ocean._om_c;
    this->_model_name = ocean._model_name;
    this->_doc = move(ocean._doc);
    this->_dpcs = ocean._dpcs;
    this->_prop = move(ocean._prop);
    this->_msk = ocean._msk;
    this->_ms = ocean._ms;
    this->_temp = ocean._temp;
    this->_temp_k = ocean._temp_k;
    this->_co2_exceed = ocean._co2_exceed;
    this->_t_exceed = ocean._t_exceed;

    // Delete previous value
    ocean._pirf = Pirf("", 0, {}, {});
    ocean._hmix = 0e0;
    ocean._cp = 0e0;
    ocean._dens = 0e0;
    ocean._dens_c = 0e0;
    ocean._aoc = 0e0;
    ocean._kg_aoc = 0e0;
    ocean._t_chem = 0e0;
    ocean._om_t = 0e0;
    ocean._om_c = 0e0;
    ocean._model_name = "";
    ocean._doc = "";
    ocean._dpcs = nullptr;
    ocean._prop = Prop(0, 0e0);
    ocean._msk = nullptr;
    ocean._ms = nullptr;
    ocean._temp = nullptr;
    ocean._temp_k = nullptr;
    ocean._co2_exceed = false;
    ocean._t_exceed = false;
}

/**
 * Move operator.
 *
 * @param ocean The object to move.
 * @return      Return the current object.
 */
Ocean &Ocean::operator=(Ocean &&ocean) noexcept {
    if (this != &ocean) {
        // Free existing resources
        delete[] this->_dpcs;
        delete[] this->_msk;
        delete[] this->_ms;
        delete[] this->_temp;
        delete[] this->_temp_k;
        // Copy data
        this->_pirf = move(ocean._pirf);
        this->_hmix = ocean._hmix;
        this->_cp = ocean._cp;
        this->_dens = ocean._dens;
        this->_dens_c = ocean._dens_c;
        this->_aoc = ocean._aoc;
        this->_kg_aoc = ocean._kg_aoc;
        this->_t_chem = ocean._t_chem;
        this->_om_t = ocean._om_t;
        this->_om_c = ocean._om_c;
        this->_model_name = ocean._model_name;
        this->_doc = move(ocean._doc);
        this->_dpcs = ocean._dpcs;
        this->_prop = move(ocean._prop);
        this->_msk = ocean._msk;
        this->_ms = ocean._ms;
        this->_temp = ocean._temp;
        this->_temp_k = ocean._temp_k;
        this->_co2_exceed = ocean._co2_exceed;
        this->_t_exceed = ocean._t_exceed;

        // Remove data from the other object
        ocean._pirf = Pirf("", 0, {}, {});
        ocean._hmix = 0e0;
        ocean._cp = 0e0;
        ocean._dens = 0e0;
        ocean._dens_c = 0e0;
        ocean._aoc = 0e0;
        ocean._kg_aoc = 0e0;
        ocean._t_chem = 0e0;
        ocean._om_t = 0e0;
        ocean._om_c = 0e0;
        ocean._model_name = "";
        ocean._doc = "";
        ocean._dpcs = nullptr;
        ocean._prop = Prop(0, 0e0);
        ocean._msk = nullptr;
        ocean._ms = nullptr;
        ocean._temp = nullptr;
        ocean._temp_k = nullptr;
        ocean._co2_exceed = false;
        ocean._t_exceed = false;
    }
    return *this;
}

/**
 * Destructor of the Ocean object.
 */
Ocean::~Ocean() {
    delete[] this->_dpcs;
    delete[] this->_msk;
    delete[] this->_ms;
    delete[] this->_temp_k;
    delete[] this->_temp;
}

/**
 * Initialize the stock of carbon.
 *
 * @param ntime The time dimension.
 */
void Ocean::initCarbonStock(int ntime) {
    this->_dpcs = new double[ntime]{0};

    // Initial ocean mixed layer C stock perturbation
    this->_msk = new double[this->_pirf.getNscale() + 1]{0};
    this->_ms = new double[ntime]{0};

    // Initial ocean temperature perturbation
    this->_temp_k = new double[this->_pirf.getNscale() + 1]{0};
    this->_temp = new double[ntime]{0};
}

/**
 * Increase the i-th value of temp_k by value.
 *
 * @param i     Index.
 * @param value New value to add at the temp_k array.
 */
void Ocean::increaseTempKElem(int i, double value) {
    this->_temp_k[i] += value;
}

/**
 * Set a new value at the index i of the dpcs array.
 *
 * @param i     The index.
 * @param value New value.
 */
void Ocean::setDpCsElem(int i, double value) {
    this->_dpcs[i] = value;
}

/**
 * Increase the i-th value of msk by value.
 *
 * @param i     Index.
 * @param value New value to add at the _msk array.
 */
void Ocean::increaseMskElem(int i, double value) {
    this->_msk[i] += value;
}

/**
 * Set a new value at the index i of the ms array.
 *
 * @param i     The index.
 * @param value New value.
 */
void Ocean::setMsElem(int i, double value) {
    this->_ms[i] = value;
}

/**
 * Reset all value of temp_k to 0.
 */
void Ocean::resetTempk() {
    delete[] this->_temp_k;
    this->_temp_k = new double[this->_pirf.getNscale() + 1]{0};

}

/**
 * Setter for the i-th element of temp.
 *
 * @param i     Index.
 * @param value Value.
 */
void Ocean::setTempElem(int i, double value) {
    this->_temp[i] = value;
}

/**
 * As the co2 exceed, we turn the variable to True.
 */
void Ocean::co2Exceed() {
    this->_co2_exceed = true;
}

/**
 * As the temperature exceeded, we turn the variable to True.
 */
void Ocean::tExceed() {
    this->_t_exceed = true;
}

/**
 * Getter for the aoc (Sea Surface Area (m²)) value.
 *
 * @return  Return the aoc.
 */
double Ocean::getAoc() const {
    return this->_aoc;
}

/**
 * Getter for the kg_aoc attribute (Gas exchange coefficient (1/yr)).
 *
 * @return  Return the kg_aoc.
 */
double Ocean::getKgAoc() const {
    return this->_kg_aoc;
}

/**
 * Getter for the _om_t value (Multiplier for heat uptake (K/(PW*yr))).
 *
 * @return  Return the om_t value.
 */
double Ocean::getOmt() const {
    return this->_om_t;
}

/**
 * Getter for the _om_c value (GtC→DIC conversion factor (umol/kg/Gt)).
 *
 * @return  Return the om_c value.
 */
double Ocean::getOmc() const {
    return this->_om_c;
}

/**
 * Getter for the documentation of the model.
 *
 * @return  Return the documentation.
 */
string Ocean::getDoc() const {
    return this->_doc;
}

/**
 * Getter for the i-th element of dpcs.
 *
 * @param i Index.
 * @return  Return the i-th value of dpcs.
 */
double Ocean::getDpCsElem(int i) const {
    return this->_dpcs[i];
}

/**
 * Getter for the propagator.
 *
 * @return  Return the propagator.
 */
Prop &Ocean::getPropagator() {
    return this->_prop;
}

/**
 * Getter for the nscale of the propagator.
 *
 * @return  Return the nscale of the propagator
 */
int Ocean::getPropNscale() const {
    return this->_prop.getNscale();
}

/**
 * Getter for the i-th element of propf attribute from the propagator.
 *
 * @param i Index.
 * @return  Return the i-th value of propf.
 */
double Ocean::getPropfElem(int i) const {
    return this->_prop.getPropfElem(i);
}

/**
 * Getter for the array msk.
 *
 * @return  Return the list msk.
 */
double *Ocean::getMsk() {
    return this->_msk;
}

/**
 * Getter for the i-th element of msk.
 *
 * @param i Index.
 * @return  Return the i-th element of msk.
 */
double Ocean::getMskElem(int i) const {
    return this->_msk[i];
}

/**
 * Getter for the i-th element of ms.
 *
 * @param i Index.
 * @return  Return the i-th value of ms.
 */
double Ocean::getMsElem(int i) const {
    return this->_ms[i];
}

/**
 * Getter for the array temp_k.
 *
 * @return  Return the list temp_k.
 */
double *Ocean::getTempK() {
    return this->_temp_k;
}

/**
 * Getter for the i-th element of temp_k.
 *
 * @param i Index.
 * @return  Return the i-th value of temp_k.
 */
double Ocean::getTempKElem(int i) const {
    return this->_temp_k[i];
}

/**
 * Getter for the i-th element of temp.
 *
 * @param i Index.
 * @return  Return the i-th value of temp.
 */
double Ocean::getTempElem(int i) const {
    return this->_temp[i];
}

/**
 * Check if the co2 exceed.
 *
 * @return  Check if it exceeded.
 */
bool Ocean::isCo2Exceed() const {
    return this->_co2_exceed;
}

/**
 * Check if the temperature exceed.
 *
 * @return  Check if it exceeded.
 */
bool Ocean::isTExceed() const {
    return this->_t_exceed;
}

/**
 * Calculate propagators.
 *
 * @param dt    Time step.
 */
void Ocean::computePropagators(double dt) {
    this->_prop.computePropagators(this->_pirf, dt);
}

/**
 * Analytical representation of the zeta-factor for a temperature
 * range from 17.7 to 18.3 degrees celsius (Tchem). Rogers chemistry model
 * was used to calculate the zeta-factor.
 *
 * @param t_dep     Temperature dependence.
 * @param co2_atm0  CO2 concentration for ocean exchange.
 * @param i         Index for the ms.
 * @param t         Global SAT change from preindustrial (℃).
 * @param deriv     Derivative dpCO2s/ddDIC.
 * @return          Ocean saturation CO2 pressure deviation from preindustrial
 *                  equilibrium (ppm), or derivative (d dpCs/d dDIC).
 */
double Ocean::computeDpCo2s(bool t_dep, double co2_atm0, int i, double t, bool deriv) {
    double d_dic = this->_ms[i] * this->_om_c, dp_co2s;
    if (!deriv) {
        // Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm)
        dp_co2s = (1.5568e0 - 1.3993e-2 * this->_t_chem) * d_dic +
                  (7.4706e0 - 0.20207e0 * this->_t_chem) * 1e-3 * pow(d_dic, 2) -
                  (1.2748e0 - 0.12015e0 * this->_t_chem) * 1e-5 * pow(d_dic, 3) +
                  (2.4491e0 - 0.12639e0 * this->_t_chem) * 1e-7 * pow(d_dic, 4) -
                  (1.5468e0 - 0.15326e0 * this->_t_chem) * 1e-10 * pow(d_dic, 5);
        if (dp_co2s > 1320 && !this->_co2_exceed) {
            cerr << "warning: CO2 parametrization range for dpCO2s exceeded (1320ppm)" << endl;
            this->_co2_exceed = true;
        }
        if (t + this->_t_chem > 25 && !this->_t_exceed) {
            cerr << "warning: temperature parametrization range for dpCO2s exceeded (25℃)" << endl;
            this->_t_exceed = true;
        }
        if (t_dep) {
            dp_co2s = (dp_co2s + co2_atm0) * exp(this->_buffer_t * t) - co2_atm0;
        }
    } else {
        double dp_co2s_dev;
        dp_co2s_dev =
                (1.5568e0 - 1.3993e-2 * this->_t_chem) + 2e0 * (7.4706 - 0.20207e0 * this->_t_chem) * 1.0e-3 * d_dic -
                3e0 * (1.2748 - 0.12015e0 * this->_t_chem) * 1.0e-5 * pow(d_dic, 2) +
                4e0 * (2.4491 - 0.12639e0 * this->_t_chem) * 1.0e-7 * pow(d_dic, 3) -
                5e0 * (1.5468 - 0.15326e0 * this->_t_chem) * 1.0e-10 * pow(d_dic, 4);
        if (t_dep) {
            dp_co2s_dev = dp_co2s_dev * exp(this->_buffer_t * t);
        }
        dp_co2s = dp_co2s_dev;
    }
    return dp_co2s;
}
