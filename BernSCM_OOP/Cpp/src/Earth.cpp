/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#include "Earth.hpp"

#include <utility>

Earth::Earth() = default;

/**
 * Constructor of the Earth object.
 *
 * @param ocean         The ocean model.
 * @param land          The land model.
 * @param atmosphere    Atmosphere of the earth.
 * @param forcing       All forcing values.
 * @param t2x           Climate sensitivity.
 * @param t_dep         Temperature dependence.
 * @param co2_dep       Co2 dependence.
 * @param scenario      The name of the scenario file without forcing_ (at the beginning) and .dat (at the end).
 * @param dt            Time step.
 */
Earth::Earth(Ocean &ocean, Land &land, Atmosphere &atmosphere, double **forcing, double t2x, bool t_dep, bool co2_dep,
             string scenario, double dt) : _ocean(move(ocean)), _land(move(land)), _atmosphere(move(atmosphere)),
                                           _forcing(forcing), _t2x(t2x), _t_dep(t_dep), _co2_dep(co2_dep),
                                           _scenario(move(scenario)), _dt(dt) {}

/**
 * Move operator.
 *
 * @param earth Object to move.
 * @return      Return the current object.
 */
Earth &Earth::operator=(Earth &&earth) noexcept {
    if (this != &earth) {
        // Free existing resources
        delete[] this->_time;
        delete[] this->_itime;
        delete[] this->_rf;
        delete[] this->_rfnc;
        delete[] this->_rfc;
        delete[] this->_rfb;
        delete[] this->_fh;
        delete[] this->_fb;
        delete[] this->_fo;
        delete[] this->_fnpp;
        // Copy data
        this->_ocean = move(earth._ocean);
        this->_land = move(earth._land);
        this->_atmosphere = move(earth._atmosphere);
        this->_forcing = earth._forcing;
        this->_t2x = earth._t2x;
        this->_t_dep = earth._t_dep;
        this->_co2_dep = earth._co2_dep;
        this->_scenario = move(earth._scenario);
        this->_dt = earth._dt;
        this->_ntime = earth._ntime;
        this->_time = earth._time;
        this->_itime = earth._itime;
        this->_rf = earth._rf;
        this->_rfnc = earth._rfnc;
        this->_rfc = earth._rfc;
        this->_rfb = earth._rfb;
        this->_fh = earth._fh;
        this->_fb = earth._fb;
        this->_fo = earth._fo;
        this->_fnpp = earth._fnpp;
        this->_f_budget = earth._f_budget;
        this->_rf_budget = earth._rf_budget;
        this->_co2_budget = earth._co2_budget;
        this->_exceed = earth._exceed;
        this->_t_exceed = earth._t_exceed;
        this->_co2_exceed = earth._co2_exceed;

        // Remove data from the other object
        earth._ocean = Ocean();
        earth._land = Land();
        earth._atmosphere = Atmosphere();
        earth._forcing = nullptr;
        earth._t2x = 0e0;
        earth._t_dep = false;
        earth._co2_dep = false;
        earth._scenario = "";
        earth._dt = 0e0;
        earth._ntime = 0;
        earth._time = nullptr;
        earth._itime = nullptr;
        earth._rf = nullptr;
        earth._rfnc = nullptr;
        earth._rfc = nullptr;
        earth._rfb = nullptr;
        earth._fh = nullptr;
        earth._fb = nullptr;
        earth._fo = nullptr;
        earth._fnpp = nullptr;
        earth._f_budget = false;
        earth._rf_budget = false;
        earth._co2_budget = false;
        earth._exceed = false;
        earth._t_exceed = false;
        earth._co2_exceed = false;
    }
    return *this;
}

/**
 * Destructor of the Earth object.
 */
Earth::~Earth() {
    delete[] this->_time;
    delete[] this->_itime;
    delete[] this->_fnpp;
    delete[] this->_rf;
    delete[] this->_rfnc;
    delete[] this->_rfc;
    delete[] this->_rfb;
    delete[] this->_fh;
    delete[] this->_fo;
}

/**
 * Initialization for the simulation.
 */
void Earth::initialize() {
    int nin = FORCING_ROWS - 1, j = 0;
    double tin[FORCING_ROWS];
    if (LINEAR && TEQUIL > 0) {
        throw "Equilibrated time steps not implemented for linear.";
    } else if (LINEAR && !(IMPLICIT_O || IMPLICIT_L) == 1) {
        throw "Linear scheme only implicit.";
    }
    for (int i = 0; i < FORCING_ROWS; i++) {
        tin[i] = this->_forcing[i][JTIME];
    }
    this->_ntime = int((tin[nin] - tin[0]) / this->_dt + 1e-6) + 1;
    this->_time = new double[this->_ntime]{0};
    this->_itime = new long int[this->_ntime]{0};

    for (int i = 0; i < this->_ntime; i++) {
        this->_time[i] = tin[0] + i * this->_dt;
        while (this->_time[i] > tin[j + 1]) {
            j += 1;
        }
        this->_itime[i] = j;
    }

    this->initAttributesToZero();
    this->_fnpp = new double[this->_ntime]{this->_land.getNpp0()};

    this->_land.computePropagators(this->_dt);
    this->_ocean.computePropagators(this->_dt);

    this->_land.initCarbonStock(this->_ntime);
    this->_ocean.initCarbonStock(this->_ntime);
}

/**
 * Init required attribute to zero.
 */
void Earth::initAttributesToZero() {
    this->_rf = new double[this->_ntime]{0};
    this->_rfnc = new double[this->_ntime]{0};
    this->_rfc = new double[this->_ntime]{0};
    this->_rfb = new double[this->_ntime]{0};
    this->_fh = new double[this->_ntime]{0};
    this->_fb = new double[this->_ntime]{0};
    this->_fo = new double[this->_ntime]{0};
    this->_atmosphere.initialize(this->_ntime);
}

/**
 * Set forcing value and verify budget cases.
 *
 * @param n Time counter.
 */
void Earth::setForcing(int n) {
    this->_ocean.setTempElem(n, this->interpol(n, JTEMP));
    this->_rfnc[n] = this->interpol(n, JRFNC);
    this->_rfb[n] = this->interpol(n, JRFB);
    this->_atmosphere.setMaElem(n, this->interpol(n, JACO2));
    this->_atmosphere.setECo2Elem(n, this->interpol(n, JECO2));
    this->_fb[n] = this->interpol(n, JFB);

    // Work out budget closure case
    if ((this->_atmosphere.getMaElem(n) != NA) || (this->_fb[n] == NA)) {
        this->_f_budget = true;
    } else {
        this->_f_budget = false;
    }
    if ((this->_ocean.getTempElem(n) != NA) || (this->_rfb[n] == NA)) {
        this->_rf_budget = true;
    } else {
        this->_rf_budget = false;
    }

    // Check budget cases (setting RF_CO2 is not implemented)
    if ((this->_atmosphere.getECo2Elem(n) == NA) || (this->_rfnc[n] == NA)) {
        throw "eCO2 and RF_nonC must always be set, use budget_RF and budget_sink to solve for RF/emissions " +
              to_string(this->_time[n]);
    }

    this->_co2_budget = false;
    if (this->_rf_budget && !this->_f_budget) {
        if (this->_ocean.getTempElem(n) == NA) {
            throw "glob_temp_dev not set when solving for budget_RF at year " + to_string(this->_time[n]);
        }
    } else if (!this->_rf_budget && this->_f_budget) {
        if (this->_atmosphere.getMaElem(n) == NA) {
            throw "CO2_atm not set when solving for budget_sink at year " + to_string(this->_time[n]);
        }
    } else if (this->_rf_budget && this->_f_budget) {
        if (this->_ocean.getTempElem(n) == NA) {
            throw "glob_temp_dev not set when solving for budget_RF at year " + to_string(this->_time[n]);
        }
        if (this->_atmosphere.getMaElem(n) == NA) {
            this->_co2_budget = true;
        }
    }
}

/**
 * Do the interpolation. It is the process of estimating unknown values that
 * fall between known values.
 *
 * @param n             Time counter.
 * @param index_forcing The index of the forcing value.
 * @return              Return the interpol value.
 */
double Earth::interpol(int n, int index_forcing) {
    double y;
    if (this->_forcing[this->_itime[n]][index_forcing] == NA) {
        if (abs(this->_time[n] - this->_forcing[this->_itime[n] + 1][JTIME]) < 1e-9) {
            y = this->_forcing[this->_itime[n] + 1][index_forcing];
        } else {
            y = NA;
        }
    } else if (this->_forcing[this->_itime[n] + 1][index_forcing] == NA) {
        if (abs(this->_time[n] - this->_forcing[this->_itime[n]][JTIME]) < 1e-9) {
            y = this->_forcing[this->_itime[n]][index_forcing];
        } else {
            y = NA;
        }
    } else {
        y = ((this->_time[n] - this->_forcing[this->_itime[n]][JTIME]) *
             this->_forcing[this->_itime[n] + 1][index_forcing] +
             (this->_forcing[this->_itime[n] + 1][JTIME] - this->_time[n]) *
             this->_forcing[this->_itime[n]][index_forcing]) /
            (this->_forcing[this->_itime[n] + 1][JTIME] - this->_forcing[this->_itime[n]][JTIME]);
    }
    return y;
}

/**
 * Set a new value at the i-th index of the budget radiative forcing (Wm⁻²).
 *
 * @param i     The index.
 * @param value The new value to set.
 */
void Earth::setRfbElem(int i, double value) {
    this->_rfb[i] = value;
}

/**
 * Setter for the budget C uptake (GtC/yr).
 *
 * @param i     The index.
 * @param value The new value to set.
 */
void Earth::setFbElem(int i, double value) {
    this->_fb[i] = value;
}

/**
 * Getter for the documentation of the ocean model.
 *
 * @return  Return the documentation.
 */
string Earth::getOceanDoc() const {
    return this->_ocean.getDoc();
}

/**
 * Getter for the i-th element of temp.
 *
 * @param i Index.
 * @return  Return the i-th value of temp.
 */
double Earth::getTempElem(int i) const {
    return this->_ocean.getTempElem(i);
}

/**
 * Getter for the i-th element of ms.
 *
 * @param i Index.
 * @return  Return the i-th value of ms.
 */
double Earth::getMsElem(int i) const {
    return this->_ocean.getMsElem(i);
}

/**
 * Getter for the documentation of the land model.
 *
 * @return  Return the documentation.
 */
string Earth::getLandDoc() const {
    return this->_land.getDoc();
}

/**
 * Getter for the npp0 attribute.
 *
 * @return  Return preindustrial net primary production (Gt/yr).
 */
double Earth::getNpp0() const {
    return this->_land.getNpp0();
}

/**
 * Getter for the i-th element of ml.
 *
 * @param i Index.
 * @return  Return the i-th value of ml.
 */
double Earth::getMlElem(int i) const {
    return this->_land.getMlElem(i);
}

/**
 * Getter for the i-th element of ma (mass of carbon in the atmosphere).
 *
 * @param i Index.
 * @return  Return the i-th value of ma.
 */
double Earth::getMaElem(int i) const {
    return this->_atmosphere.getMaElem(i);
}

/**
 * Getter for the preindustrial equilibrium CO2 concentration.
 *
 * @return  Return the co2_atm0 attribute.
 */
double Earth::getCo2Atm0() const {
    return this->_atmosphere.getCo2Atm0();
}

/**
 * Getter for the t2x (Climate sensitivity).
 *
 * @return  Return the Climate sensitivity.
 */
double Earth::getT2x() const {
    return this->_t2x;
}

/**
 * Check if it is temperature dependent.
 *
 * @return  Return true or false according to the temperature dependence.
 */
bool Earth::isTDep() const {
    return this->_t_dep;
}

/**
 * Check if it is co2 dependent.
 *
 * @return  Return true or false according to the co2 dependence.
 */
bool Earth::isCo2Dep() const {
    return this->_co2_dep;
}

/**
 * Getter for the scenario of the simulation.
 *
 * @return  Return the name of the scenario.
 */
string Earth::getScenario() const {
    return this->_scenario;
}

/**
 * Getter for the time step.
 *
 * @return  Return the time step.
 */
double Earth::getDt() const {
    return this->_dt;
}

/**
 * Getter for the time dimension.
 *
 * @return  Return the time dimension.
 */
int Earth::getNtime() const {
    return this->_ntime;
}

/**
 * Getter for the element i of the time.
 *
 * @param i Index.
 * @return  Return the i-th value of time.
 */
double Earth::getTimeElem(int i) const {
    return this->_time[i];
}

/**
 * Getter for the total radiative forcing (Wm⁻²) at the index i.
 *
 * @param i Index.
 * @return  Return the total radiative forcing (Wm⁻²).
 */
double Earth::getRfElem(int i) const {
    return this->_rf[i];
}

/**
 * Getter for the Non-CO₂ radiative forcing (Wm⁻²) at index i.
 *
 * @param i Index.
 * @return  Return the i-th Non-CO₂ radiative forcing (Wm⁻²).
 */
double Earth::getRfncElem(int i) const {
    return this->_rfnc[i];
}

/**
 * Getter for the CO₂ radiative forcing (Wm⁻²) at index i.
 *
 * @param i Index.
 * @return  Return the i-th CO₂ radiative forcing (Wm⁻²).
 */
double Earth::getRfcElem(int i) const {
    return this->_rfc[i];
}

/**
 * Getter for the Budget radiative forcing (Wm⁻²) at the index i.
 *
 * @param i Index.
 * @return  Return the i-th Budget radiative forcing (Wm⁻²).
 */
double Earth::getRfbElem(int i) const {
    return this->_rfb[i];
}

/**
 * Getter for the Air-sea heat flux (PW) at index i.
 *
 * @param i Index.
 * @return  Return the i-th Air-sea heat flux (PW).
 */
double Earth::getFhElem(int i) const {
    return this->_fh[i];
}

/**
 * Getter for the Budget C uptake (GtC/yr) at index i.
 *
 * @param i Index.
 * @return  Return the i-th Budget C uptake (GtC/yr).
 */
double Earth::getFbElem(int i) const {
    return this->_fb[i];
}

/**
 * Getter for the Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm) at index i.
 *
 * @param i Index.
 * @return  Return the i-th Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm).
 */
double Earth::getDpcsElem(int i) const {
    return this->_ocean.getDpCsElem(i);
}

/**
 * Getter for Air-sea C flux (GtC/yr) at index i.
 *
 * @param i Index.
 * @return  Return the i-th Air-sea C flux (GtC/yr).
 */
double Earth::getFoElem(int i) const {
    return this->_fo[i];
}

/**
 * Getter for list of NPP (GtC/yr).
 *
 * @return  Return the list of NPP.
 */
double *Earth::getFnpp() const {
    return this->_fnpp;
}

/**
 * Getter for the NPP (GtC/yr) at the index i.
 *
 * @param i Index.
 * @return  Return the i-th NPP (GtC/yr).
 */
double Earth::getFnppElem(int i) const {
    return this->_fnpp[i];
}

/**
 * Getter for CO₂ emissions (GtC/yr) at the index i.
 *
 * @param i Index.
 * @return  Return the i-th CO₂ emissions (GtC/yr).
 */
double Earth::getECo2Elem(int i) const {
    return this->_atmosphere.getECo2Elem(i);
}

/**
 * Main loop for the carbon cycle simulation.
 */
void Earth::carbonCycleClimateSimulation() {
    this->initialize();
    this->setForcing(0);
    this->_atmosphere.initializeCo2Atm0();
    for (int n = 1; n < this->_ntime; n++) {
        this->setForcing(n);
        this->timeStep(n);
    }
}

/**
 * Progress in the simulation step by step.
 *
 * @param n Correspond to the number of the current step.
 */
void Earth::timeStep(int n) {
    double nenner_u = 0e0, nenner_v = 0e0, nenner_w = 0e0, df_npp_dma = 0e0, ma_eq = 0e0, dm_ao = 0e0, e = 0e0, t_com;
    double tmp_sum, ml_com, ms_com, df_npp;
    int a; //counter
    if (IMPLICIT_O) {
        this->_fh[n] = this->_fh[n - 1];
        t_com = this->stepPulse(n, this->_fh, this->_ocean.getTempK(), this->_ocean.getPropagator(),
                                this->_ocean.getOmt());
    } else {
        t_com = this->stepPulse(n - 1, this->_fh, this->_ocean.getTempK(), this->_ocean.getPropagator(),
                                this->_ocean.getOmt());
    }

    // Ocean heat uptake
    if (this->_rf_budget) {
        // Update heat uptake (const. flux commitment)
        tmp_sum = 0;
        for (a = 0; a < this->_ocean.getPropNscale() + 1; a++) {
            tmp_sum += this->_ocean.getPropfElem(a);
        }
        this->_fh[n] =
                this->_fh[n - 1] + (this->_ocean.getTempElem(n) - t_com) / (this->_ocean.getOmt() * tmp_sum);
        // Update Tempk!
        for (a = 0; a < this->_ocean.getPropNscale() + 1; a++) {
            this->_ocean.increaseTempKElem(a, (this->_fh[n] - this->_fh[n - 1]) * this->_ocean.getPropfElem(a) *
                                              this->_ocean.getOmt());
        }
        // Current RF (W/m²)
        this->_rf[n] =
                this->_fh[n] / (this->_ocean.getAoc() / OFRAC / PETA) + this->_ocean.getTempElem(n) / this->_t2x * RF2X;
        if (this->_co2_budget) {
            // Calculate equivalent atmospheric CO2 (in GtC) from RF
            // It gives the atmospheric CO₂ (Gt).
            this->_atmosphere.setMaElem(n, exp((this->_rf[n] - this->_rfnc[n]) / RECO2) * CO2PREIND * PPMTOGT);
        }
    }
    if (LAND_T_DEP) {
        if (this->_t_dep) {
            if (this->_rf_budget) {
                this->_land.computePropTempDependent(
                        (this->_ocean.getTempElem(n - 1) + this->_ocean.getTempElem(n)) / 2e0, this->_dt);
            } else {
                this->_land.computePropTempDependent((this->_ocean.getTempElem(n - 1) + t_com) / 2e0, this->_dt);
            }
        }
        if (this->_land.isExceed() && !this->_exceed) {
            this->_exceed = true;
        }
    }
    // land C exchange
    if (this->_f_budget) {
        // solve for net C emissions
        if (LINEAR) {
            // (endyear value)
            if (this->_rf_budget) {
                // Use actual T
                this->_fnpp[n] = this->npp(this->_atmosphere.getMaElem(n), this->_ocean.getTempElem(n), false);
            } else {
                // Use committed T
                this->_fnpp[n] = this->npp(this->_atmosphere.getMaElem(n), t_com, false);
            }
        } else {
            // (Midyear value)
            if (this->_rf_budget) {
                this->_fnpp[n] = this->npp((this->_atmosphere.getMaElem(n) + this->_atmosphere.getMaElem(n - 1)) / 2e0,
                                           (this->_ocean.getTempElem(n) + this->_ocean.getTempElem(n - 1)) / 2e0,
                                           false);
            } else {
                this->_fnpp[n] = this->npp((this->_atmosphere.getMaElem(n) + this->_atmosphere.getMaElem(n - 1)) / 2e0,
                                           this->_ocean.getTempElem(n - 1), false);
            }
        }
    } else {
        // Anthro emissions
        e = (this->_atmosphere.getECo2Elem(n) + this->_atmosphere.getECo2Elem(n - 1)) / 2e0 -
            (this->_fb[n] + this->_fb[n - 1]) / 2e0;
        if (IMPLICIT_L) {
            // Auxiliary parameters
            // Commitment step with previous flux
            this->_fnpp[n] = this->_fnpp[n - 1];
            df_npp_dma = this->npp(this->_atmosphere.getMaElem(n - 1), this->_ocean.getTempElem(n - 1), true);
            tmp_sum = 0;
            for (a = 0; a < this->_land.getPropNscale() + 1; a++) {
                tmp_sum += this->_land.getPropfElem(a);
            }
            nenner_v = (df_npp_dma * tmp_sum + 1e0);
        }
    }
    if (IMPLICIT_L) {
        ml_com = this->stepPulse(n, this->_fnpp, this->_land.getMlk(), this->_land.getPropagator(), 1e0);
    } else {
        if (this->_f_budget) {
            ml_com = this->stepPulse(n, this->_fnpp, this->_land.getMlk(), this->_land.getPropagator(), 1e0);
        } else {
            ml_com = this->stepPulse(n - 1, this->_fnpp, this->_land.getMlk(), this->_land.getPropagator(), 1e0);
        }
    }
    if (IMPLICIT_O) {
        // Commitment step with current flux=0
        this->_fo[n] = 0;
        dm_ao = this->computeDpCo2s(n - 1, t_com, true) * PPMTOGT * this->_ocean.getOmc();
        ma_eq = (this->computeDpCo2s(n - 1, t_com, false) + this->_atmosphere.getCo2Atm0()) * PPMTOGT;
        tmp_sum = 0;
        for (a = 0; a < this->_ocean.getPropNscale() + 1; a++) {
            tmp_sum += this->_ocean.getPropfElem(a);
        }
        nenner_u = (this->_ocean.getKgAoc() * dm_ao * tmp_sum + 1e0);
        nenner_w = this->_dt * this->_ocean.getKgAoc();
    } else {
        this->_fo[n] = this->fasC(n - 1);
    }

    ms_com = this->stepPulse(n, this->_fo, this->_ocean.getMsk(), this->_ocean.getPropagator(), 1e0);

    if (IMPLICIT_L) {
        if (this->_f_budget) {
            this->_land.setMlElem(n, ml_com);
        } else {
            // Implicit step for flux change (zeroE commitment for ocean)
            tmp_sum = 0;
            for (a = 0; a < this->_ocean.getPropNscale() + 1; a++) {
                    tmp_sum += this->_ocean.getPropfElem(a);
            }
            df_npp = df_npp_dma / (nenner_u * nenner_v + nenner_w) *
                     (this->_land.getMlElem(n - 1) - ml_com + this->_dt * e + this->_dt * this->_ocean.getKgAoc() *
                                                                              (+ma_eq -
                                                                               this->_atmosphere.getMaElem(n - 1) +
                                                                               dm_ao *
                                                                               (ms_com - this->_ocean.getMsElem(n - 1) +
                                                                                tmp_sum *
                                                                                ((this->_land.getMlElem(n - 1) -
                                                                                  ml_com) / this->_dt + e))));
            tmp_sum = 0;
            for (a = 0; a < this->_land.getPropNscale() + 1; a++) {
                tmp_sum += this->_land.getPropfElem(a);
            }
            this->_land.setMlElem(n, df_npp * tmp_sum + ml_com);
        }
    } else {
        this->_land.setMlElem(n, ml_com);
    }

    if (IMPLICIT_O) {
        // Implicit ocean step
        if (this->_f_budget) {
            this->_fo[n] = this->_ocean.getKgAoc() *
                           (this->_atmosphere.getMaElem(n) - ma_eq - dm_ao * (ms_com - this->_ocean.getMsElem(n - 1))) /
                           nenner_u;
        } else {
            this->_fo[n] = this->_ocean.getKgAoc() / (nenner_u + nenner_w) *
                           (this->_atmosphere.getMaElem(n - 1) - ma_eq -
                            dm_ao * (ms_com - this->_ocean.getMsElem(n - 1)) -
                            (this->_land.getMlElem(n) - this->_land.getMlElem(n - 1)) + this->_dt * e);
        }
        tmp_sum = 0;
        for (a = 0; a < this->_ocean.getPropNscale() + 1; a++) {
            this->_ocean.increaseMskElem(a, this->_fo[n] * this->_ocean.getPropfElem(a));
            tmp_sum += this->_ocean.getMskElem(a);
        }
        this->_ocean.setMsElem(n, tmp_sum);
    } else {
        this->_ocean.setMsElem(n, ms_com);
    }

    // Total C budget
    if (!this->_f_budget) {
        if (LINEAR) {
            this->_atmosphere.setMaElem(n, this->_atmosphere.getMaElem(n - 1) +
                                           this->_dt * (e - (this->_fo[n - 1] + this->_fo[n]) / 2e0) -
                                           (this->_land.getMlElem(n) - this->_land.getMlElem(n - 1)));
        } else {
            this->_atmosphere.setMaElem(n, this->_atmosphere.getMaElem(n - 1) + this->_dt * (e - this->_fo[n]) -
                                           (this->_land.getMlElem(n) - this->_land.getMlElem(n - 1)));
        }
    }
    // update CO₂ RF
    this->_rfc[n] = this->_atmosphere.rfCo2(n);

    if (!this->_rf_budget) {
        // Update total RF (W/m²)
        this->_rf[n] = this->_rfc[n] + this->_rfnc[n] + this->_rfb[n];
        if (this->_t2x > 0e0) {
            if (IMPLICIT_O) {
                // Const flux commitment
                tmp_sum = 0;
                for (a = 0; a < this->_ocean.getPropNscale() + 1; a++) {
                    tmp_sum += this->_ocean.getPropfElem(a);
                }
                this->_fh[n] = (this->_rf[n] - RF2X * t_com / this->_t2x +
                                this->_fh[n - 1] * this->_ocean.getOmt() * tmp_sum * RF2X / this->_t2x) /
                               (RF2X / this->_t2x * this->_ocean.getOmt() * tmp_sum +
                                OFRAC * PETA / this->_ocean.getAoc());
                tmp_sum = 0;
                for (a = 0; a < this->_ocean.getPropNscale() + 1; a++) {
                    this->_ocean.increaseTempKElem(a, (this->_fh[n] - this->_fh[n - 1]) * this->_ocean.getPropfElem(a) *
                                                      this->_ocean.getOmt());
                    tmp_sum += this->_ocean.getTempKElem(a);
                }
                this->_ocean.setTempElem(n, tmp_sum);
            } else {
                this->_ocean.setTempElem(n, t_com);
                this->_fh[n] = this->fasT(n);
            }
        } else {
            this->_fh[n] = 0e0;
            this->_ocean.setTempElem(n, 0e0);
            this->_ocean.resetTempk();
        }
    }
    this->_fnpp[n] = this->npp(this->_atmosphere.getMaElem(n), this->_ocean.getTempElem(n), false);
    if (!this->_f_budget) {
        //Update mLk with updated NPP
        for (a = 0; a < this->_land.getPropNscale() + 1; a++) {
            this->_land.increaseMlkElem(a, (this->_fnpp[n] - this->_fnpp[n - 1]) * this->_land.getPropfElem(a));
        }
    }
    this->_ocean.setDpCsElem(n, this->computeDpCo2s(n, this->_ocean.getTempElem(n), false));
}

/**
 * Function to advance the integration of a tracer.
 *
 * @param n     Time index.
 * @param f     Flux to mixed layer.
 * @param mk    Boxes/tracer pools (input/output).
 * @param q     Propagators.
 * @param x     Variable-specific multiplier/conversion factor.
 * @return      Return the committed temperature change (K).
 */
double Earth::stepPulse(int n, const double *f, double *mk, Prop &q, double x) {
    double sum = 0;
    for (int j = 0; j < q.getNscale() + 1; j++) {
        mk[j] = mk[j] * q.getPropmElem(j) + f[n] * q.getPropfElem(j) * x;
        if (LINEAR) {
            mk[j] = mk[j] + f[n - 1] * q.getPropfoElem(j) * x;
        }
        sum += mk[j];
    }
    return sum;
}

/**
 * Compute the atmosphere-ocean CO2 flux (Gt/yr).
 *
 * @param i Index.
 * @return  Return the atmosphere-ocean CO2 flux (Gt/yr)
 */
double Earth::fasC(int i) {
    return this->_ocean.getKgAoc() * ((this->_atmosphere.getMaElem(i) - this->_atmosphere.getCo2Atm0() * PPMTOGT) -
                                      this->_ocean.getDpCsElem(i) * PPMTOGT);
}

/**
 * Compute the air-sea heat flux.
 *
 * @param i Index.
 * @return  Return the air-sea heat flux (PW).
 */
double Earth::fasT(int i) {
    double fas_t;
    if (this->_t2x > 0e0) {
        fas_t = (this->_ocean.getAoc() / OFRAC / PETA) *
                (this->_rf[i] - (this->_ocean.getTempElem(i) / this->_t2x) * RF2X);
    } else {
        fas_t = 0e0;
    }
    return fas_t;
}

/**
 * Compute the net primary production and update exceed values if needed.
 *
 * @param ma    Mass of C in atmosphere.
 * @param t     Global ΔSAT (℃).
 * @param deriv Derivative dNPP/dm.
 * @return      Return the net primary production.
 */
double Earth::npp(double ma, double t, bool deriv) {
    double npp;
    npp = this->_land.computeNpp(ma, t, this->_atmosphere.getCo2Atm0(), this->_co2_dep, this->_t_dep, deriv);
    this->updateExceedValues(this->_land, this->_ocean);
    return npp;
}

/**
 * Compute the DpCo2s and update exceed values.
 *
 * @param i     Index for the ms.
 * @param t     Global SAT change from preindustrial (℃).
 * @param deriv Derivative dpCO2s/ddDIC.
 * @return      Ocean saturation CO2 pressure deviation from preindustrial
 *              equilibrium (ppm), or derivative (d dpCs/d dDIC).
 */
double Earth::computeDpCo2s(int i, double t, bool deriv) {
    double dp_co2s;
    dp_co2s = this->_ocean.computeDpCo2s(this->_t_dep, this->_atmosphere.getCo2Atm0(), i, t, deriv);
    this->updateExceedValues(this->_ocean, this->_land);
    return dp_co2s;
}

/**
 * Update exceed values if needed.
 *
 * @tparam S                This type is Land or Ocean.
 * @tparam T                This type is Land or Ocean.
 * @param model             The model that have maybe modify the exceed values.
 * @param model_to_update   The model where we need to update the exceed values.
 */
template<typename S, typename T>
void Earth::updateExceedValues(S &model, T &model_to_update) {
    if (model.isCo2Exceed() and !this->_co2_exceed) {
        this->_co2_exceed = true;
        if (!model_to_update.isCo2Exceed()) {
            model_to_update.co2Exceed();
        }
    }
    if (model.isTExceed() and !this->_t_exceed) {
        this->_t_exceed = true;
        if (!model_to_update.isTExceed()) {
            model_to_update.tExceed();
        }
    }
}
