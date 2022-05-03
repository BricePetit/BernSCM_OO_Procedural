"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""
import numpy as np

import Globals


class Earth:
    def __init__(self, ocean, land, atmosphere, forcing, t2x, t_dep, co2_dep, scenario, dt):
        """
        Constructor of the earth.

        :param ocean:       The ocean model.
        :type ocean:        Ocean.
        :param land:        The land model.
        :type land:         Land.
        :param atmosphere:  Atmosphere of the earth.
        :type atmosphere:   Atmosphere.
        :param forcing:     All forcing values.
        :type forcing:      Numpy array of float.
        :param t2x:         Climate sensitivity.
        :type t2x:          Float.
        :param t_dep:       Temperature dependence.
        :type t_dep:        Boolean.
        :param co2_dep:     Co2 dependence.
        :type co2_dep:      Boolean.
        :param scenario:    The name of the scenario file without forcing_ (at the beginning) and .dat (at the end).
        :type scenario:     Str.
        :param dt:          Time step.
        :type dt:           Float.
        """
        self._ocean = ocean
        self._land = land
        self._atmosphere = atmosphere
        # forcing is a matrix where the column represents:
        # [Year, glob_temp_dev, RF_nonCO2, RF_budget, co2_atm, fossil_CO2_em, budget_C_uptake]
        self._forcing = forcing
        self._t2x = t2x
        self._t_dep = t_dep
        self._co2_dep = co2_dep
        self._scenario = scenario
        self._dt = dt
        # Time dimension
        self._ntime = 0
        # Time model
        self._time = np.array([], dtype=np.longdouble)
        # Forcing time series
        self._itime = np.array([], dtype=np.longdouble)
        # Total radiative forcing (Wm⁻²)
        self._rf = np.array([], dtype=np.longdouble)
        # Non-CO₂ radiative forcing (Wm⁻²)
        self._rfnc = np.array([], dtype=np.longdouble)
        # CO₂ radiative forcing (Wm⁻²)
        self._rfc = np.array([], dtype=np.longdouble)
        # Budget radiative forcing (Wm⁻²)
        self._rfb = np.array([], dtype=np.longdouble)
        # Air-sea heat flux (PW)
        self._fh = np.array([], dtype=np.longdouble)
        # Budget C uptake (GtC/yr)
        self._fb = np.array([], dtype=np.longdouble)
        # Air-sea C flux (GtC/yr)
        self._fo = np.array([], dtype=np.longdouble)
        # NPP (GtC/yr)
        self._fnpp = np.array([], dtype=np.longdouble)
        # Budget C uptake boolean
        self._f_budget = False
        # Budget radiative forcing boolean
        self._rf_budget = False
        # CO2 budget boolean
        self._co2_budget = False

        self._exceed = False
        self._t_exceed = False
        self._co2_exceed = False

    def initialize(self):
        """
        Initialization for the simulation.
        """
        nin = self._forcing.shape[0] - 1
        if Globals.LINEAR and Globals.TEQUIL > 0:
            raise NameError("Equilibrated time steps not implemented for linear.")
        elif Globals.LINEAR and not ((Globals.IMPLICIT_O or Globals.IMPLICIT_L) == 1):
            raise NameError("Linear scheme only implicit.")
        tin = self._forcing[:, Globals.JTIME]
        self._ntime = int((tin[nin] - tin[0]) / self._dt + 1e-6) + 1
        j = 0
        self._time = np.zeros(self._ntime, dtype=np.longdouble)
        self._itime = np.zeros(self._ntime, dtype="int64")
        for i in range(self._ntime):
            self._time[i] = tin[0] + i * self._dt
            while self._time[i] > tin[j + 1]:
                j += 1
            self._itime[i] = j

        self.initAttributesToZero()
        self._fnpp = np.full(self._ntime, self._land.getNpp0(), dtype=np.longdouble)

        self._land.computePropagators(self._dt)
        self._ocean.computePropagators(self._dt)

        self._land.initCarbonStock(self._ntime)
        self._ocean.initCarbonStock(self._ntime)

    def initAttributesToZero(self):
        """
        Init required attribute to zero.
        """
        self._rf = np.zeros(self._ntime, dtype=np.longdouble)
        self._rfnc = np.zeros(self._ntime, dtype=np.longdouble)
        self._rfc = np.zeros(self._ntime, dtype=np.longdouble)
        self._rfb = np.zeros(self._ntime, dtype=np.longdouble)
        self._fh = np.zeros(self._ntime, dtype=np.longdouble)
        self._fb = np.zeros(self._ntime, dtype=np.longdouble)
        self._fo = np.zeros(self._ntime, dtype=np.longdouble)
        self._atmosphere.initialize(self._ntime)

    def setForcing(self, n):
        """
        Set forcing value and verify budget cases.

        :param n:   Time counter.
        :type n:    Int.
        """
        self._ocean.setTempElem(n, self.interpol(n, Globals.JTEMP))
        self._rfnc[n] = self.interpol(n, Globals.JRFNC)
        self._rfb[n] = self.interpol(n, Globals.JRFB)
        self._atmosphere.setMaElem(n, self.interpol(n, Globals.JACO2))
        self._atmosphere.setECo2Elem(n, self.interpol(n, Globals.JECO2))
        self._fb[n] = self.interpol(n, Globals.JFB)

        # Work out budget closure case
        if (self._atmosphere.getMaElem(n) != Globals.NA) or (self._fb[n] == Globals.NA):
            self._f_budget = True
        else:
            self._f_budget = False

        if (self._ocean.getTempElem(n) != Globals.NA) or (self._rfb[n] == Globals.NA):
            self._rf_budget = True
        else:
            self._rf_budget = False

        # Check budget cases (setting RF_CO2 is not implemented)
        if (self._atmosphere.getECo2Elem(n) == Globals.NA) or (self._rfnc[n] == Globals.NA):
            raise Exception(
                "eCO2 and RF_nonC must always be set, use budget_RF and budget_sink to solve for RF/emissions {}".format(
                    self._time[n]))

        self._co2_budget = False
        if self._rf_budget and not self._f_budget:
            if self._ocean.getTempElem(n) == Globals.NA:
                raise Exception("glob_temp_dev not set when solving for budget_RF at year {}".format(self._time[n]))
        elif not self._rf_budget and self._f_budget:
            if self._atmosphere.getMaElem(n) == Globals.NA:
                raise Exception("CO2_atm not set when solving for budget_sink at year {}".format(self._time[n]))
        elif self._rf_budget and self._f_budget:
            if self._ocean.getTempElem(n) == Globals.NA:
                raise Exception("glob_temp_dev not set when solving for budget_RF at year {}".format(self._time[n]))
            if self._atmosphere.getMaElem(n) == Globals.NA:
                self._co2_budget = True

    def interpol(self, n, index_forcing):
        """
        Do the interpolation. It is the process of estimating unknown values that 
        fall between known values.

        :param n:               Time counter.
        :type n:                Int.
        :param index_forcing:   The index of the forcing value.
        :type index_forcing:    Int.

        :return:                Return the interpol value.
        :rtype:                 Float.
        """
        if self._forcing[self._itime[n]][index_forcing] == Globals.NA:
            if abs(self._time[n] - self._forcing[self._itime[n] + 1][Globals.JTIME]) < 1e-9:
                y = self._forcing[self._itime[n] + 1][index_forcing]
            else:
                y = Globals.NA
        elif self._forcing[self._itime[n] + 1][index_forcing] == Globals.NA:
            if abs(self._time[n] - self._forcing[self._itime[n]][Globals.JTIME]) < 1e-9:
                y = self._forcing[self._itime[n]][index_forcing]
            else:
                y = Globals.NA
        else:
            y = ((self._time[n] - self._forcing[self._itime[n]][Globals.JTIME]) *
                 self._forcing[self._itime[n] + 1][index_forcing] +
                 (self._forcing[self._itime[n] + 1][Globals.JTIME] - self._time[n]) *
                 self._forcing[self._itime[n]][index_forcing]) / \
                (self._forcing[self._itime[n] + 1][Globals.JTIME] - self._forcing[self._itime[n]][Globals.JTIME])
        return y

    def setRfbElem(self, i, value):
        """
        Set a new value at the i-th index of the budget radiative forcing (Wm⁻²).

        :param i:       The index.
        :type i:        Int.
        :param value:   The new value to set.
        :type value:    Float.
        """
        self._rfb[i] = value

    def setFbElem(self, i, value):
        """
        Setter for the budget C uptake (GtC/yr).

        :param i:       The index.
        :type i:        Int.
        :param value:   The new value to set.
        :type value:    Float.
        """
        self._fb[i] = value

    def getOceanDoc(self):
        """
        Getter for the documentation of the ocean model.

        :return:    Return the documentation.
        :rtype:     Str.
        """
        return self._ocean.getDoc()

    def getTempElem(self, i):
        """
        Getter for the i-th element of temp.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of temp.
        :rtype:     Float.
        """
        return self._ocean.getTempElem(i)

    def getMsElem(self, i):
        """
        Getter for the i-th element of ms.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of ms.
        :rtype:     Float.
        """
        return self._ocean.getMsElem(i)

    def getLandDoc(self):
        """
        Getter for the documentation of the land model.

        :return:    Return the documentation.
        :rtype:     Str.
        """
        return self._land.getDoc()

    def getNpp0(self):
        """
        Getter for the npp0 attribute.

        :return:    Return preindustrial net primary production (Gt/yr).
        :rtype:     Float.
        """
        return self._land.getNpp0()

    def getMlElem(self, i):
        """
        Getter for the i-th element of ml.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of ml.
        :rtype:     Float.
        """
        return self._land.getMlElem(i)

    def getMaElem(self, i):
        """
        Getter for the i-th element of ma (mass of carbon in the atmosphere).

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of ma.
        :rtype:     Float.
        """
        return self._atmosphere.getMaElem(i)

    def getCo2Atm0(self):
        """
        Getter for the preindustrial equilibrium CO2 concentration.

        :return:    Return the co2_atm0 attribute.
        :rtype:     Float.
        """
        return self._atmosphere.getCo2Atm0()

    def getT2x(self):
        """
        Getter for the t2x (Climate sensitivity).

        :return:    Return the Climate sensitivity.
        :rtype:     Float.
        """
        return self._t2x

    def isTDep(self):
        """
        Check if it is temperature dependent.

        :return:    Return true or false according to the temperature dependence.
        :rtype:     Boolean.
        """
        return self._t_dep

    def isCo2Dep(self):
        """
        Check if it is co2 dependent.

        :return:    Return true or false according to the co2 dependence.
        :rtype:     Boolean.
        """
        return self._co2_dep

    def getScenario(self):
        """
        Getter for the scenario of the simulation.

        :return:    Return the name of the scenario.
        :rtype:     Str.
        """
        return self._scenario

    def getDt(self):
        """
        Getter for the time step.

        :return:    Return the time step.
        :rtype:     Float.
        """
        return self._dt

    def getNtime(self):
        """
        Getter for the time dimension.

        :return:    Return the time dimension.
        :rtype:     Int.
        """
        return self._ntime

    def getTimeElem(self, i):
        """
        Getter for the element i of the time.
 
        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of time.
        :rtype:     Float.
        """
        return self._time[i]

    def getRfElem(self, i):
        """
        Getter for the total radiative forcing (Wm⁻²) at the index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the total radiative forcing (Wm⁻²).
        :rtype:     Float.
        """
        return self._rf[i]

    def getRfncElem(self, i):
        """
        Getter for the Non-CO₂ radiative forcing (Wm⁻²) at index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th Non-CO₂ radiative forcing (Wm⁻²).
        :rtype:     Float.
        """
        return self._rfnc[i]

    def getRfcElem(self, i):
        """
        Getter for the CO₂ radiative forcing (Wm⁻²) at index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th CO₂ radiative forcing (Wm⁻²).
        :rtype:     Float.
        """
        return self._rfc[i]

    def getRfbElem(self, i):
        """
        Getter for the Budget radiative forcing (Wm⁻²) at the index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th Budget radiative forcing (Wm⁻²).
        :rtype:     Float.
        """
        return self._rfb[i]

    def getFhElem(self, i):
        """
        Getter for the Air-sea heat flux (PW) at index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th Air-sea heat flux (PW).
        :rtype:     Float.
        """
        return self._fh[i]

    def getFbElem(self, i):
        """
        Getter for the Budget C uptake (GtC/yr) at index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th Budget C uptake (GtC/yr).
        :rtype:     Float.
        """
        return self._fb[i]

    def getDpcsElem(self, i):
        """
        Getter for the Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm) at index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm).
        :rtype:     Float.
        """
        return self._ocean.getDpCsElem(i)

    def getFoElem(self, i):
        """
        Getter for Air-sea C flux (GtC/yr) at index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th Air-sea C flux (GtC/yr).
        :rtype:     Float.
        """
        return self._fo[i]

    def getFnpp(self):
        """
        Getter for list of NPP (GtC/yr).

        :return:    Return the list of NPP.
        :rtype:     Numpy array of float.
        """
        return self._fnpp

    def getFnppElem(self, i):
        """
        Getter for the NPP (GtC/yr) at the index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th NPP (GtC/yr).
        :rtype:     Float.
        """
        return self._fnpp[i]

    def getECo2Elem(self, i):
        """
        Getter for CO₂ emissions (GtC/yr) at the index i.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th CO₂ emissions (GtC/yr).
        :rtype:     Float.
        """
        return self._atmosphere.getECo2Elem(i)

    def carbonCycleClimateSimulation(self):
        """
        Main loop for the carbon cycle simulation.
        """
        self.initialize()
        self.setForcing(0)
        self._atmosphere.initializeCo2Atm0()
        for n in range(1, self._ntime):
            self.setForcing(n)
            self.timeStep(n)
        

    def timeStep(self, n):
        """
        Progress in the simulation step by step.

        :param n:   Correspond to the number of the current step.
        :type n:    Int.
        """
        nenner_u = 0e0
        nenner_v = 0e0
        nenner_w = 0e0
        df_npp_dma = 0e0
        ma_eq = 0e0
        dm_ao = 0e0
        e = 0e0
        if Globals.IMPLICIT_O:
            self._fh[n] = self._fh[n - 1]
            t_com = self.stepPulse(n, self._fh, self._ocean.getTempK(), self._ocean.getPropagator(),
                                              self._ocean.getOmt())
        else:
            t_com = self.stepPulse(n - 1, self._fh, self._ocean.getTempK(), self._ocean.getPropagator(),
                                              self._ocean.getOmt())

        # Ocean heat uptake
        if self._rf_budget:
            # Update heat uptake (const. flux commitment)
            tmp_sum = 0
            for a in range(self._ocean.getPropNscale() + 1):
                tmp_sum += self._ocean.getPropfElem(a)
            self._fh[n] = self._fh[n - 1] + (self._ocean.getTempElem(n) - t_com) / (self._ocean.getOmt() * tmp_sum)
            # Update Tempk!
            for a in range(self._ocean.getPropNscale() + 1):
                self._ocean.increaseTempKElem(a, (self._fh[n] - self._fh[n - 1]) * self._ocean.getPropfElem(a) * self._ocean.getOmt())
            # Current RF (W/m²)
            self._rf[n] = self._fh[n] / (self._ocean.getAoc() / Globals.OFRAC
                                         / Globals.PETA) + self._ocean.getTempElem(n) / self._t2x * Globals.RF2X

            if self._co2_budget:
                # Calculate equivalent atmospheric CO2 (in GtC) from RF
                # It gives the atmospheric CO₂ (Gt).
                self._atmosphere.setMaElem(n, np.exp((self._rf[n] - self._rfnc[n]) / Globals.RECO2, dtype=np.longdouble)
                                           * Globals.CO2PREIND * Globals.PPMTOGT)

        if Globals.LAND_T_DEP:
            if self._t_dep:
                if self._rf_budget:
                    self._land.computePropTempDependent(
                        (self._ocean.getTempElem(n - 1) + self._ocean.getTempElem(n)) / 2e0, self._dt)
                else:
                    self._land.computePropTempDependent((self._ocean.getTempElem(n - 1) + t_com) / 2e0, self._dt)
            if self._land.isExceed() and not self._exceed:
                self._exceed = True

        # land C exchange
        if self._f_budget:
            # solve for net C emissions
            if Globals.LINEAR:
                # (endyear value)
                if self._rf_budget:
                    # Use actual T
                    self._fnpp[n] = self.npp(self._atmosphere.getMaElem(n), self._ocean.getTempElem(n), False)
                else:
                    # Use committed T
                    self._fnpp[n] = self.npp(self._atmosphere.getMaElem(n), t_com, False)
            else:
                # (Midyear value)
                if self._rf_budget:
                    self._fnpp[n] = self.npp((self._atmosphere.getMaElem(n) + self._atmosphere.getMaElem(n - 1)) / 2e0,
                                             (self._ocean.getTempElem(n) + self._ocean.getTempElem(n - 1)) / 2e0, False)
                else:
                    self._fnpp[n] = self.npp((self._atmosphere.getMaElem(n) + self._atmosphere.getMaElem(n - 1)) / 2e0,
                                             self._ocean.getTempElem(n - 1), False)
        else:
            # Anthro emissions
            e = (self._atmosphere.getECo2Elem(n) + self._atmosphere.getECo2Elem(n - 1)) / 2e0 \
                - (self._fb[n] + self._fb[n - 1]) / 2e0
            if Globals.IMPLICIT_L:
                # Auxiliary parameters
                # Commitment step with previous flux
                self._fnpp[n] = self._fnpp[n - 1]
                df_npp_dma = self.npp(self._atmosphere.getMaElem(n - 1), self._ocean.getTempElem(n - 1), True)
                tmp_sum = 0
                for a in range(self._land.getPropNscale() + 1):
                    tmp_sum += self._land.getPropfElem(a)
                nenner_v = (df_npp_dma * tmp_sum + 1e0)

        if Globals.IMPLICIT_L:
            ml_com = self.stepPulse(n, self._fnpp, self._land.getMlk(), self._land.getPropagator(), 1e0)
        else:
            if self._f_budget:
                ml_com = self.stepPulse(n, self._fnpp, self._land.getMlk(), self._land.getPropagator(), 1e0)
            else:
                ml_com = self.stepPulse(n - 1, self._fnpp, self._land.getMlk(), self._land.getPropagator(), 1e0)

        if Globals.IMPLICIT_O:
            # Commitment step with current flux=0
            self._fo[n] = 0
            dm_ao = self.computeDpCo2s(n - 1, t_com, True) * Globals.PPMTOGT * self._ocean.getOmc()
            ma_eq = (self.computeDpCo2s(n - 1, t_com, False) + self._atmosphere.getCo2Atm0()) * Globals.PPMTOGT
            tmp_sum = 0
            for a in range(self._ocean.getPropNscale() + 1):
                tmp_sum += self._ocean.getPropfElem(a)
            nenner_u = (self._ocean.getKgAoc() * dm_ao * tmp_sum + 1e0)
            nenner_w = self._dt * self._ocean.getKgAoc()
        else:
            self._fo[n] = self.fasC(n - 1)

        ms_com = self.stepPulse(n, self._fo, self._ocean.getMsk(), self._ocean.getPropagator(), 1e0)


        if Globals.IMPLICIT_L:
            if self._f_budget:
                self._land.setMlElem(n, ml_com)
            else:
                # Implicit step for flux change (zeroE commitment for ocean)
                tmp_sum = 0
                for a in range(self._ocean.getPropNscale() + 1):
                    tmp_sum += self._ocean.getPropfElem(a)
                df_npp = df_npp_dma / (nenner_u * nenner_v + nenner_w) * (
                        self._land.getMlElem(n - 1) - ml_com + self._dt * e
                        + self._dt * self._ocean.getKgAoc()
                        * (+ma_eq - self._atmosphere.getMaElem(n - 1) + dm_ao
                           * (ms_com - self._ocean.getMsElem(n - 1)
                              + tmp_sum * ((self._land.getMlElem(n - 1) - ml_com) / self._dt + e))))
                tmp_sum = 0
                for a in range(self._land.getPropNscale() + 1):
                    tmp_sum += self._land.getPropfElem(a)
                self._land.setMlElem(n, df_npp * tmp_sum + ml_com)
        else:
            self._land.setMlElem(n, ml_com)

        if Globals.IMPLICIT_O:
            # Implicit ocean step
            if self._f_budget:
                self._fo[n] = self._ocean.getKgAoc() * (self._atmosphere.getMaElem(n) - ma_eq - dm_ao
                                                        * (ms_com - self._ocean.getMsElem(n - 1))) / nenner_u
            else:
                self._fo[n] = self._ocean.getKgAoc() / (nenner_u + nenner_w) \
                              * (self._atmosphere.getMaElem(n - 1)
                                 - ma_eq - dm_ao * (ms_com - self._ocean.getMsElem(n - 1)) - (
                                         self._land.getMlElem(n) - self._land.getMlElem(n - 1)) + self._dt * e)
            tmp_sum = 0 
            for a in range(self._ocean.getPropNscale() + 1):
                self._ocean.increaseMskElem(a, self._fo[n] * self._ocean.getPropfElem(a))
                tmp_sum += self._ocean.getMskElem(a)
            self._ocean.setMsElem(n, tmp_sum)
        else:
            self._ocean.setMsElem(n, ms_com)

        # Total C budget
        if not self._f_budget:
            if Globals.LINEAR:
                self._atmosphere.setMaElem(n, self._atmosphere.getMaElem(n - 1) + self._dt
                                           * (e - (self._fo[n - 1] + self._fo[n]) / 2e0)
                                           - (self._land.getMlElem(n) - self._land.getMlElem(n - 1)))
            else:
                self._atmosphere.setMaElem(n, self._atmosphere.getMaElem(n - 1) + self._dt * (e - self._fo[n])
                                           - (self._land.getMlElem(n) - self._land.getMlElem(n - 1)))
        # update CO₂ RF
        self._rfc[n] = self._atmosphere.rfCo2(n)

        if not self._rf_budget:
            # Update total RF (W/m²)
            self._rf[n] = self._rfc[n] + self._rfnc[n] + self._rfb[n]
            if self._t2x > 0e0:
                if Globals.IMPLICIT_O:
                    # Const flux commitment
                    tmp_sum = 0
                    for a in range(self._ocean.getPropNscale() + 1):
                        tmp_sum += self._ocean.getPropfElem(a)
                    self._fh[n] = (self._rf[n] - Globals.RF2X * t_com / self._t2x + self._fh[n - 1]
                                   * self._ocean.getOmt()
                                   * tmp_sum
                                   * Globals.RF2X / self._t2x) \
                                  / (Globals.RF2X / self._t2x * self._ocean.getOmt()
                                     * tmp_sum
                                     + Globals.OFRAC * Globals.PETA / self._ocean.getAoc())
                    tmp_sum = 0
                    for a in range(self._ocean.getPropNscale() + 1):
                        self._ocean.increaseTempKElem(a, (self._fh[n] - self._fh[n - 1]) * self._ocean.getPropfElem(a) 
                                              * self._ocean.getOmt())
                        tmp_sum += self._ocean.getTempKElem(a)
                    self._ocean.setTempElem(n, tmp_sum)
                else:
                    self._ocean.setTempElem(n, t_com)
                    self._fh[n] = self.fasT(n)
            else:
                self._fh[n] = 0e0
                self._ocean.setTempElem(n, 0e0)
                self._ocean.resetTempk()
        self._fnpp[n] = self.npp(self._atmosphere.getMaElem(n), self._ocean.getTempElem(n), False)
        if not self._f_budget:
            # Update mLk with updated NPP
            for a in range(self._land.getPropNscale() + 1):
                self._land.increaseMlkElem(a, (self._fnpp[n] - self._fnpp[n - 1]) * self._land.getPropfElem(a))
        self._ocean.setDpCsElem(n, self.computeDpCo2s(n, self._ocean.getTempElem(n), False))

    def stepPulse(self, n, f, mk, q, x):
        """
        Function to advance the integration of a tracer.

        :param n:   Time index.
        :type n:    Int.
        :param f:   Flux to mixed layer.
        :type f:    Numpy array of float.
        :param mk:  Boxes/tracer pools (input/output).
        :type mk:   Numpy array of float.
        :param q:   Propagators.
        :type q:    Prop.
        :param x:   Variable-specific multiplier/conversion factor.
        :type x:    Float.

        :return:    Return the committed temperature change (K).
        :rtype:     Float.
        """
        tot_sum = 0
        for j in range(q.getNscale() + 1):
            mk[j] = mk[j] * q.getPropmElem(j) + f[n] * q.getPropfElem(j) * x
            if Globals.LINEAR:
                mk[j] = mk[j] + f[n - 1] * q.getPropfoElem(j) * x
            tot_sum += mk[j]
        return tot_sum

    def fasC(self, i):
        """
        Compute the atmosphere-ocean CO2 flux (Gt/yr).

        :param i:   Index.
        :type i:    Int.

        :return:    Return the atmosphere-ocean CO2 flux (Gt/yr)
        :rtype:     Float.
        """
        return self._ocean.getKgAoc() * ((self._atmosphere.getMaElem(i)
                                          - self._atmosphere.getCo2Atm0() * Globals.PPMTOGT)
                                         - self._ocean.getDpCsElem(i) * Globals.PPMTOGT)

    def fasT(self, i):
        """
        Compute the air-sea heat flux.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the air-sea heat flux (PW).
        :rtype:     Float.
        """
        if self._t2x > 0e0:
            fas_t = (self._ocean.getAoc() / Globals.OFRAC / Globals.PETA) \
                    * (self._rf[i] - (self._ocean.getTempElem(i) / self._t2x) * Globals.RF2X)
        else:
            fas_t = 0e0
        return fas_t

    def npp(self, ma, t, deriv):
        """
        Compute the net primary production and update exceed values if needed.

        :param ma:      Mass of C in atmosphere.
        :type ma:       Float.
        :param t:       Global ΔSAT (℃).
        :type t:        Float.
        :param deriv:   Derivative dNPP/dm.
        :type deriv:    Boolean.

        :return:        Return the net primary production.
        :rtype:         Float.
        """
        npp = self._land.computeNpp(ma, t, self._atmosphere.getCo2Atm0(), self._co2_dep, self._t_dep, deriv)
        self.updateExceedValues(self._land, self._ocean)
        return npp

    def computeDpCo2s(self, i, t, deriv):
        """
        Compute the DpCo2s and update exceed values.

        :param i:       Index for the ms.
        :type i:        Int.
        :param t:       Global SAT change from preindustrial (℃).
        :type t:        Float.
        :param deriv:   Derivative dpCO2s/ddDIC.
        :type deriv:    Boolean.

        :return:        Ocean saturation CO2 pressure deviation from preindustrial
                        equilibrium (ppm), or derivative (d dpCs/d dDIC).
        :rtype:         Float.
        """
        dp_co2s = self._ocean.computeDpCo2s(self._t_dep, self._atmosphere.getCo2Atm0(), i, t, deriv)
        self.updateExceedValues(self._ocean, self._land)
        return dp_co2s

    def updateExceedValues(self, model, model_to_update):
        """
        Update exceed values if needed.

        :param model:           The model that have maybe modify the exceed values.
        :type model:            Land or Ocean.
        :param model_to_update: The model where we need to update the exceed values.
        :type model_to_update:  Land or Ocean.
        """
        if model.isCo2Exceed() and not self._co2_exceed:
            self._co2_exceed = True
            if not model_to_update.isCo2Exceed():
                model_to_update.co2Exceed()
        if model.isTExceed() and not self._t_exceed:
            self._t_exceed = True
            if not model_to_update.isTExceeded():
                model_to_update.tExceed()
