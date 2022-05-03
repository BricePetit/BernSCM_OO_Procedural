"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""
import sys

import numpy as np

import Globals
import Prop


class Ocean:
    """
    This class represents the ocean.
    """

    def __init__(self, pirf, hmix, cp, dens, dens_c, aoc, kg_aoc, t_chem, model_name, doc):
        """
        Constructor for the class ocean

        :param pirf:        IRF coefficients refitted with 6 finite time scales.
        :type pirf:         Pirf.
        :param hmix:        Mixed layer depth (m).
        :type hmix:         Float.
        :param cp:          Heat capacity (J/kg/K).
        :type cp:           Float.
        :param dens:        Water density (kg/m³).
        :type dens:         Float.
        :param dens_c:      Water density (value used for DIC).
        :type dens_c:       Float.
        :param aoc:         Sea Surface Area (m²).
        :type aoc:          Float.
        :param kg_aoc:      Gas exchange coefficient (1/yr).
        :type kg_aoc:       Float.
        :param t_chem:      SST for calculating pCO2~DIC dependence (K); fixed, not=Temp.
        :type t_chem:       Float.
        :param model_name:  The name of the model.
        :type model_name:   Str.
        :param doc:         The documentation about the model.
        :type doc:          Str.
        """
        self._pirf = pirf
        self._hmix = hmix
        self._cp = cp
        self._dens = dens
        self._dens_c = dens_c
        self._aoc = aoc
        self._kg_aoc = kg_aoc
        self._t_chem = t_chem
        # om_t: Multiplier for heat uptake (K/(PW*yr)).
        self._om_t = Globals.PETA * Globals.SECTOYEAR / (self._hmix * self._cp * self._dens * self._aoc)
        # om_c: GtC→DIC conversion factor (umol/kg/Gt).
        self._om_c = Globals.PETA / Globals.MUMOL / self._dens_c / self._hmix / self._aoc
        self._model_name = model_name
        self._doc = doc
        # ocean chemistry T-dependence (Takahashi 1993)
        self._buffer_t = 0.0423e0
        # Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm)
        self._dpcs = np.array([], dtype=np.longdouble)
        self._prop = Prop.Prop(0, 0e0)
        # Ocean C surface pools (GtC)
        self._msk = np.array([], dtype=np.longdouble)
        # Mass of C in ocean surface layer (GtC)
        self._ms = np.array([], dtype=np.longdouble)
        # Global mean SAT deviation from preindustrial (K)
        self._temp = np.array([], dtype=np.longdouble)
        # Global mean SAT box components
        self._temp_k = np.array([], dtype=np.longdouble)

        self._co2_exceed = False
        self._t_exceed = False

    def initCarbonStock(self, ntime):
        """
        Initialize the stock of carbon.

        :param ntime:   The time dimension.
        :type ntime:    Int.
        """
        self._dpcs = np.zeros(ntime, dtype=np.longdouble)

        # Initial ocean mixed layer C stock perturbation
        self._msk = np.zeros(self._pirf.getNscale() + 1, dtype=np.longdouble)
        self._ms = np.zeros(ntime, dtype=np.longdouble)

        # Initial ocean temperature perturbation
        self._temp_k = np.zeros(self._pirf.getNscale() + 1, dtype=np.longdouble)
        self._temp = np.zeros(ntime, dtype=np.longdouble)

    def increaseTempKElem(self, index, value):
        """
        Increase the i-th value of temp_k by value.

        :param index:   Index.
        :type index:    Int.
        :param value:   New value to add at the temp_k array.
        :type value:    Float.
        """
        self._temp_k[index] += value

    def setDpCsElem(self, i, value):
        """
        Set a new value at the index i of the dpcs array.

        :param i:       The index.
        :type i:        Int.
        :param value:   New value.
        :type value:    Float.
        """
        self._dpcs[i] = value

    def increaseMskElem(self, index, value):
        """
        Increase the i-th value of msk by value.

        :param index:   Index.
        :type bound:    Int.
        :param values:  New value.
        :type values:   Numpy array of float.
        """
        self._msk[index] += value

    def setMsElem(self, i, value):
        """
        Set a new value at the index i of the ms array.

        :param i:       The index.
        :type i:        Int.
        :param value:   New value.
        :type value:    Float.
        """
        self._ms[i] = value

    def resetTempk(self):
        """
        Reset all value of temp_k to 0.
        """
        self._temp_k = np.zeros(self._temp_k.size, dtype=np.longdouble)

    def setTempElem(self, i, value):
        """
        Setter for the i-th element of temp.

        :param i:       Index.
        :type i:        Int.
        :param value:   Value.
        :type value:    Float.
        """
        self._temp[i] = value

    def co2Exceed(self):
        """
        As the co2 exceed, we turn the variable to True.
        """
        self._co2_exceed = True

    def tExceed(self):
        """
        As the temperature exceeded, we turn the variable to True.
        """
        self._t_exceed = True

    def getAoc(self):
        """
        Getter for the aoc (Sea Surface Area (m²)) value.

        :return:    Return the aoc.
        :rtype:     Float.
        """
        return self._aoc

    def getKgAoc(self):
        """
        Getter for the kg_aoc attribute (Gas exchange coefficient (1/yr)).

        :return:    Return the kg_aoc.
        :rtype:     Float.
        """
        return self._kg_aoc

    def getOmt(self):
        """
        Getter for the _om_t value (Multiplier for heat uptake (K/(PW*yr))).

        :return:    Return the om_t value.
        :rtype:     Float.
        """
        return self._om_t

    def getOmc(self):
        """
        Getter for the _om_c value (GtC→DIC conversion factor (umol/kg/Gt)).

        :return:    Return the om_c value.
        :rtype:     Float.
        """
        return self._om_c

    def getDoc(self):
        """
        Getter for the documentation of the model.

        :return:    Return the documentation.
        :rtype:     Str.
        """
        return self._doc

    def getDpCsElem(self, i):
        """
        Getter for the i-th element of dpcs.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of dpcs.
        :rtype:     Float.
        """
        return self._dpcs[i]

    def getPropagator(self):
        """
        Getter for the propagator.

        :return:    Return the propagator.
        :rtype:     Prop.
        """
        return self._prop

    def getPropNscale(self):
        """
        Getter for the nscale of the propagator.

        :return:    Return the nscale of the propagator
        :rtype:     Int.
        """
        return self._prop.getNscale()

    def getPropfElem(self, index):
        """
        Getter for the i-th element of propf attribute from the propagator.s
        The propf array is an attribute of the propagator prop.

        :param index:   Index.
        :type index:    Int.

        :return:        Return the i-th value of propf.
        :rtype:         Float.
        """
        return self._prop.getPropfElem(index)

    def getMsk(self):
        """
        Getter for the array msk.

        :return:    Return the list msk.
        :rtype:     Numpy array of float.
        """
        return self._msk

    def getMskElem(self, index):
        """
        Getter for the i-th element of msk.

        :param index:   Index.
        :type index:    Int.

        :return:        Return the i-th element of msk.
        :rtype:         Float.
        """
        return self._msk[index]

    def getMsElem(self, i):
        """
        Getter for the i-th element of ms.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of ms.
        :rtype:     Float.
        """
        return self._ms[i]

    def getTempK(self):
        """
        Getter for the array temp_k.

        :return:    Return the list temp_k.
        :rtype:     Numpy array of float.
        """
        return self._temp_k

    def getTempKElem(self, index):
        """
        Getter for the first elements of the temp k array.

        :param index:   Index.
        :type index:    Int.

        :return:        Return the i-th value of temp_k.
        :rtype:         Float.
        """
        return self._temp_k[index]

    def getTempElem(self, i):
        """
        Getter for the i-th element of temp.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of temp.
        :rtype:     Float.
        """
        return self._temp[i]

    def isCo2Exceed(self):
        """
        Check if the co2 exceed.

        :return:    Check if it exceeded.
        :rtype:     Boolean.
        """
        return self._co2_exceed

    def isTExceed(self):
        """
        Check if the temperature exceed.

        :return:    Check if it exceeded.
        :rtype:     Boolean.
        """
        return self._t_exceed

    def computePropagators(self, dt):
        """
        Calculate propagators.

        :param dt:  Time step.
        :type dt:   Float.
        """
        self._prop.computePropagators(self._pirf, dt)

    def computeDpCo2s(self, t_dep, co2_atm0, i, t, deriv):
        """
        Analytical representation of the zeta-factor for a temperature 
        range from 17.7 to 18.3 degrees celsius (Tchem). Rogers chemistry model
        was used to calculate the zeta-factor.

        :param t_dep:       Temperature dependence.
        :type t_dep:        Boolean.
        :param co2_atm0:    CO2 concentration for ocean exchange.
        :type co2_atm0:     Float.
        :param i:           Index for the ms.
        :type i:            Int.
        :param t:           Global SAT change from preindustrial (℃).
        :type t:            Float.
        :param deriv:       Derivative dpCO2s/ddDIC.
        :type deriv:        Boolean.

        :return:            Ocean saturation CO2 pressure deviation from preindustrial 
                            equilibrium (ppm), or derivative (d dpCs/d dDIC).
        :rtype:             Float.
        """
        d_dic = self._ms[i] * self._om_c
        if not deriv:
            # Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm)
            dp_co2s = (1.5568e0 - 1.3993e-2 * self._t_chem) * d_dic \
                      + (7.4706e0 - 0.20207e0 * self._t_chem) * 1e-3 * d_dic ** 2 \
                      - (1.2748e0 - 0.12015e0 * self._t_chem) * 1e-5 * d_dic ** 3 \
                      + (2.4491e0 - 0.12639e0 * self._t_chem) * 1e-7 * d_dic ** 4 \
                      - (1.5468e0 - 0.15326e0 * self._t_chem) * 1e-10 * d_dic ** 5
            if dp_co2s > 1320 and not self._co2_exceed:
                print("warning: CO2 parametrization range for dpCO2s exceeded (1320ppm)", file=sys.stdout)
                self._co2_exceed = True
            if t + self._t_chem > 25 and not self._t_exceed:
                print("warning: temperature parametrization range for dpCO2s exceeded (25℃)", file=sys.stdout)
                self._t_exceed = True
            if t_dep:
                dp_co2s = (dp_co2s + co2_atm0) * np.exp(self._buffer_t * t, dtype=np.longdouble) - co2_atm0
        else:
            dp_co2s_dev = (1.5568e0 - 1.3993e-2 * self._t_chem) \
                          + 2e0 * (7.4706 - 0.20207e0 * self._t_chem) * 1.0e-3 * d_dic \
                          - 3e0 * (1.2748 - 0.12015e0 * self._t_chem) * 1.0e-5 * d_dic ** 2 \
                          + 4e0 * (2.4491 - 0.12639e0 * self._t_chem) * 1.0e-7 * d_dic ** 3 \
                          - 5e0 * (1.5468 - 0.15326e0 * self._t_chem) * 1.0e-10 * d_dic ** 4
            if t_dep:
                dp_co2s_dev = dp_co2s_dev * np.exp(self._buffer_t * t, dtype=np.longdouble)
            dp_co2s = dp_co2s_dev
        return dp_co2s
