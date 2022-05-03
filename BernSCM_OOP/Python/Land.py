"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""

import copy as cp
import sys

import numpy as np

import Globals
import Prop
from UtilsMath import UtilsMath


class Land:
    """
    This class represents the land.
    """

    def __init__(self, pirf, sirf, npp0, fert, model_name, doc):
        """
        Constructor of the land class.

        :param pirf:        IRF coefficients.
        :type pirf:         Pirf.
        :param sirf:        IRF sensitivites. Can be none according to the model.
        :type sirf:         Sirf.
        :param npp0:        Preindustrial (potential natural) NPP (net primary production) (Gt/yr).
        :type npp0:         Float.
        :param fert:        Fertilization Coefficient, only in the npp_4box. Can be none according to the model.
        :type fert:         Float.
        :param model_name:  The name of the model.
        :type model_name:   Str.
        :param doc:         The documentation about the model.
        :type doc:          Str.
        """
        self._pirf = pirf
        self._sirf = sirf
        self._npp0 = npp0
        self._fert = fert
        self._model_name = model_name
        self._doc = doc
        self._prop = Prop.Prop(0, 0e0)
        # Land biosphere C pools (GtC)
        self._mlk = np.array([], dtype=np.longdouble)
        # Mass of C in land biosphere (GtC)
        self._ml = np.array([], dtype=np.longdouble)

        self._exceed = False
        self._co2_exceed = False
        self._t_exceed = False

    def initCarbonStock(self, ntime):
        """
        Initialize the stock of carbon.

        :param ntime:   The time dimension.
        :type ntime:    Int.
        """
        # Initial Land C stock
        self._mlk = np.zeros(self._pirf.getNscale() + 1, dtype=np.longdouble)
        self._ml = np.zeros(ntime, dtype=np.longdouble)
        for i in range(self._pirf.getNscale()):
            self._mlk[i] = (self._npp0 * self._pirf.getWeightElem(i) * self._pirf.getTScaleElem(i))
            self._ml[0] += self._mlk[i]

    def increaseMlkElem(self, index, value):
        """
        Increase mlk at the index by value.

        :param index:   Upper bound for the array.
        :type index:    Int.
        :param value:   New value to add at the mlk array.
        :type value:    Float.
        """
        self._mlk[index] += value

    def setMlElem(self, i, value):
        """
        Setter for the i-th element of ml.

        :param i:       Index.
        :type i:        Int.
        :param value:   Value.
        :type value:    Float.
        """
        self._ml[i] = value

    def co2Exceed(self):
        """
        As the co2 exceeded, we turn the variable to True.
        """
        self._co2_exceed = True

    def tExceed(self):
        """
        As the temperature exceeded, we turn the variable to True.
        """
        self._t_exceed = True

    def getNpp0(self):
        """
        Getter for the npp0 attribute.

        :return:    Return preindustrial net primary production (Gt/yr).
        :rtype:     Float.
        """
        return self._npp0

    def getDoc(self):
        """
        Getter for the documentation of the model.

        :return:    Return the documentation.
        :rtype:     Str.
        """
        return self._doc

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
        Getter for the i-th element of propf.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of propf.
        :rtype:     Float.
        """
        return self._prop.getPropfElem(index)

    def getMlk(self):
        """
        Getter for the mlk attribute (Land biosphere C pools (GtC)).

        :return:    Return the mlk array.
        :rtype:     Numpy array of float.
        """
        return self._mlk

    def getMlElem(self, i):
        """
        Getter for the i-th element of ml.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of ml.
        :rtype:     Float.
        """
        return self._ml[i]

    def isExceed(self):
        """
        Check if the parametrization range of IRF exceed.

        :return:    Check if it exceeded.
        :rtype:     Boolean.
        """
        return self._exceed

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

    def computePropagators(self, dt, pirf=None):
        """
        Calculate propagators.

        :param dt:      Time step.
        :type dt:       Float.
        :param pirf:    IRF coefficients.
        :type pirf:     Pirf.
        """
        if pirf is None:
            self._prop.computePropagators(self._pirf, dt)
        else:
            self._prop.computePropagators(pirf, dt)

    def computePropTempDependent(self, t, dt):
        """
        Wrapper that updates propagators for T-dependent IRF coefficients.

        :param t:   Temperature perturbation.
        :type t:    Float.
        :param dt:  Time step.
        :type dt:   Float.
        """
        p = cp.deepcopy(self._pirf)
        tmp_weights = np.zeros(p.getNscale() + 1, dtype=np.longdouble)
        tmp_t_scales = np.zeros(p.getNscale(), dtype=np.longdouble)
        tmp_sum = 0
        if t > self._sirf.getTMax() and not self._exceed:
            print("warning: temperature parametrization range of IRF ", p.getName(), " exceeded (",
                  self._sirf.getTMax(), "K).", file=sys.stdout)
            self._exceed = True
        for i in range(p.getNscale() + 1):
            tmp_weights[i] = np.exp(self._sirf.getWeightElem(i) * t, dtype=np.longdouble)
        p.multiplyWeights(p.getNscale() + 1, tmp_weights)
        for i in range(p.getNscale() + 1):
            tmp_sum += p.getWeightElem(i)
        p.divideWeights(p.getNscale() + 1, tmp_sum)
        for j in range(p.getNscale()):
            tmp_t_scales[j] = np.exp(-self._sirf.getTscaleElem(j) * t, dtype=np.longdouble)
        p.multiplyTscales(p.getNscale(), tmp_t_scales)
        self.computePropagators(dt, p)

    def computeNpp(self, ma, t, co2_atm0, co2_dep, t_dep, deriv):
        """
        Do the computation of the npp according to the model.

        :param ma:          Mass of C in atmosphere.
        :type ma:           Float.
        :param t:           Global ΔSAT (℃)
        :type t:            Float.
        :param co2_atm0:    Preindustrial equilibrium CO2 concentration.
        :type co2_atm0:     Float.
        :param co2_dep:     Co2 dependence.
        :type co2_dep:      Boolean.
        :param t_dep:       Temperature dependence.
        :type t_dep:        Boolean.
        :param deriv:       Derivative dNPP/dm.
        :type deriv:        Boolean.

        :return:            Return the npp.
        :rtype:             Float.
        """
        if self._model_name == "HRBM":
            return self.computeHrbmNpp(ma, t, co2_dep, t_dep, deriv)
        elif self._model_name == "4box":
            return self.compute4boxNpp(ma, co2_atm0, co2_dep, deriv)

    def computeHrbmNpp(self, ma, t, co2_dep, t_dep, deriv):
        """
        Compute the net primary production for the model HRBM.

        :param ma:      Mass of C in atmosphere.
        :type ma:       Float.
        :param t:       Global ΔSAT (℃)
        :type t:        Float.
        :param co2_dep: Co2 dependence.
        :type co2_dep:  Boolean.
        :param t_dep:   Temperature dependence.
        :type t_dep:    Boolean.
        :param deriv:   Derivative dNPP/dm.
        :type deriv:    Boolean.

        :return:        Return the npp.
        :rtype:         Float.
        """
        if ma / Globals.PPMTOGT > 1120 and not self._co2_exceed:
            print("warning: CO2 parametrization range for NPP in HRBM exceeded (1120ppm)", file=sys.stdout)
            self._co2_exceed = True
        if t > 5 and not self._t_exceed:
            print("warning: temperature parametrization range for NPP in HRBM exceeded (5K)", file=sys.stdout)
            self._t_exceed = True
        if not deriv:
            if co2_dep:
                npp = -np.exp(3.672801e0, dtype=np.longdouble) \
                      + np.exp(-0.430818e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 1, dtype=np.longdouble) \
                      - np.exp(-6.145559e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 2, dtype=np.longdouble) \
                      + np.exp(-12.353878e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 3, dtype=np.longdouble) \
                      - np.exp(-19.010800e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 4, dtype=np.longdouble) \
                      + np.exp(-26.183752e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 5, dtype=np.longdouble) \
                      - np.exp(-34.317488e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 6, dtype=np.longdouble) \
                      - np.exp(-41.553715e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 7, dtype=np.longdouble) \
                      + np.exp(-48.265138e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 8, dtype=np.longdouble) \
                      - np.exp(-56.056095e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 9, dtype=np.longdouble) \
                      + np.exp(-64.818185e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 10, dtype=np.longdouble)
            else:
                npp = self._npp0
            if t_dep:
                npp *= (1 + 0.11780208e+0 * np.tanh(t / 0.509312421e+02,
                                                    dtype=np.longdouble) + 0.24305130e-02 * np.tanh(
                    t / 0.885326739e+01, dtype=np.longdouble))
        else:
            if co2_dep:
                npp_dev = (+ np.exp(-0.430818e0, dtype=np.longdouble)
                           - 2 * np.exp(
                            -6.145559e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 1, dtype=np.longdouble)
                           + 3 * np.exp(
                            -12.353878e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 2, dtype=np.longdouble)
                           - 4 * np.exp(
                            -19.010800e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 3, dtype=np.longdouble)
                           + 5 * np.exp(
                            -26.183752e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 4, dtype=np.longdouble)
                           - 6 * np.exp(
                            -34.317488e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 5, dtype=np.longdouble)
                           - 7 * np.exp(
                            -41.553715e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 6, dtype=np.longdouble)
                           + 8 * np.exp(
                            -48.265138e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 7, dtype=np.longdouble)
                           - 9 * np.exp(
                            -56.056095e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 8, dtype=np.longdouble)
                           + 10 * np.exp(
                            -64.818185e0 + UtilsMath.myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 9, dtype=np.longdouble)
                           ) / Globals.PPMTOGT
            else:
                npp_dev = 0e0
            if t_dep:
                npp_dev *= (1e0 + 0.11780208e+0 * np.tanh(t / 0.509312421e+02, dtype=np.longdouble)
                            + 0.24305130e-02 * np.tanh(t / 0.885326739e+01, dtype=np.longdouble))
            npp = npp_dev
        return npp

    def compute4boxNpp(self, ma, co2_atm0, co2_dep, deriv):
        """
        Compute the net primary production for the 4box model.

        :param ma:          Mass of C in atmosphere.
        :type ma:           Float.
        :param co2_atm0:    Preindustrial equilibrium CO2 concentration.
        :type co2_atm0:     Float.
        :param co2_dep:     Co2 dependence.
        :type co2_dep:      Boolean.
        :param deriv:       Derivative dNPP/dm.
        :type deriv:        Boolean.

        :return:            Return the npp.
        :rtype:             Float.
        """
        if not deriv:
            if co2_dep:
                m0 = co2_atm0 * Globals.PPMTOGT
                d_npp = self._npp0 * self._fert * UtilsMath.myLog(ma / m0)
                npp = d_npp + self._npp0
            else:
                npp = self._npp0
        else:
            if co2_dep:
                npp_dev = self._npp0 * self._fert * (1 / ma)
            else:
                npp_dev = 0e0
            npp = npp_dev
        return npp
