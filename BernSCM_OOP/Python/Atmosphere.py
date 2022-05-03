"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""
import numpy as np

import Globals
from UtilsMath import UtilsMath


class Atmosphere:
    """
    This class represents the atmosphere.
    """

    def __init__(self):
        # Mass of C in atmosphere (GtC)
        self._ma = np.array([], dtype=np.longdouble)
        # CO₂ emissions (GtC/yr)
        self._e_co2 = np.array([], dtype=np.longdouble)
        # Preindustrial equilibrium CO2 concentration
        self._co2_atm0 = 0e0

    def initialize(self, ntime):
        """
        Initialize the mass of carbon in atmosphere according to the time dimension.

        :param ntime:   Time dimension.
        :type ntime:    Int.
        """
        self._ma = np.zeros(ntime, dtype=np.longdouble)
        self._e_co2 = np.zeros(ntime, dtype=np.longdouble)

    def initializeCo2Atm0(self):
        """
        Initialize the preindustrial equilibrium CO2 concentration.
        """
        self._co2_atm0 = self._ma[0] / Globals.PPMTOGT

    def setMaElem(self, i, value):
        """
        Setter for the i-th element of ma.

        :param i:       Index.
        :type i:        Int.
        :param value:   Value.
        :type value:    Float.
        """
        self._ma[i] = value

    def setECo2Elem(self, i, value):
        """
        Setter for the i-th element of e_co2.

        :param i:       Index.
        :type i:        Int.
        :param value:   Value.
        :type value:    Float.
        """
        self._e_co2[i] = value

    def getMaElem(self, i):
        """
        Getter for the i-th element of ma (mass of carbon in the atmosphere).

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of ma.
        :rtype:     Float.
        """
        return self._ma[i]

    def getECo2Elem(self, i):
        """
        Getter for the i-th element of eCo2 (CO₂ emissions (GtC/yr)).

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of e_co2.
        :rtype:     Float.
        """
        return self._e_co2[i]

    def getCo2Atm0(self):
        """
        Getter for the preindustrial equilibrium CO2 concentration.

        :return:    Return the co2_atm0 attribute.
        :rtype:     Float.
        """
        return self._co2_atm0

    def rfCo2(self, i):
        """
        Compute the RF of the atmospheric CO2.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the RF (radiation forcing) of atmospheric CO2 (Wm⁻²).
        :rtype:     Float.
        """
        return Globals.RECO2 * UtilsMath.myLog((self._ma[i] / Globals.PPMTOGT) / Globals.CO2PREIND)
