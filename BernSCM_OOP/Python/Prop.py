"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""
import numpy as np

import Globals


class Prop:
    """
    This class represents propagators.
    """

    def __init__(self, nscale, x):
        """
        Constructor of the propagator.

        :param nscale:  Number of finite timescales.
        :type nscale:   Int.
        :param x:       X value.
        :type x:        Float.
        """
        self._nscale = nscale
        self._x = x
        self._propm = np.zeros(Globals.NSCALEMAX + 1, dtype=np.longdouble)
        self._propf = np.zeros(Globals.NSCALEMAX + 1, dtype=np.longdouble)
        self._propfo = np.zeros(Globals.NSCALEMAX + 1, dtype=np.longdouble)

    def getNscale(self):
        """
        Getter for nscale value.

        :return:    Return the nscale value.
        :rtype:     Int.
        """
        return self._nscale

    def getPropmElem(self, i):
        """
        Getter for the i-th element of propm.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of propm.
        :rtype:     Float.
        """
        return self._propm[i]

    def getPropfElem(self, i):
        """
        Getter for the i-th element of propf.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of propf.
        :rtype:     Float.
        """
        return self._propf[i]

    def getPropfoElem(self, i):
        """
        Getter for the i-th element of propfo.

        :param i:   Index.
        :type i:    Int.

        :return:    Return the i-th value of propfo.
        :rtype:     Float.
        """
        return self._propfo[i]

    def computePropagators(self, pirf, dt):
        """
        Calculate propagators.

        :param pirf:    IRF coefficients.
        :type pirf:     Pirf.
        :param dt:      Time step.
        :type dt:       Float.
        """
        if not Globals.LINEAR:
            for i in range(pirf.getNscale()):
                if pirf.getTScaleElem(i) <= Globals.TEQUIL:
                    self._propf[i] = pirf.getWeightElem(i) * pirf.getTScaleElem(i)
                    self._propm[i] = 0e0
                else:
                    self._propf[i] = pirf.getWeightElem(i) * pirf.getTScaleElem(i) \
                                     * (1e0 - np.exp(-dt / pirf.getTScaleElem(i), dtype=np.longdouble))
                    self._propm[i] = np.exp(-dt / pirf.getTScaleElem(i), dtype=np.longdouble)
            self._propf[pirf.getNscale()] = dt * pirf.getWeightElem(pirf.getNscale())
            self._propm[pirf.getNscale()] = 1e0
        elif Globals.LINEAR:
            for i in range(pirf.getNscale()):
                self._propf[i] = pirf.getWeightElem(i) * pirf.getTScaleElem(i) \
                                 * (dt - pirf.getTScaleElem(i)
                                    * (1e0 - np.exp(-dt / pirf.getTScaleElem(i), dtype=np.longdouble))) / dt
                self._propfo[i] = pirf.getWeightElem(i) * pirf.getTScaleElem(i) * (
                        1e0 - np.exp(-dt / pirf.getTScaleElem(i), dtype=np.longdouble)) - self._propf[i]
            self._propf[pirf.getNscale()] = pirf.getWeightElem(pirf.getNscale()) * dt / 2e0
            self._propfo[pirf.getNscale()] = pirf.getWeightElem(pirf.getNscale()) * dt \
                - self._propf[pirf.getNscale()]

            for i in range(pirf.getNscale()):
                self._propm[i] = np.exp(-dt / pirf.getTScaleElem(i), dtype=np.longdouble)
            self._propm[pirf.getNscale()] = 1e0
        self._nscale = pirf.getNscale()
