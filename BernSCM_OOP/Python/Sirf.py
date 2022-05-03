"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""


class Sirf:
    """
    This class represents IRF temperature sensitivities.
    """

    def __init__(self, t_max, weight, t_scale):
        """
        Constructor of the SIRF.

        :param t_max:   Parametrization range.
        :type t_max:    Float.
        :param weight:  Temperature sensitivity of weights.
        :type weight:   Numpy array of float.
        :param t_scale: Temperature sensitivity of tscales.
        :type t_scale:  Numpy array of float.
        """
        self._t_max = t_max
        self._weight = weight
        self._t_scale = t_scale

    def getTMax(self):
        """
        Return the t_max (Temperature max).

        :return:    Return the t_max value.
        :rtype:     Float.
        """
        return self._t_max

    def getWeightElem(self, index):
        """
        Getter for the i-th element of the weight array.

        :param index:   Index.
        :type index:    Int.

        :return:        The i-th value of the array of weights.
        :rtype:         Float
        """
        return self._weight[index]

    def getTscaleElem(self, index):
        """
        Getter for the i-th element of the t_scale array.

        :param index:   Index.
        :type index:    Int.

        :return:        The i-th value of the array of t_scale.
        :rtype:         Float
        """
        return self._t_scale[index]
