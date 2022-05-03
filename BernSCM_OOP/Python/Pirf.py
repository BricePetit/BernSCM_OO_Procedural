"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""


class Pirf:
    """
    This class represents IRF parameter set.
    """

    def __init__(self, name, nscale, weight, t_scale):
        """
        Constructor of the PIRF.

        :param name:    Name of the IRF.
        :type name:     Str.
        :param nscale:  Number of finite timescales.
        :type nscale:   Int.
        :param weight:  Numpy array of weight or fraction of input for each box.
        :type weight:   Numpy array of float.
        :param t_scale: Numpy array of time scales of each box (yr).
        :type t_scale:  Numpy array of float.
        """
        self._name = name
        self._nscale = nscale
        self._weight = weight
        self._t_scale = t_scale

    def getName(self):
        """
        Getter for the name of the PIRF.

        :return:    Return the name of the PIRF.
        :rtype:     Str.
        """
        return self._name

    def getNscale(self):
        """
        Getter for the nscale value.

        :return:    Return nscale.
        :rtype:     Int.
        """
        return self._nscale

    def getWeightElem(self, i):
        """
        Getter for the i-th element of weight.

        :param i:   The element i of the weight.
        :type i:    Int.

        :return:    Return the i-th weight value.
        :rtype:     Float.
        """
        return self._weight[i]

    def getTScaleElem(self, i):
        """
        Getter for the i-th element of t_scale.

        :param i:   The element i of the t_scale.
        :type i:    Int.

        :return:    Return the i-th t_scale value.
        :rtype:     Float.
        """
        return self._t_scale[i]

    def multiplyWeights(self, bound, values):
        """
        Multiply weights in the list starting from the first element up to the bound by the values
        in the array.

        :param bound:   The index used to define the number of value that we want to multiply.
        :type bound:    Integer.
        :param value:   Values that we want to use to multiply the weight.
        :type value:    Numpy array of float.
        """
        for i in range(bound):
            self._weight[i] *= values[i]

    def divideWeights(self, bound, value):
        """
        Divide weights in the list starting from the first element up to the bound by the value.

        :param bound:   Upper bound of the list.
        :type bound:    Int.
        :param value:   Value by which we will divide weights.
        :type value:    Float.
        """
        for i in range(bound):
            self._weight[i] /= value

    def multiplyTscales(self, bound, values):
        """
        Multiply t_scale in the list starting from the first element up to the index i by the by the values
        in the array.

        :param bound:   Upper bound of the list.
        :type bound:    Int.
        :param value:   Values by which we will multiply t_scale values.
        :type value:    Numpy array of float.
        """
        for i in range(bound):
            self._t_scale[i] *= values[i]
