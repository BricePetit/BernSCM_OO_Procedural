"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""

import numpy as np

class UtilsMath:
    """
    A class where we will put all util function for math application.
    """

    @staticmethod
    def myLog(value):
        """
        We create our specific log because the owner use a certain flag in order to do the following.

        :param value:   The value to use in the logarithm.
        :type value:    Float.

        :return:        Return the computed value.
        :rtype:         Float.
        """
        if value > 0e0:
            res = np.log(value, dtype=np.longdouble)
        else:
            res = value
        return res

