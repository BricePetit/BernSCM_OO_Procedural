"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM in procedural
Python version
"""
import os

import numpy as np

import Globals


def getPath(folder_name):
    """
    Get the complete path and go to the folder corresponding to the folder name.

    :param folder_name: The name of the folder where we want to go.
    :type folder_name:  Str.

    :return:            Return the path.
    :rtype:             Str
    """
    path = os.path.abspath(os.getcwd())[:-6] + folder_name
    if os.name == "nt":
        separator = "\\"
    else:
        separator = "/"
    return path + separator


def readForcing():
    """
    Read forcing file.
    """
    with open(getPath("forcing") + "forcing_" + Globals.scenario + ".dat", "r") as file:
        data = file.readlines()
        for i in range(len(data)):
            splited_tmp = np.array(data[i].split())
            if splited_tmp[0][0] != "#":
                tmp_res = np.array([], dtype=np.longdouble)
                for j in range(Globals.NFORC):
                    if abs(float(splited_tmp[j]) - Globals.NA) < 1e-3:
                        tmp_res = np.append(tmp_res, Globals.NA)
                    elif j == Globals.JACO2:
                        tmp_res = np.append(tmp_res, float(splited_tmp[j]) * Globals.PPMTOGT)
                    else:
                        tmp_res = np.append(tmp_res, float(splited_tmp[j]))
                Globals.forcing = np.append(Globals.forcing, [tmp_res], axis=0)
