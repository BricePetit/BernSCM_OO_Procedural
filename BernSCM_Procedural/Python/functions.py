"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM in procedural
Python version
"""
import sys

import numpy as np

import Globals


def interpol(x, x0, x1, y0, y1):
    """
    Do the interpolation. It is the process of estimating unknown values that 
    fall between known values.

    :param x:   X-value for interpolation.
    :type x:    Float.
    :param x0:  X-value for interpolation.
    :type x0:   Float.
    :param x1:  X-value for interpolation.
    :type x1:   Float.
    :param y0:  Y(x) interpol value.
    :type y0:   Float.
    :param y1:  Y(x) interpol value.
    :type y1:   Float.

    :return:    Return the interpol value.
    :rtype:     Float.
    """
    if y0 == Globals.NA:
        if abs(x - x1) < 1e-9:
            y = y1
        else:
            y = Globals.NA
    elif y1 == Globals.NA:
        if abs(x - x0) < 1e-9:
            y = y0
        else:
            y = Globals.NA
    else:
        y = ((x - x0) * y1 + (x1 - x) * y0) / (x1 - x0)
    return y


def npp(ma, t, deriv):
    """
    Do the computation of the npp according to the model.

    :param ma:          Mass of C in atmosphere.
    :type ma:           Float.
    :param t:           Global ΔSAT (℃)
    :type t:            Float.           
    :param deriv:       Derivative dNPP/dm.
    :type deriv:        Boolean.      

    :return:            Return the npp.
    :rtype:             Float.
    """
    if Globals.LAND_MODEL == "HRBM":
        return computeHrbmNpp(ma, t, deriv)
    elif Globals.LAND_MODEL == "4box":
        return compute4boxNpp(ma, deriv)


def computeHrbmNpp(ma, t, deriv):
    """
    Compute the net primary production for the model HRBM.

    :param ma:      Mass of C in atmosphere.
    :type ma:       Float.
    :param t:       Global ΔSAT (℃)
    :type t:        Float.
    :param deriv:   Derivative dNPP/dm.
    :type deriv:    Boolean.

    :return:        Return the npp.
    :rtype:         Float.
    """
    if ma / Globals.PPMTOGT > 1120 and not Globals.co2_exceed:
        print("warning: CO2 parametrization range for NPP in HRBM exceeded (1120ppm)", file=sys.stdout)
        Globals.co2_exceed = True
    if t > 5 and not Globals.t_exceed:
        print("warning: temperature parametrization range for NPP in HRBM exceeded (5K)", file=sys.stdout)
        Globals.t_exceed = True
    if not deriv:
        if Globals.co2_dep:
            npp = -np.exp(3.672801e0, dtype=np.longdouble) \
                  + np.exp(-0.430818e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 1, dtype=np.longdouble) \
                  - np.exp(-6.145559e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 2, dtype=np.longdouble) \
                  + np.exp(-12.353878e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 3, dtype=np.longdouble) \
                  - np.exp(-19.010800e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 4, dtype=np.longdouble) \
                  + np.exp(-26.183752e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 5, dtype=np.longdouble) \
                  - np.exp(-34.317488e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 6, dtype=np.longdouble) \
                  - np.exp(-41.553715e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 7, dtype=np.longdouble) \
                  + np.exp(-48.265138e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 8, dtype=np.longdouble) \
                  - np.exp(-56.056095e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 9, dtype=np.longdouble) \
                  + np.exp(-64.818185e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 10, dtype=np.longdouble)
        else:
            npp = Globals.NPP0
        if Globals.t_dep:
            npp *= (1 + 0.11780208e+0 * np.tanh(t / 0.509312421e+02, dtype=np.longdouble) + 0.24305130e-02 * np.tanh(
                t / 0.885326739e+01, dtype=np.longdouble))
    else:
        if Globals.co2_dep:
            npp_dev = (+ np.exp(-0.430818e0, dtype=np.longdouble)
                       - 2 * np.exp(
                        -6.145559e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 1, dtype=np.longdouble)
                       + 3 * np.exp(
                        -12.353878e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 2, dtype=np.longdouble)
                       - 4 * np.exp(
                        -19.010800e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 3, dtype=np.longdouble)
                       + 5 * np.exp(
                        -26.183752e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 4, dtype=np.longdouble)
                       - 6 * np.exp(
                        -34.317488e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 5, dtype=np.longdouble)
                       - 7 * np.exp(
                        -41.553715e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 6, dtype=np.longdouble)
                       + 8 * np.exp(
                        -48.265138e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 7, dtype=np.longdouble)
                       - 9 * np.exp(
                        -56.056095e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 8, dtype=np.longdouble)
                       + 10 * np.exp(
                        -64.818185e0 + myLog(min(ma / Globals.PPMTOGT, 1274e0)) * 9, dtype=np.longdouble)
                       ) / Globals.PPMTOGT
        else:
            npp_dev = 0e0
        if Globals.t_dep:
            npp_dev *= (1e0 + 0.11780208e+0 * np.tanh(t / 0.509312421e+02, dtype=np.longdouble)
                        + 0.24305130e-02 * np.tanh(t / 0.885326739e+01, dtype=np.longdouble))
        npp = npp_dev
    return npp


def compute4boxNpp(ma, deriv):
    """
    Compute the net primary production for the 4box model.

    :param ma:          Mass of C in atmosphere.
    :type ma:           Float.
    :param deriv:       Derivative dNPP/dm.
    :type deriv:        Boolean.

    :return:            Return the npp.
    :rtype:             Float.
    """
    if not deriv:
        if Globals.co2_dep:
            m0 = Globals.co2_atm0 * Globals.PPMTOGT
            d_npp = Globals.NPP0 * Globals.FERT * myLog(ma / m0)
            npp = d_npp + Globals.NPP0
        else:
            npp = Globals.NPP0
    else:
        if Globals.co2_dep:
            npp_dev = Globals.NPP0 * Globals.FERT * (1 / ma)
        else:
            npp_dev = 0e0
        npp = npp_dev
    return npp


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


def fasC(ca, dp_co2s):
    """
    Compute the atmosphere-ocean CO2 flux (Gt/yr).

    :param ca:      Atmospheric CO2 concentration (Gt).
    :type ca:       Float.
    :param dp_co2s: Ocean saturation CO2 pressure deviation from preindustrial equilibrium (ppm).
    :type dp_co2s:  Float.

    :return:        Return the atmosphere-ocean CO2 flux (Gt/yr)
    :rtype:         Float.
    """
    return Globals.kg_aoc * ((ca - Globals.co2_atm0 * Globals.PPMTOGT) - dp_co2s * Globals.PPMTOGT)


def fasT(rf_tot, t):
    """
    Compute the air-sea heat flux.

    :param rf_tot:  Radiative forcing (Wm⁻²).
    :type rf_tot:   Float.
    :param t:       Global near-surface atmospheric temperature deviation (℃).
    :type t:        Float.

    :return:        Return the air-sea heat flux (PW).
    :rtype:         Float.
    """
    if Globals.t2x > 0e0:
        return (Globals.aoc / Globals.OFRAC / Globals.PETA) * (rf_tot - (t / Globals.t2x) * Globals.RF2X)
    else:
        # T is always in equilibrium
        return 0e0


def rfCo2(ma):
    """
    Compute the RF of the atmospheric CO2.

    :param ma:  Atmospheric CO₂ (Gt).
    :type ma:   Float.

    :return:    Return RF of atmospheric CO2 (Wm⁻²). 
    :rtype:     Float.
    """
    return Globals.RECO2 * myLog((ma / Globals.PPMTOGT) / Globals.CO2PREIND)


def rfEqC02Ma(rf_co2):
    """
    Calculate equivalent atmospheric CO2 (in GtC) from RF.

    :param rf_co2:  Co2 RF (Wm⁻²).
    :type rf_co2:   Float.

    :return:        Return the atmospheric CO₂ (Gt).
    :rtype:         Float.
    """
    return np.exp(rf_co2 / Globals.RECO2, dtype=np.longdouble) * Globals.CO2PREIND * Globals.PPMTOGT


def computeDpCo2s(d_dic, t, deriv):
    """
    Analytical representation of the zeta-factor for a temperature 
    range from 17.7 to 18.3 degrees celsius (Tchem). Rogers chemistry model
    was used to calculate the zeta-factor.

    :param d_dic:       Change in ocean surface DIC (μmol:kg).
    :type d_dic:        Float.
    :param t:           Global SAT change from preindustrial (℃).
    :type t:            Float.
    :param deriv:       Derivative dpCO2s/ddDIC.
    :type deriv:        Boolean.

    :return:            Ocean saturation CO2 pressure deviation from preindustrial 
                        equilibrium (ppm), or derivative (d dpCs/d dDIC).
    :rtype:             Float.
    """
    if not deriv:
        # Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm)
        dp_co2s = (1.5568e0 - 1.3993e-2 * Globals.t_chem) * d_dic \
                  + (7.4706e0 - 0.20207e0 * Globals.t_chem) * 1e-3 * d_dic ** 2 \
                  - (1.2748e0 - 0.12015e0 * Globals.t_chem) * 1e-5 * d_dic ** 3 \
                  + (2.4491e0 - 0.12639e0 * Globals.t_chem) * 1e-7 * d_dic ** 4 \
                  - (1.5468e0 - 0.15326e0 * Globals.t_chem) * 1e-10 * d_dic ** 5
        if dp_co2s > 1320 and not Globals.co2_exceed:
            print("warning: CO2 parametrization range for dpCO2s exceeded (1320ppm)", file=sys.stdout)
            Globals.co2_exceed = True
        if t + Globals.t_chem > 25 and not Globals.t_exceed:
            print("warning: temperature parametrization range for dpCO2s exceeded (25℃)", file=sys.stdout)
            Globals.t_exceed = True
        if Globals.t_dep:
            dp_co2s = (dp_co2s + Globals.co2_atm0) * np.exp(Globals.buffer_t * t,
                                                            dtype=np.longdouble) - Globals.co2_atm0
    else:
        dp_co2s_dev = (1.5568e0 - 1.3993e-2 * Globals.t_chem) \
                      + 2e0 * (7.4706 - 0.20207e0 * Globals.t_chem) * 1.0e-3 * d_dic \
                      - 3e0 * (1.2748 - 0.12015e0 * Globals.t_chem) * 1.0e-5 * d_dic ** 2 \
                      + 4e0 * (2.4491 - 0.12639e0 * Globals.t_chem) * 1.0e-7 * d_dic ** 3 \
                      - 5e0 * (1.5468 - 0.15326e0 * Globals.t_chem) * 1.0e-10 * d_dic ** 4
        if Globals.t_dep:
            dp_co2s_dev = dp_co2s_dev * np.exp(Globals.buffer_t * t, dtype=np.longdouble)
        dp_co2s = dp_co2s_dev
    return dp_co2s
