"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM in procedural
Python version
"""
from datetime import datetime
import math

import numpy as np

import Globals
import bernSCM_readforcing


def getFileName():
    """
    Compute the name of the file.

    :return:    Return the name of the file where we want to write the solution.
    :rtype:     Str.
    """
    x_imp = ""
    x_lin = ""
    if Globals.IMPLICIT_O:
        if Globals.IMPLICIT_L:
            x_imp = "I"  # Implicit step
        else:
            x_imp = "Io"  # Ocean-only implicit step
    if Globals.LINEAR:
        x_lin = "Q"

    vork = max(1, int(math.log10(1e-5 + Globals.DELTA_T)) + 1)
    nachk = - math.floor(math.log10(1e-5 + Globals.DELTA_T - int(Globals.DELTA_T)))
    if nachk > 4:
        xdt = "{0:{1}d}".format(int(Globals.DELTA_T), int(vork))
    else:
        xdt = "{0:{1}.{2}f}".format(int(Globals.DELTA_T), int(vork + nachk + 1), nachk)

    if Globals.t_dep:
        x_sens = "_t"
    else:
        x_sens = "_t0"

    if Globals.co2_dep:
        x_sens += "_f"
    else:
        x_sens += "_f0"

    vork = int(math.log10(1e-1 + Globals.t2x * 10e0)) + 1
    x_sens += "_CS" + "{0:{1}d}".format(int(Globals.t2x * 10 + 0.5e0), int(vork))
    if Globals.ID is not None:
        Globals.ID = "_" + Globals.ID + "_"
    else:
        Globals.ID = "_"
    return Globals.scenario + "_D" + xdt + x_lin + x_imp + Globals.ID + "BernSCM" + x_sens + ".dat"


def output():
    """
    Function to write the output of our simulation.
    """
    for i in range(1, Globals.ntime):
        # (RF(n)==NA.or.RFnC(n)==NA.or.RFC(n)==NA)
        if Globals.rf[i] == Globals.NA or Globals.rfnc[i] == Globals.NA or Globals.rfc[i] == Globals.NA:
            Globals.rfb[i] = Globals.NA
        else:
            Globals.rfb[i] = Globals.rf[i] - Globals.rfnc[i] - Globals.rfc[i]

    fa = np.zeros(Globals.ntime, dtype=np.longdouble)
    frh = np.full(Globals.ntime, Globals.NPP0, dtype=np.longdouble)
    npp_out = np.full(Globals.ntime, Globals.NPP0, dtype=np.longdouble)
    land_c_uptake = np.zeros(Globals.ntime, dtype=np.longdouble)
    ocean_c_uptake = np.zeros(Globals.ntime, dtype=np.longdouble)
    fossil_emissions = np.zeros(Globals.ntime, dtype=np.longdouble)
    f_deep = np.zeros(Globals.ntime, dtype=np.longdouble)
    midtime = np.zeros(Globals.ntime, dtype=np.longdouble)

    # C budget as in code
    # The following fluxes correspond to midtime (box-centered).
    for i in range(1, Globals.ntime):
        fa[i] = (Globals.ma[i] - Globals.ma[i - 1]) / Globals.DELTA_T
        fossil_emissions[i] = (Globals.e_co2[i] + Globals.e_co2[i - 1]) / 2e0
        if Globals.LINEAR:
            ocean_c_uptake[i] = (Globals.fo[i] + Globals.fo[i - 1]) / 2e0
            npp_out[i] = (Globals.fnpp[i] + Globals.fnpp[i - 1]) / 2e0
        else:
            ocean_c_uptake[i] = Globals.fo[i]
        land_c_uptake[i] = (Globals.ml[i] - Globals.ml[i - 1]) / Globals.DELTA_T
        frh[i] = Globals.fnpp[i] - land_c_uptake[i]
        Globals.fb[i] = fossil_emissions[i] - fa[i] - ocean_c_uptake[i] - land_c_uptake[i]
    if not Globals.LINEAR:
        npp_out = Globals.fnpp

    if Globals.DO_INTERPOL:
        for i in range(Globals.ntime - 1):
            fa[i] = (fa[i] + fa[i + 1]) / 2
            ocean_c_uptake[i] = (ocean_c_uptake[i] + ocean_c_uptake[i + 1]) / 2
            land_c_uptake[i] = (land_c_uptake[i] + land_c_uptake[i + 1]) / 2
            fossil_emissions[i] = (fossil_emissions[i] + fossil_emissions[i + 1]) / 2
        Globals.fb = fossil_emissions - fa - land_c_uptake - ocean_c_uptake
        frh = land_c_uptake - npp_out
    else:
        midtime = Globals.time - Globals.DELTA_T / 2

    desc_file = "# Creation date: " + datetime.now().strftime("%d/%m/%Y %H:%M:%S") \
                + "\n" \
                  "#\n" \
                  "# Bern Simple Climate Model (BernSCM) version 1.0\n" \
                  "# = = = = = = = = = = = = = = = = = = = = = = = =\n" \
                  "#\n" \
                  "# Authors:\n" \
                  "# Kuno Strassmann, Fortunat Joos\n" \
                  "# Climate and Environmental Physics\n" \
                  "# Sidlerstr 5\n" \
                  "# CH-3012 Bern\n" \
                  "# kuno.strassmann@alumni.ethz.ch, joos@climate.unibe.ch\n" \
                  "# Tel.: 0041-31-631 44 61\n" \
                  "#\n" \
                  "# Description:\n" \
                  "#\n" \
                  "# The model is based on impulse response representation (Joos et al., 1996) of the ocean (surface-to-deep\n" \
                  "# transport) and the terrestrial biosphere (accumulation and decay of NPP carbon). Nonlinear ocean surface carbon\n" \
                  "# chemistry and NPP changes are treated by separate parametrizations (see description of C cycle components below).\n" \
                  "#\n" \
                  "# References:\n" \
                  "#\n" \
                  "# - K. Strassmann and F. Joos, 2017: The Bern Simple Climate Model: an extensible and fully\n" \
                  "# documented open source reimplementation of the Bern reduced form model for global carbon\n" \
                  "# cycle-climate simulations. Submitted to Geophysical Model Development.\n" \
                  "# - F. Joos, M. Bruno, R. Fink, T. F. Stocker, U. Siegenthaler, C. Le Quere, and\n" \
                  "# J. L. Sarmiento, 1996: An efficient and accurate representation of complex oceanic\n" \
                  "# and biospheric models of anthropogenic carbon uptake. Tellus, 48B:397-417\n" \
                  "# - Joos, F. and M. Bruno, 1996: Pulse response functions are cost-efficient tools to model the link\n" \
                  "# between carbon emissions, atmospheric CO2 and global warming. Physics and Chemistry of the Earth\n" \
                  "# 21:471-476.\n" \
                  "#\n" \
                  "#\n" \
                  "#\n" \
                  "# Numerical solution\n" \
                  "# = = = = = = = = =\n" \
                  "# Time step: " + str(
        Globals.DELTA_T) + "yr\n"
    if Globals.IMPLICIT_L or Globals.IMPLICIT_O:
        desc_file += "# Implicite step: \n"
        if Globals.IMPLICIT_L:
            desc_file += "# - Land C exchange\n"
        if Globals.IMPLICIT_O:
            desc_file += "# - Ocean C and heat exchange\n"
    if Globals.LINEAR:
        desc_file += "# Discretization: piecewise linear\n"
    else:
        desc_file += "# Discretization: piecewise constant\n"
    desc_file += "#\n" \
                 "#\n" \
                 "#\n" \
                 "# Simulation setup\n" \
                 "# = = = = = = = = =\n" \
                 "#\n" \
                 "# Forcing scenario: " \
                 + Globals.scenario \
                 + "\n" \
                   "#\n" \
                   "# Climate Sensitivity: " \
                 + str(Globals.t2x) \
                 + " degrees C per doubling of atm. CO2 \n" \
                   "#\n" \
                   "# Carbon Cycle:\n" \
                   "# Process sensitivity to atmospheric CO2:\n" \
                   "# - Ocean CO2 uptake (see ocean component)\n"
    if Globals.co2_dep:
        desc_file += "# - Land C exchange (see land component)\n"
    desc_file += "# Process sensitivity to global mean SAT:\n"
    if Globals.t_dep:
        desc_file += "# - Ocean CO2 uptake (see ocean component)\n"
        if Globals.LAND_T_DEP:
            desc_file += "# - Land C exchange (see land component)\n"
        desc_file += "#\n#\n"
    else:
        desc_file += "# none\n#\n#\n"
    desc_file += Globals.OCEAN_DOC
    desc_file += Globals.LAND_DOC

    desc_file += "#time_(yr) glob_temp_dev_(°C) RF_tot_(W/m²) RF_CO2_(W/m²) RF_nonCO2_(W/m²) RF_budget_(W/m²) " \
                 "ocean_heat_uptake_(PW) co2_atm_(ppm) co2_seasurf_(ppm) atm_CO2_change_(GtC/yr) " \
                 "fossil_CO2_em_(GtC/yr) budget_C_uptake_(GtC/yr) ocean_C_uptake_(GtC/yr) land_C_uptake_(GtC/yr) " \
                 "NPP_(GtC/yr) RH_(GtC/yr) LandC_(GtC) dDIC_(µmol/kg) fdeep_(GtC/yr)"
    if not Globals.DO_INTERPOL:
        desc_file += " midtime_(yr)"
    desc_file += "\n"

    desc_file += "#time glob_temp_dev RF_tot RF_CO2 RF_nonCO2 RF_budget ocean_heat_uptake co2_atm co2_seasurf " \
                 "atm_CO2_change fossil_CO2_em budget_C_uptake ocean_C_uptake land_C_uptake NPP RH LandC dDIC fdeep"
    if not Globals.DO_INTERPOL:
        desc_file += " midtime"
    desc_file += "\n"

    for i in range(Globals.ntime):
        desc_file += "{:20.10f}".format(Globals.time[i])
        desc_file += "{:20.10f}".format(Globals.temp[i])
        desc_file += "{:20.10f}".format(Globals.rf[i])
        desc_file += "{:20.10f}".format(Globals.rfc[i])
        desc_file += "{:20.10f}".format(Globals.rfnc[i])
        desc_file += "{:20.10f}".format(Globals.rfb[i])
        desc_file += "{:20.10f}".format(Globals.fh[i])
        desc_file += "{:20.10f}".format(Globals.ma[i] / Globals.PPMTOGT)
        desc_file += "{:20.10f}".format(Globals.co2_atm0 + Globals.dpcs[i])
        desc_file += "{:20.10f}".format(fa[i])
        desc_file += "{:20.10f}".format(fossil_emissions[i])
        desc_file += "{:20.10f}".format(Globals.fb[i])
        desc_file += "{:20.10f}".format(ocean_c_uptake[i])
        desc_file += "{:20.10f}".format(land_c_uptake[i])
        desc_file += "{:20.10f}".format(npp_out[i])
        desc_file += "{:20.10f}".format(frh[i])
        desc_file += "{:20.10f}".format(Globals.ml[i])
        desc_file += "{:20.10f}".format(Globals.ms[i])
        desc_file += "{:20.10f}".format(f_deep[i])
        if not Globals.DO_INTERPOL:
            desc_file += "{:20.10f}".format(midtime[i])
        desc_file += "\n"

    output_path = bernSCM_readforcing.getPath("output")
    file_name = getFileName()
    with open(output_path + file_name, "w+") as file:
        file.write(desc_file)
