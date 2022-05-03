"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""
import math
import os
from datetime import datetime

import numpy as np

import Globals


class UtilsFile:
    """
    A class where we will put all util function.
    """

    def __init__(self, ID):
        """
        Init this class with utils values.

        :param ID:  Additional simulation identifier.
        :type ID:   Str.
        """
        self._id = ID
        self._do_interpol = 1

    def getPath(self, folder_name):
        """
        Get the complete path and go to the folder corresponding to the folder name.

        :param folder_name: The name of the folder where we want to go.
        :type folder_name:  Str.

        :return:            Return the path.
        :rtype:             Str.
        """
        path = os.path.abspath(os.getcwd())[:-6] + folder_name
        if os.name == "nt":
            separator = "\\"
        else:
            separator = "/"
        return path + separator

    def getFileName(self, t2x, t_dep, co2_dep, scenario, dt):
        """
        Compute the name of the file.

        :param t2x:         Climate sensitivity.
        :type t2x:          Float.
        :param t_dep:       Temperature dependence.
        :type t_dep:        Boolean.
        :param co2_dep:     Co2 dependence.
        :type co2_dep:      Boolean.
        :param scenario:    The name of the scenario file without forcing_ (at the beginning) and .dat (at the end).
        :type scenario:     Str.
        :param dt:          Time step.
        :type dt:           Float.

        :return:            Return the name of the file where we want to write the solution.
        :rtype:             Str.
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

        vork = max(1, int(math.log10(1e-5 + dt)) + 1)
        nachk = - math.floor(math.log10(1e-5 + dt - int(dt)))
        if nachk > 4:
            xdt = "{0:{1}d}".format(int(dt), int(vork))
        else:
            xdt = "{0:{1}.{2}f}".format(int(dt), int(vork + nachk + 1), nachk)

        if t_dep:
            x_sens = "_t"
        else:
            x_sens = "_t0"

        if co2_dep:
            x_sens += "_f"
        else:
            x_sens += "_f0"

        vork = int(math.log10(1e-1 + t2x * 10e0)) + 1
        x_sens += "_CS" + "{0:{1}d}".format(int(t2x * 10 + 0.5e0), int(vork))
        if self._id is not None:
            self._id = "_" + self._id + "_"
        else:
            self._id = "_"
        return scenario + "_D" + xdt + x_lin + x_imp + self._id + "BernSCM" + x_sens + ".dat"

    def readForcing(self, scenario):
        """
        Read forcing file.

        :param scenario:    The name of the scenario file without forcing_ (at the beginning) and .dat (at the end).
        :type scenario:     Str.

        :return:            Return all forcing values.
        :rtype:             Numpy array of float.
        """
        forcing = np.empty([0, 7], dtype=np.longdouble)
        with open(self.getPath("forcing") + "forcing_" + scenario + ".dat", "r") as file:
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
                    forcing = np.append(forcing, [tmp_res], axis=0)
                    # forcing is a matrix where the column represents:
                    # [Year, glob_temp_dev, RF_nonCO2, RF_budget, co2_atm, fossil_CO2_em, budget_C_uptake]
        return forcing

    def output(self, earth):
        """
        Function to write the output of our simulation.

        :param earth:  Object Earth containing the atmosphere, land and ocean model.
        :type earth:   Earth.
        """
        for i in range(1, earth.getNtime()):
            if earth.getRfElem(i) == Globals.NA or earth.getRfncElem(i) == Globals.NA or earth.getRfcElem(
                    i) == Globals.NA:
                earth.setRfbElem(i, Globals.NA)
            else:
                earth.setRfbElem(i, earth.getRfElem(i) - earth.getRfncElem(i) - earth.getRfcElem(i))

        fa = np.zeros(earth.getNtime(), dtype=np.longdouble)
        frh = np.full(earth.getNtime(), earth.getNpp0(), dtype=np.longdouble)
        npp_out = np.full(earth.getNtime(), earth.getNpp0(), dtype=np.longdouble)
        land_c_uptake = np.zeros(earth.getNtime(), dtype=np.longdouble)
        ocean_c_uptake = np.zeros(earth.getNtime(), dtype=np.longdouble)
        fossil_emissions = np.zeros(earth.getNtime(), dtype=np.longdouble)
        f_deep = np.zeros(earth.getNtime(), dtype=np.longdouble)
        midtime = np.zeros(earth.getNtime(), dtype=np.longdouble)

        # C budget as in code
        # The following fluxes correspond to midtime (box-centered).
        for i in range(1, earth.getNtime()):
            fa[i] = (earth.getMaElem(i) - earth.getMaElem(i - 1)) / earth.getDt()
            fossil_emissions[i] = (earth.getECo2Elem(i) + earth.getECo2Elem(i - 1)) / 2e0
            if Globals.LINEAR:
                ocean_c_uptake[i] = (earth.getFoElem(i) + earth.getFoElem(i - 1)) / 2e0
                npp_out[i] = (earth.getFnppElem(i) + earth.getFnppElem(i - 1)) / 2e0
            else:
                ocean_c_uptake[i] = earth.getFoElem(i)
            land_c_uptake[i] = (earth.getMlElem(i) - earth.getMlElem(i - 1)) / earth.getDt()
            frh[i] = earth.getFnppElem(i) - land_c_uptake[i]
            earth.setFbElem(i, fossil_emissions[i] - fa[i] - ocean_c_uptake[i] - land_c_uptake[i])
        if not Globals.LINEAR:
            npp_out = earth.getFnpp()

        if self._do_interpol:
            for i in range(earth.getNtime() - 1):
                fa[i] = (fa[i] + fa[i + 1]) / 2
                ocean_c_uptake[i] = (ocean_c_uptake[i] + ocean_c_uptake[i + 1]) / 2
                land_c_uptake[i] = (land_c_uptake[i] + land_c_uptake[i + 1]) / 2
                fossil_emissions[i] = (fossil_emissions[i] + fossil_emissions[i + 1]) / 2

                earth.setFbElem(i, fossil_emissions[i] - fa[i] - land_c_uptake[i] - ocean_c_uptake[i])
                frh[i] = land_c_uptake[i] - npp_out[i]
            tmp_index = earth.getNtime() - 1
            earth.setFbElem(tmp_index, fossil_emissions[tmp_index] - fa[tmp_index] - land_c_uptake[tmp_index] - 
                            ocean_c_uptake[tmp_index])
            frh[tmp_index] = land_c_uptake[tmp_index] - npp_out[tmp_index]
        else:
            for i in range(earth.getNtime()):
                midtime[i] = earth.getTimeElem(i) - earth.getDt() / 2

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
                      "# Time step: " + str(earth.getDt()) + "yr\n"
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
                     + earth.getScenario() \
                     + "\n" \
                       "#\n" \
                       "# Climate Sensitivity: " \
                     + str(earth.getT2x()) \
                     + " degrees C per doubling of atm. CO2 \n" \
                       "#\n" \
                       "# Carbon Cycle:\n" \
                       "# Process sensitivity to atmospheric CO2:\n" \
                       "# - Ocean CO2 uptake (see ocean component)\n"
        if earth.isCo2Dep():
            desc_file += "# - Land C exchange (see land component)\n"
        desc_file += "# Process sensitivity to global mean SAT:\n"
        if earth.isTDep():
            desc_file += "# - Ocean CO2 uptake (see ocean component)\n"
            if Globals.LAND_T_DEP:
                desc_file += "# - Land C exchange (see land component)\n"
            desc_file += "#\n#\n"
        else:
            desc_file += "# none\n#\n#\n"
        desc_file += earth.getOceanDoc()
        desc_file += earth.getLandDoc()

        desc_file += "#time_(yr) glob_temp_dev_(°C) RF_tot_(W/m²) RF_CO2_(W/m²) RF_nonCO2_(W/m²) RF_budget_(W/m²) " \
                     "ocean_heat_uptake_(PW) co2_atm_(ppm) co2_seasurf_(ppm) atm_CO2_change_(GtC/yr) " \
                     "fossil_CO2_em_(GtC/yr) budget_C_uptake_(GtC/yr) ocean_C_uptake_(GtC/yr) land_C_uptake_(GtC/yr) " \
                     "NPP_(GtC/yr) RH_(GtC/yr) LandC_(GtC) dDIC_(µmol/kg) fdeep_(GtC/yr)"
        if not self._do_interpol:
            desc_file += " midtime_(yr)"

        desc_file += "\n#time glob_temp_dev RF_tot RF_CO2 RF_nonCO2 RF_budget ocean_heat_uptake co2_atm co2_seasurf " \
                     "atm_CO2_change fossil_CO2_em budget_C_uptake ocean_C_uptake land_C_uptake NPP RH LandC dDIC fdeep"
        if not self._do_interpol:
            desc_file += " midtime"
        desc_file += "\n"

        for i in range(earth.getNtime()):
            desc_file += "{:20.10f}".format(earth.getTimeElem(i))
            desc_file += "{:20.10f}".format(earth.getTempElem(i))
            desc_file += "{:20.10f}".format(earth.getRfElem(i))
            desc_file += "{:20.10f}".format(earth.getRfcElem(i))
            desc_file += "{:20.10f}".format(earth.getRfncElem(i))
            desc_file += "{:20.10f}".format(earth.getRfbElem(i))
            desc_file += "{:20.10f}".format(earth.getFhElem(i))
            desc_file += "{:20.10f}".format(earth.getMaElem(i) / Globals.PPMTOGT)
            desc_file += "{:20.10f}".format(earth.getCo2Atm0() + earth.getDpcsElem(i))
            desc_file += "{:20.10f}".format(fa[i])
            desc_file += "{:20.10f}".format(fossil_emissions[i])
            desc_file += "{:20.10f}".format(earth.getFbElem(i))
            desc_file += "{:20.10f}".format(ocean_c_uptake[i])
            desc_file += "{:20.10f}".format(land_c_uptake[i])
            desc_file += "{:20.10f}".format(npp_out[i])
            desc_file += "{:20.10f}".format(frh[i])
            desc_file += "{:20.10f}".format(earth.getMlElem(i))
            desc_file += "{:20.10f}".format(earth.getMsElem(i))
            desc_file += "{:20.10f}".format(f_deep[i])
            if not self._do_interpol:
                desc_file += "{:20.10f}".format(midtime[i])
            desc_file += "\n"

        output_path = self.getPath("output")
        file_name = self.getFileName(earth.getT2x(), earth.isTDep(), earth.isCo2Dep(), earth.getScenario(),
                                     earth.getDt())
        with open(output_path + file_name, "w+") as file:
            file.write(desc_file)
