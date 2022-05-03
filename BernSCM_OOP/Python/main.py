"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""
import os
import time

import numpy as np

import Atmosphere
import Earth
import Globals
import Land
import Ocean
import Pirf
import Sirf
import UtilsFile


def createOceanModel():
    """
    Function to create the Ocean object.

    :return:    Return the ocean model.
    :rtype:     Ocean.
    """
    if Globals.OCEAN_MODEL == "HILDA":
        ocean_weight = np.array(
            [0.27830433e0, 0.23337218e0, 0.13732822e0, 0.05154051e0, 0.03503318e0, 0.24013944e0, .022936e0],
            dtype=np.longdouble)
        ocean_t_scale = np.array(
            [0.45253504e0, 2.19901724e0, 12.03837102e0, 59.58359820e0, 237.30651757e0, 0.03855458e0, Globals.NA],
            dtype=np.longdouble)
        ocean_pirf = Pirf.Pirf(Globals.OCEAN_MODEL, 6, ocean_weight, ocean_t_scale)
        doc = "# HILDA ocean mixed layer pulse response\n" \
              "# = = = = = = = = = = = = = = = = = = = =\n" \
              "#\n" \
              "# Description:\n" \
              "#\n" \
              "# The near-linear transport of heat and carbon content perturbations from the preindustrial\n" \
              "# equilibrium of the HILDA box-diffusion ocean model is captured by a mixed-layer pulse\n" \
              "# response function (Joos 1996). HILDA is calibrated using observed radiocarbon tracers.\n" \
              "#\n" \
              "# References:\n" \
              "#\n" \
              "# - Joos, F., 1992. Modellierung der Verteilung von Spurenstoffen im Ozean und\n" \
              "# des globalen Kohlenstoffkreislaufes. Phd thesis, Universität Bern.\n" \
              "# - Siegenthaler, U. and Joos, F., 1992. Use of a simple model for studying\n" \
              "# oceanic tracer distributions and the global carbon cycle. Tellus B 44B, 186-207.\n" \
              "# - F. Joos, M. Bruno, R. Fink, T. F. Stocker, U. Siegenthaler, C. Le Quere, and\n" \
              "# J. L. Sarmiento, 1996: An efficient and accurate representation of complex oceanic\n" \
              "# and biospheric models of anthropogenic carbon uptake. Tellus, 48B:397-417\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "# Ocean chemistry parametrization\n" \
              "# = = = = = = = = = = = = = = = = =\n" \
              "#\n" \
              "# Description: \n" \
              "#\n" \
              "# Nonlinear dependence of buffer factor on DIC after Joos et al. (1996).\n" \
              "# Impact of sea surface warming on carbonate chemistry after Takahashi (1993).\n" \
              "#\n" \
              "# References:\n" \
              "#\n" \
              "# - F. Joos, M. Bruno, R. Fink, T. F. Stocker, U. Siegenthaler, C. Le Quere, and\n" \
              "# J. L. Sarmiento, 1996: An efficient and accurate representation of complex oceanic\n" \
              "# and biospheric models of anthropogenic carbon uptake. Tellus, 48B:397-417\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "#\n"
        ocean = Ocean.Ocean(ocean_pirf, 75e0, 4000e0, 1028e0, 1026.5e0, 3.62e14, 1e0 / 9.06e0, 18.1716e0, 
                            Globals.OCEAN_MODEL, doc)
    elif Globals.OCEAN_MODEL == "B2D":
        ocean_weight = np.array(
            [0.09467125e0, 0.1029231e0, 0.03928349e0, 0.4593721e0, 0.0129862e0, 0.2702249, 1.3691e-2],
            dtype=np.longdouble)
        ocean_t_scale = np.array(
            [2.690038e0, 13.61728e0, 86.79685e0, 0.5762091e0, 337.2983e0, 0.07027151e0, Globals.NA],
            dtype=np.longdouble)
        ocean_pirf = Pirf.Pirf(Globals.OCEAN_MODEL, 6, ocean_weight, ocean_t_scale)
        doc = "# Bern2.5D ocean mixed layer pulse response\n" \
              "# = = = = = = = = = = = = = = = = = = = = = =\n" \
              "#\n" \
              "# Description:\n" \
              "#\n" \
              "# The near-linear transport of heat and carbon content perturbations from the preindustrial\n" \
              "# equilibrium of the Bern2.5D ocean model (Stocker et al., 1992) is captured by a mixed-layer\n" \
              "# pulse-response function (Joos 1996).\n" \
              "#\n" \
              "# References:\n" \
              "#\n" \
              "# - T.F. Stocker, D.G. Wright, L.A. Mysak, 1992: A zonally averaged, coupled ocean-atmosphere\n" \
              "# model for paleoclimate studies. J. Clim. 5, 773-797.\n" \
              "# - F. Joos, M. Bruno, R. Fink, T. F. Stocker, U. Siegenthaler, C. Le Quere, and\n" \
              "# J. L. Sarmiento, 1996: An efficient and accurate representation of complex oceanic\n" \
              "# and biospheric models of anthropogenic carbon uptake. Tellus, 48B:397-417\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "# Ocean chemistry parametrization\n" \
              "# = = = = = = = = = = = = = = = = =\n" \
              "#\n" \
              "# Description: \n" \
              "#\n" \
              "# Nonlinear dependence of buffer factor on DIC after Joos et al. (1996).\n" \
              "# Impact of sea surface warming on carbonate chemistry after Takahashi (1993).\n" \
              "#\n" \
              "# References:\n" \
              "#\n" \
              "# - F. Joos, M. Bruno, R. Fink, T. F. Stocker, U. Siegenthaler, C. Le Quere, and\n" \
              "# J. L. Sarmiento, 1996: An efficient and accurate representation of complex oceanic\n" \
              "# and biospheric models of anthropogenic carbon uptake. Tellus, 48B:397-417\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "#\n"
        ocean = Ocean.Ocean(ocean_pirf, 50e0, 4000e0, 1028e0, 1026.5e0, 3.5375e14, 1e0 / 7.46e0, 18.2997e0,
                            Globals.OCEAN_MODEL, doc)
    elif Globals.OCEAN_MODEL == "Princeton3D":
        ocean_weight = np.array([2.27446514e0, 0.06161763e0, 0.03726494e0, 1.28186186e0, 0.01956537e0, -2.70925536e0,
                                 0.01481883e0], dtype=np.longdouble)
        ocean_t_scale = np.array(
            [1.19761983e0, 16.67585709e0, 65.10188851e0, 2.00904478e0, 347.58378832e0, 1.55213441e0, Globals.NA],
            dtype=np.longdouble)
        ocean_pirf = Pirf.Pirf(Globals.OCEAN_MODEL, 6, ocean_weight, ocean_t_scale)
        doc = "# Princeton AOGCM ocean effective mixed layer pulse response\n" \
              "# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n" \
              "#\n" \
              "# Description:\n" \
              "#\n" \
              "# An effective mixed-layer pulse response is derived by deconvolving an atmospheric pulse response\n" \
              "# experiment with the Princeton AOGCM (pulse size 265 ppm), approximating the effect of nonlinearities\n" \
              "# arising from the spatial variability of air-sea exchange and local transport. Representation error\n" \
              "# in cumulative ocean uptake is <=6% for stabilization scenarios S450 and S750 (better for future\n" \
              "# projection.\n" \
              "#\n" \
              "# References:\n" \
              "#\n" \
              "# - Sarmiento, J., Orr, J. and Siegenthaler, U., 1992. A perturbation simulation of CO2 uptake in an ocean\n" \
              "# general circulation model. Journal of Geophysical Research 97, 3621-3645. doi:10.1029/91JC02849.\n" \
              "# - F. Joos, M. Bruno, R. Fink, T. F. Stocker, U. Siegenthaler, C. Le Quere, and\n" \
              "# J. L. Sarmiento, 1996: An efficient and accurate representation of complex oceanic\n" \
              "# and biospheric models of anthropogenic carbon uptake. Tellus, 48B:397-417\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "# Ocean chemistry parametrization\n" \
              "# = = = = = = = = = = = = = = = = = \n" \
              "#\n" \
              "# Description: \n" \
              "#\n" \
              "# Nonlinear dependence of buffer factor on DIC after Joos et al. (1996).\n" \
              "# Impact of sea surface warming on carbonate chemistry after Takahashi (1993).\n" \
              "#\n" \
              "# References:\n" \
              "#\n" \
              "# - F. Joos, M. Bruno, R. Fink, T. F. Stocker, U. Siegenthaler, C. Le Quere, and\n" \
              "# J. L. Sarmiento, 1996: An efficient and accurate representation of complex oceanic\n" \
              "# and biospheric models of anthropogenic carbon uptake. Tellus, 48B:397-417\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "#\n"
        ocean = Ocean.Ocean(ocean_pirf, 50.9e0, 4000e0, 1028e0, 1026.5e0, 3.55e14, 1e0 / 7.66e0, 17.7e0, 
                            Globals.OCEAN_MODEL, doc)
    else:
        raise NameError("Unknown ocean model.")
    return ocean


def createLandModel():
    """
    Function to create Land object.

    :return:    Return the land model.
    :rtype:     Land.
    """
    if Globals.LAND_MODEL == "HRBM":
        land_weight = np.array([-0.15431665e0, 0.56172938e0, 0.07486977e0, 0.41365999e0, 0.10405757e0, 0e0, Globals.NA],
                               dtype=np.longdouble)
        land_t_scale = np.array(
            [0.2010730e0, 1.4753908e0, 8.889799331e0, 74.098061812e0, 2.538054e02, Globals.NA, Globals.NA],
            dtype=np.longdouble)
        land_pirf = Pirf.Pirf(Globals.LAND_MODEL, 5, land_weight, land_t_scale)
        land_weight_sirf = np.array([0.14e0, 0.056e0, 0.072e0, 0.044e0, 0.069e0, 0e0, Globals.NA], dtype=np.longdouble)
        land_t_scale_sirf = np.array([0.056e0, 0.079e0, 0.057e0, 0.053e0, 0.036e0, Globals.NA, Globals.NA],
                                dtype=np.longdouble)
        land_sirf = Sirf.Sirf(5e0, land_weight_sirf, land_t_scale_sirf)
        doc = "# HRBM pulse response/multibox substitute model\n" \
              "# = = = = = = = = = = = = = = = = = = = = = = = =\n" \
              "#\n" \
              "# Description:\n" \
              "#\n" \
              "# The fate of assimilated C as simulated by HRBM is represented by several pools.\n" \
              "# Each pool receives a fraction of global NPP as input, and loses C by heterotrophic\n" \
              "# respiration according to a characteristic turnover time. Both input fractions and\n" \
              "# Turnover times are dependent on global mean SAT (dynamic vegetation changes are\n" \
              "# not captured.).\n" \
              "# The box coefficients are derived from a decay-IRF by integrating to an IRF\n" \
              "# for the current stored fraction (as in ocean IRFs).\n" \
              "#\n" \
              "# References:\n" \
              "#\n" \
              "# - Meyer et al. 1999: The substitution of high-resolution terrestrial biosphere\n" \
              "# models and carbon sequestration in response to changing CO2 and climate,\n" \
              "# Global Biogeochemical Cycles 13, 785)\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "# HRBM NPP parametrization\n" \
              "# = = = = = = = = = = = = =\n" \
              "#\n" \
              "# Description:\n" \
              "#\n" \
              "# NPP is parametrized as product of two functions fitted to\n" \
              "# a series of HRBM simulations, one for temperature dependence,\n" \
              "# and one for CO2 dependence.\n" \
              "#\n" \
              "#\n"
        land = Land.Land(land_pirf, land_sirf, 41.73282e0, None, Globals.LAND_MODEL, doc)
    elif Globals.LAND_MODEL == "4box":
        land_weight = np.array([2.006028671730e0, 0.26828e0, -1.56754896028e0, 0.29323e0, 0e0, Globals.NA, Globals.NA],
                               dtype=np.longdouble)
        land_t_scale = np.array([2.857143e0, 20e0, 2.181818e0, 100e0, Globals.NA, Globals.NA, Globals.NA],
                                dtype=np.longdouble)
        land_pirf = Pirf.Pirf("4box", 4, land_weight, land_t_scale)
        doc = "# 4-box Biosphere of Bern Model\n" \
              "# = = = = = = = = = = = = = = = = \n" \
              "#\n" \
              "# Description:\n" \
              "#\n" \
              "# Global box model consisting of four pools representing ground vegetation,\n" \
              "# wood, detritus and soil organic carbon, with each reservoir having a distinct\n" \
              "# turnover time independent of environmental conditions (SAT).\n" \
              "#\n" \
              "# References:\n" \
              "#\n" \
              "# Ref: F. Joos, M. Bruno, R. Fink, T. F. Stocker,\n" \
              "# U. Siegenthaler and F. Joos, 1992: Use of a simple model\n" \
              "# for studying oceanic tracer distributions and the global carbon cycle,\n" \
              "# Tellus, Ser. B, 44, 186-207, 1992.\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "#\n" \
              "# NPP function for Bern 4 Box biosphere\n" \
              "# = = = = = = = = = = = = = = = = = = = =\n" \
              "#\n" \
              "# Description:\n" \
              "#\n" \
              "# NPP depends logarithmically on atmospheric CO2 (Enting et al., 1994) with the\n" \
              "# scale factor beta=0.287 determined from C budget closure (Schimel et al. 1996) for 1980-1990.\n" \
              "#\n" \
              "# References:\n" \
              "#\n" \
              "# - I.G. Enting, T.M.L. Wigley, and M. Heimann, 1994: Future emissions and concentrations of\n" \
              "# carbon dioxide: Key ocean/atmosphere/land analyses, Tech. Rep. 31, CSIRO\n" \
              "# Atmos. Res., Melbourne, Victoria.\n" \
              "# - D.J. Schimel, D. Alves, I.G. Enting, M. Heimann, F. Joos, D. Raynaud, T. Wigley, 1996:\n" \
              "# CO2 and the carbon cycle, in IPCC Second Scientific Assessment of Climate Change,\n" \
              "# ed. by Houghton, pp. 76-86, Cambridge Univ. Press, NY.\n" \
              "#\n" \
              "#\n"
        land = Land.Land(land_pirf, None, 60e0, 0.287e0, Globals.LAND_MODEL, doc)
    else:
        raise NameError("Unknown land model.")
    return land


if __name__ == "__main__":
    dt = Globals.DELTA_T
    t2x = float(input("Climate sensitivity (K): "))
    t_dep = True if input("T dependence (True/False): ") == "True" else False
    co2_dep = True if input("Co2 dependence (True/False): ") == "True" else False
    scenario = input("scenario (string as in file name forcing_<scenario>.dat): ")
    ID = input("additional simulation identifier (string, press enter to skip): ")
    print()
    if ID.strip() <= "":
        ID = os.getenv("ID")
    utils = UtilsFile.UtilsFile(ID)
    forcing = utils.readForcing(scenario)
    tot = 0
    for i in range(100):
        start = time.time()
        earth = Earth.Earth(createOceanModel(), createLandModel(), Atmosphere.Atmosphere(),
                            forcing, t2x, t_dep, co2_dep, scenario, dt)
        earth.carbonCycleClimateSimulation()
        tot += time.time() - start
    print("Time = ", str(tot / 100))
    utils.output(earth)
