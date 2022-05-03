/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */
#include "bernSCM_init.hpp"
#include "bernSCM_readforcing.hpp"
#include "bernSCM_output.hpp"
#include "Globals.hpp"

#include <iostream>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

/**
 * Set global values for the ocean model.
 *
 * @param _hmix     Mixed layer depth (m).
 * @param _cp       Heat capacity (J/kg/K).
 * @param _dens     Water density (kg/m³).
 * @param _dens_c   Water density (value used for DIC).
 * @param _aoc      Sea Surface Area (m²).
 * @param _kg_aoc   Gas exchange coefficient (1/yr).
 * @param _t_chem   SST for calculating pCO2~DIC dependence (K); fixed, not=Temp.
 */
void setOceanValues(double _hmix, double _cp, double _dens, double _dens_c, double _aoc, double _kg_aoc,
                    double _t_chem) {
    hmix = _hmix;
    cp = _cp;
    dens = _dens;
    dens_c = _dens_c;
    aoc = _aoc;
    kg_aoc = _kg_aoc;
    t_chem = _t_chem;
    om_t = PETA * SECTOYEAR / (_hmix * _cp * _dens * _aoc);
    om_c = PETA / MUMOL / _dens_c / _hmix / _aoc;
}

/**
 * Function to create the Ocean representation.
 */
void createOceanModel() {
    if (OCEAN_MODEL == "HILDA") {
        double ocean_weight[WEIGHT_SIZE] = {0.27830433e0, 0.23337218e0, 0.13732822e0, 0.05154051e0, 0.03503318e0,
                                            0.24013944e0, .022936e0};
        double ocean_t_scale[T_SCALE_SIZE] = {0.45253504e0, 2.19901724e0, 12.03837102e0, 59.58359820e0, 237.30651757e0,
                                              0.03855458e0, NA};
        ocean_pirf.name = OCEAN_MODEL;
        ocean_pirf.nscale = 6;
        ocean_pirf.weight = new double[WEIGHT_SIZE];
        memcpy(ocean_pirf.weight, ocean_weight, sizeof(double) * WEIGHT_SIZE);
        ocean_pirf.t_scale = new double[T_SCALE_SIZE];
        memcpy(ocean_pirf.t_scale, ocean_t_scale, sizeof(double) * T_SCALE_SIZE);
        setOceanValues(75e0, 4000e0, 1028e0, 1026.5e0, 3.62e14, 1e0 / 9.06e0, 18.1716e0);
        OCEAN_DOC = "# HILDA ocean mixed layer pulse response\n" \
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
              "#\n";
    } else if (OCEAN_MODEL == "B2D") {
        double ocean_weight[WEIGHT_SIZE] = {0.09467125e0, 0.1029231e0, 0.03928349e0, 0.4593721e0, 0.0129862e0,
                                            0.2702249, 1.3691e-2};
        double ocean_t_scale[T_SCALE_SIZE] = {2.690038e0, 13.61728e0, 86.79685e0, 0.5762091e0, 337.2983e0, 0.07027151e0,
                                              NA};
        ocean_pirf.name = OCEAN_MODEL;
        ocean_pirf.nscale = 6;
        ocean_pirf.weight = new double[WEIGHT_SIZE];
        memcpy(ocean_pirf.weight, ocean_weight, sizeof(double) * WEIGHT_SIZE);
        ocean_pirf.t_scale = new double[T_SCALE_SIZE];
        memcpy(ocean_pirf.t_scale, ocean_t_scale, sizeof(double) * T_SCALE_SIZE);
        setOceanValues(50e0, 4000e0, 1028e0, 1026.5e0, 3.5375e14, 1e0 / 7.46e0, 18.2997e0);
        OCEAN_DOC = "# Bern2.5D ocean mixed layer pulse response\n" \
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
              "#\n";
    } else if (OCEAN_MODEL == "Princeton3D") {
        double ocean_weight[WEIGHT_SIZE] = {2.27446514e0, 0.06161763e0, 0.03726494e0, 1.28186186e0, 0.01956537e0,
                                            -2.70925536e0, 0.01481883e0};
        double ocean_t_scale[T_SCALE_SIZE] = {1.19761983e0, 16.67585709e0, 65.10188851e0, 2.00904478e0, 347.58378832e0,
                                              1.55213441e0, NA};
        ocean_pirf.name = OCEAN_MODEL;
        ocean_pirf.nscale = 6;
        ocean_pirf.weight = new double[WEIGHT_SIZE];
        memcpy(ocean_pirf.weight, ocean_weight, sizeof(double) * WEIGHT_SIZE);
        ocean_pirf.t_scale = new double[T_SCALE_SIZE];
        memcpy(ocean_pirf.t_scale, ocean_t_scale, sizeof(double) * T_SCALE_SIZE);
        setOceanValues(50.9e0, 4000e0, 1028e0, 1026.5e0, 3.55e14, 1e0 / 7.66e0, 17.7e0);
        OCEAN_DOC = "# Princeton AOGCM ocean effective mixed layer pulse response\n" \
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
              "#\n";
    } else {
        throw "Unknown ocean model.";
    }
}

/**
 * Function to create Land representation.
 */
void createLandModel() {
    if (LAND_MODEL == "HRBM") {
        double land_weight[WEIGHT_SIZE] = {-0.15431665e0, 0.56172938e0, 0.07486977e0, 0.41365999e0, 0.10405757e0, 0e0,
                                           NA};
        double land_t_scale[T_SCALE_SIZE] = {0.2010730e0, 1.4753908e0, 8.889799331e0, 74.098061812e0, 2.538054e02, NA,
                                             NA};
        land_pirf.name = LAND_MODEL;
        land_pirf.nscale = 5;
        land_pirf.weight = new double[WEIGHT_SIZE];
        memcpy(land_pirf.weight, land_weight, sizeof(double) * WEIGHT_SIZE);
        land_pirf.t_scale = new double[T_SCALE_SIZE];
        memcpy(land_pirf.t_scale, land_t_scale, sizeof(double) * T_SCALE_SIZE);
        double land_weight_sirf[WEIGHT_SIZE] = {0.14e0, 0.056e0, 0.072e0, 0.044e0, 0.069e0, 0e0, NA};
        double land_t_scale_sirf[T_SCALE_SIZE] = {0.056e0, 0.079e0, 0.057e0, 0.053e0, 0.036e0, NA, NA};
        land_sirf.t_max = 5e0;
        land_sirf.weight = new double[WEIGHT_SIZE];
        memcpy(land_sirf.weight, land_weight_sirf, sizeof(double) * WEIGHT_SIZE);
        land_sirf.t_scale = new double[T_SCALE_SIZE];
        memcpy(land_sirf.t_scale, land_t_scale_sirf, sizeof(double) * T_SCALE_SIZE);
        NPP0 = 41.73282e0;
        LAND_DOC = "# HRBM pulse response/multibox substitute model\n" \
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
              "#\n";
    } else if (LAND_MODEL == "4box") {
        double land_weight[WEIGHT_SIZE] = {2.006028671730e0, 0.26828e0, -1.56754896028e0, 0.29323e0, 0e0, NA, NA};
        double land_t_scale[T_SCALE_SIZE] = {2.857143e0, 20e0, 2.181818e0, 100e0, NA, NA, NA};
        land_pirf.name = LAND_MODEL;
        land_pirf.nscale = 4;
        land_pirf.weight = new double[WEIGHT_SIZE];
        memcpy(land_pirf.weight, land_weight, sizeof(double) * WEIGHT_SIZE);
        land_pirf.t_scale = new double[T_SCALE_SIZE];
        memcpy(land_pirf.t_scale, land_t_scale, sizeof(double) * T_SCALE_SIZE);
        NPP0 = 60e0;
        FERT = 0.287e0;
        LAND_DOC = "# 4-box Biosphere of Bern Model\n" \
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
              "#\n";
    } else {
        throw "Unknown land model.";
    }
}

int main() {
    double tot_time = 0;
    string id, input;
    cout << "Climate sensitivity (K): ";
    cin >> t2x;
    cout << "T dependence (True/False): ";
    cin >> input;
    t_dep = (input == "True");
    cout << "Co2 dependence (True/False): ";
    cin >> input;
    co2_dep = (input == "True");
    cout << "scenario (string as in file name forcing_<scenario>.dat): ";
    cin >> scenario;
    cout << "additional simulation identifier (string, press enter to skip): " << endl;
    getline(cin, id);
    readForcing();
    createOceanModel();
    createLandModel();
    for (int i = 0; i < 100; i++) {
        auto start = chrono::high_resolution_clock::now();
        initialize();
        setForcing(0);
        co2_atm0 = ma[0] / PPMTOGT;
        for (int n = 1; n < ntime; n++) {
            setForcing(n);
            timeStep(n);
        }
        auto end = chrono::high_resolution_clock::now();
        tot_time += double(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
    }
    cout << "Time = " << (tot_time / 100) * 1e-9 << setprecision(20) << endl;
    output();
    // Delete arrays
    delete[] _time;
    delete[] itime;
    delete[] rf;
    delete[] rfnc;
    delete[] rfc;
    delete[] rfb;
    delete[] fh;
    delete[] fb;
    delete[] fo;
    delete[] ma;
    delete[] e_co2;
    delete[] fnpp;
    delete[] mlk;
    delete[] ml;
    delete[] msk;
    delete[] ms;
    delete[] tempk;
    delete[] temp;
    // Delete the matrix of forcing values.
    for (int i = 0; i < FORCING_ROWS; i++) {
        delete[] forcing[i];
    }
    delete[] forcing;
    return 0;
}
