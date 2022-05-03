/**
 * Brice Petit 000439675
 * Master Thesis
 * Transformation of BernSCM using UML and OOP
 * CPP version
 */

#ifndef CPP_GLOBALS_HPP
#define CPP_GLOBALS_HPP

#include <string>

// In this file, we will find all global variable used for the configuration of our model,
// variable parametrized and other useful global variables.

// ---------------------------Configuration---------------------------- //

inline constexpr double DELTA_T = 10e0;
inline constexpr int LINEAR = 1;
inline constexpr int IMPLICIT_L = 1;
inline constexpr int IMPLICIT_O = 1;
inline constexpr int TEQUIL = 0;
inline constexpr int LAND_T_DEP = 1;
inline constexpr int DO_INTERPOL = 1;
inline constexpr int WEIGHT_SIZE = 7;
inline constexpr int T_SCALE_SIZE = 7;
// ----------------------------Input/Output---------------------------- //
// Number of forces
inline constexpr int NFORC = 7;
// Column index of each forcing
// Timed (yr)
inline constexpr int JTIME = 0;
// Temperature (K)
inline constexpr int JTEMP = 1;
// Non-CO₂ RF (Wm⁻²)
inline constexpr int JRFNC = 2;
// Budget RF (Wm⁻²)
inline constexpr int JRFB = 3;
// Atmospheric CO₂ (ppm)
inline constexpr int JACO2 = 4;
// Anthropogenic C emissions (Gt/yr)
inline constexpr int JECO2 = 5;
// Budget C uptake (Gt/yr)
inline constexpr int JFB = 6;

// --------------------------General Constants------------------------- //

// C mass equivalent for Atmospheric CO₂ concentration (ppm/GtC)
inline constexpr double PPMTOGT = 2.123e0;
// Micromole mass of carbon (g/μmol)
inline constexpr double MUMOL = 12.0107e-6;
// Ocean fraction of earth surface
inline constexpr double OFRAC = 0.71e0;
// Factor 10E15 for scaling g/Pg, W/PW
inline constexpr double PETA = 1e15;
// Rad. efficiency parameter for CO₂
inline constexpr double RECO2 = 5.35e0;
// Log(2)
inline constexpr double LOG2 = 0.6931471805599453;
// Seconds per year (s/yr, joos PR)
inline constexpr double SECTOYEAR = 365e0 * 24e0 * 3600e0;
// Reference preindustrial CO₂ for RF calculation (ppm)
inline constexpr double CO2PREIND = 278.6822074e0;
// RF for CO₂ doubling
inline constexpr double RF2X = RECO2 * LOG2;
// Undefined
inline constexpr double NA = -9999.9999e0;
// Maximal number of finite timescales of boxes
inline constexpr int NSCALEMAX = 10;
// The number of row in the forcing matrix
inline int FORCING_ROWS = 0;
// -------------------------------Struct------------------------------- //
struct Pirf {
    std::string name;
    int nscale = 0;
    double *weight{};
    double *t_scale{};
};

struct Prop {
    int nscale = 0;
    double x = 0e0;
    double propm[NSCALEMAX + 1] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    double propf[NSCALEMAX + 1] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
    double propfo[NSCALEMAX + 1] = {0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0};
};

struct Sirf {
    double t_max = 0e0;
    double *weight{};
    double *t_scale{};
};
// --------------------------Global Variables-------------------------- //
// We use global variables to do like inc file in fortran

// ####################
// # Input Variables. #
// ####################
// Climate sensitivity
inline double t2x = 0.0;
// Temperature dependence
inline bool t_dep = false;
// Co2 dependence
inline bool co2_dep = false;
// The scenario of the simulation
inline std::string scenario = "";
// Additional simulation identifier.
inline std::string ID = "";

// ######################
// # General Variables. #
// ######################

// Budget C uptake boolean
inline bool f_budget = false;
// Budget radiative forcing boolean
inline bool rf_budget = false;
// CO2 budget boolean
inline bool co2_budget = false;
// Preindustrial equilibrium CO2 concentration
inline double co2_atm0 = 0e0;

// forcing is a matrix where the column represents:
// [Year, glob_temp_dev, RF_nonCO2, RF_budget, co2_atm, fossil_CO2_em, budget_C_uptake]
inline double **forcing{};

// Number of time steps
inline int ntime = 0;

// Time model
inline double *_time{};
// Forcing time series
inline long int *itime{};
// Global mean SAT deviation from preindustrial (K)
inline double *temp{};
// Global mean SAT deviation from preindustrial (K)
inline double *tempk{};
// Non-CO₂ radiative forcing (Wm⁻²)
inline double *rfnc{};
// CO₂ radiative forcing (Wm⁻²)
inline double *rfc{};
// Budget radiative forcing (Wm⁻²)
inline double *rfb{};
// Total radiative forcing (Wm⁻²)
inline double *rf{};
// Air-sea heat flux (PW)
inline double *fh{};
// Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm)
inline double *dpcs{};
// Mass of C in atmosphere (GtC)
inline double *ma{};
// Mass of C in ocean surface layer (GtC)
inline double *ms{};
// Ocean C surface pools (GtC)
inline double *msk{};
// Mass of C in land biosphere (GtC)
inline double *ml{};
// Land biosphere C pools (GtC)
inline double *mlk{};
// CO₂ emissions (GtC/yr)
inline double *e_co2{};
// Budget C uptake (GtC/yr)
inline double *fb{};
// NPP (GtC/yr)
inline double *fnpp{};
// Air-sea C flux (GtC/yr)
inline double *fo{};

// ###################
// # Land variables. #
// ###################
inline constexpr std::string_view LAND_MODEL = "HRBM";
inline struct Pirf land_pirf;
inline struct Sirf land_sirf;
inline std::string LAND_DOC = "";
// Propagator for the land
inline struct Prop land_prop;
inline double NPP0 = 0e0;
inline double FERT = 0e0;

// ####################
// # Ocean variables. #
// ####################
inline constexpr std::string_view OCEAN_MODEL = "HILDA";
inline struct Pirf ocean_pirf;
inline std::string OCEAN_DOC = "";
// Propagator for the ocean
inline struct Prop ocean_prop;

// Mixed layer depth (m)
inline double hmix = 0e0;
// Heat capacity (J/kg/K)
inline double cp = 0e0;
// Water density (kg/m³)
inline double dens = 0e0;
// Water density (value used for DIC)
inline double dens_c = 0e0;
// Sea Surface Area (m²)
inline double aoc = 0e0;
// Gas exchange coefficient (1/yr).
inline double kg_aoc = 0e0;
// SST for calculating pCO2~DIC dependence (K); fixed, not=Temp
inline double t_chem = 0e0;
// om_t: Multiplier for heat uptake (K/(PW*yr)).
inline double om_t = 0e0;
// om_c: GtC→DIC conversion factor (umol/kg/Gt).
inline double om_c = 0e0;
// Ocean chemistry T-dependence (Takahashi 1993)
inline constexpr double buffer_t = 0.0423e0;

// #########################
// # Exceed verifications. #
// #########################
inline bool exceed = false;
inline bool co2_exceed = false;
inline bool t_exceed = false;

#endif //CPP_GLOBALS_HPP
