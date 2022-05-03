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
inline constexpr std::string_view OCEAN_MODEL = "HILDA";
inline constexpr std::string_view LAND_MODEL = "HRBM";
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

#endif //CPP_GLOBALS_HPP
