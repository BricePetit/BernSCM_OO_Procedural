"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM using UML and OOP
Python version
"""
# In this file, we will find all global variable used for the configuration of our model,
# variable parametrized and other useful global variables.

# ---------------------------Configuration----------------------------#

DELTA_T = 10e0
LINEAR = 1
IMPLICIT_L = 1
IMPLICIT_O = 1
TEQUIL = 0
LAND_T_DEP = 1
OCEAN_MODEL = "HILDA"
LAND_MODEL = "HRBM"

# ----------------------------Input/Output----------------------------#
# Number of forces
NFORC = 7
# Column index of each forcing
# Timed (yr)
JTIME = 0
# Temperature (K)
JTEMP = 1
# Non-CO₂ RF (Wm⁻²)
JRFNC = 2
# Budget RF (Wm⁻²)
JRFB = 3
# Atmospheric CO₂ (ppm)
JACO2 = 4
# Anthropogenic C emissions (Gt/yr)
JECO2 = 5
# Budget C uptake (Gt/yr)
JFB = 6

# --------------------------General Constants-------------------------#

# C mass equivalent for Atmospheric CO₂ concentration (ppm/GtC)
PPMTOGT = 2.123e0
# Micromole mass of carbon (g/μmol)
MUMOL = 12.0107e-6
# Ocean fraction of earth surface
OFRAC = 0.71e0
# Factor 10E15 for scaling g/Pg, W/PW
PETA = 1e15
# Rad. efficiency parameter for CO₂
RECO2 = 5.35e0
# Log(2)
LOG2 = 0.6931471805599453
# Seconds per year (s/yr, joos PR)
SECTOYEAR = 365e0 * 24e0 * 3600e0
# Reference preindustrial CO₂ for RF calculation (ppm)
CO2PREIND = 278.6822074e0
# RF for CO₂ doubling
RF2X = RECO2 * LOG2
# Undefined
NA = -9999.9999e0
# Maximal number of finite timescales of boxes
NSCALEMAX = 10
