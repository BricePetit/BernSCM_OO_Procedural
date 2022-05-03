"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM in procedural
Python version
"""
import numpy as np

# In this file, we will find all global variable used for the configuration of our model,
# variable parametrized and other useful global variables.

# ---------------------------Configuration----------------------------#

DELTA_T = 10e0
LINEAR = 1
IMPLICIT_L = 1
IMPLICIT_O = 1
TEQUIL = 0
LAND_T_DEP = 1
DO_INTERPOL = 1

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

# --------------------------Global Variables--------------------------#
# We use global variables to do like inc file in fortran

####################
# Input Variables. #
####################
# Climate sensitivity
t2x = 0.0
# Temperature dependence
t_dep = False
# Co2 dependence
co2_dep = False
# The scenario of the simulation
scenario = ""
# Additional simulation identifier.
ID = ""

######################
# General Variables. #
######################

# Budget C uptake boolean
f_budget = False
# Budget radiative forcing boolean
rf_budget = False
# CO2 budget boolean
co2_budget = False
# Preindustrial equilibrium CO2 concentration
co2_atm0 = 0e0

# forcing is a matrix where the column represents:
# [Year, glob_temp_dev, RF_nonCO2, RF_budget, co2_atm, fossil_CO2_em, budget_C_uptake]
forcing = np.empty([0, 7], dtype=np.longdouble)

# Number of time steps
ntime = 0

# Time model
time = np.array([], dtype=np.longdouble)
# Forcing time series
itime = np.array([], dtype=np.longdouble)
# Global mean SAT deviation from preindustrial (K)
temp = np.array([], dtype=np.longdouble)
# Global mean SAT deviation from preindustrial (K)
tempk = np.array([], dtype=np.longdouble)
# Non-CO₂ radiative forcing (Wm⁻²)
rfnc = np.array([], dtype=np.longdouble)
# CO₂ radiative forcing (Wm⁻²)
rfc = np.array([], dtype=np.longdouble)
# Budget radiative forcing (Wm⁻²)
rfb = np.array([], dtype=np.longdouble)
# Total radiative forcing (Wm⁻²)
rf = np.array([], dtype=np.longdouble)
# Air-sea heat flux (PW)
fh = np.array([], dtype=np.longdouble)
# Ocean surface CO₂ saturation partial pressure deviation from preindustrial (ppm)
dpcs = np.array([], dtype=np.longdouble)
# Mass of C in atmosphere (GtC)
ma = np.array([], dtype=np.longdouble)
# Mass of C in ocean surface layer (GtC)
ms = np.array([], dtype=np.longdouble)
# Ocean C surface pools (GtC)
msk = np.array([], dtype=np.longdouble)
# Mass of C in land biosphere (GtC)
ml = np.array([], dtype=np.longdouble)
# Land biosphere C pools (GtC)
mlk = np.array([], dtype=np.longdouble)
# CO₂ emissions (GtC/yr)
e_co2 = np.array([], dtype=np.longdouble)
# Budget C uptake (GtC/yr)
fb = np.array([], dtype=np.longdouble)
# NPP (GtC/yr)
fnpp = np.array([], dtype=np.longdouble)
# Air-sea C flux (GtC/yr)
fo = np.array([], dtype=np.longdouble)

###################
# Land variables. #
###################
LAND_MODEL = "HRBM"
land_pirf = None
land_sirf = None
LAND_DOC = ""
# Propagator for the land
land_prop = None
NPP0 = 0e0
FERT = 0e0

####################
# Ocean variables. #
####################
OCEAN_MODEL = "HILDA"
ocean_pirf = None
OCEAN_DOC = ""
# Propagator for the ocean
ocean_prop = None

# Mixed layer depth (m)
hmix = 0e0
# Heat capacity (J/kg/K)
cp = 0e0
# Water density (kg/m³)
dens = 0e0
# Water density (value used for DIC)
dens_c = 0e0
# Sea Surface Area (m²)
aoc = 0e0
# Gas exchange coefficient (1/yr).
kg_aoc = 0e0
# SST for calculating pCO2~DIC dependence (K); fixed, not=Temp
t_chem = 0e0
# om_t: Multiplier for heat uptake (K/(PW*yr)).
om_t = 0e0
# om_c: GtC→DIC conversion factor (umol/kg/Gt).
om_c = 0e0
# Ocean chemistry T-dependence (Takahashi 1993)
buffer_t = 0.0423e0

#########################
# Exceed verifications. #
#########################
exceed = False
co2_exceed = False
t_exceed = False
