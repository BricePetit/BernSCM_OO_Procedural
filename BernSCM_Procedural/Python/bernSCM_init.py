"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM in procedural
Python version
"""
import numpy as np

import Globals
import Prop
import bernSCM


def initialize():
    """
    Initialization for the simulation.
    """
    nin = Globals.forcing.shape[0] - 1
    if Globals.LINEAR and Globals.TEQUIL > 0:
        raise NameError("Equilibrated time steps not implemented for linear.")
    elif Globals.LINEAR and not ((Globals.IMPLICIT_O or Globals.IMPLICIT_L) == 1):
        raise NameError("Linear scheme only implicit.")
    tin = Globals.forcing[:, Globals.JTIME]
    Globals.ntime = int((tin[nin] - tin[0]) / Globals.DELTA_T + 1e-6) + 1
    j = 0
    Globals.time = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.itime = np.zeros(Globals.ntime, dtype="int64")
    for i in range(Globals.ntime):
        Globals.time[i] = tin[0] + i * Globals.DELTA_T
        while Globals.time[i] > tin[j + 1]:
            j += 1
        Globals.itime[i] = j

    Globals.rf = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.rfnc = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.rfc = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.rfb = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.fh = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.fb = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.fo = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.ms = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.fnpp = np.full(Globals.ntime, Globals.NPP0, dtype=np.longdouble)
    Globals.ml = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.ma = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.e_co2 = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.dpcs = np.zeros(Globals.ntime, dtype=np.longdouble)

    Globals.land_prop = Prop.Prop(0, 0e0, np.zeros(Globals.NSCALEMAX + 1, dtype=np.longdouble),
                                  np.zeros(Globals.NSCALEMAX + 1, dtype=np.longdouble), 
                                  np.zeros(Globals.NSCALEMAX + 1, dtype=np.longdouble))
    Globals.ocean_prop = Prop.Prop(0, 0e0, np.zeros(Globals.NSCALEMAX + 1, dtype=np.longdouble),
                                   np.zeros(Globals.NSCALEMAX + 1, dtype=np.longdouble), 
                                   np.zeros(Globals.NSCALEMAX + 1, dtype=np.longdouble))

    bernSCM.computePropagators(Globals.land_pirf, Globals.land_prop)
    bernSCM.computePropagators(Globals.ocean_pirf, Globals.ocean_prop)

    Globals.mlk = np.zeros(Globals.land_pirf.nscale + 1, dtype=np.longdouble)

    for i in range(Globals.land_pirf.nscale):
        Globals.mlk[i] = (Globals.NPP0 * Globals.land_pirf.weight[i] * Globals.land_pirf.t_scale[i])
        Globals.ml[0] += Globals.mlk[i]

    Globals.msk = np.zeros(Globals.ocean_pirf.nscale + 1, dtype=np.longdouble)
    Globals.ms = np.zeros(Globals.ntime, dtype=np.longdouble)
    Globals.tempk = np.zeros(Globals.ocean_pirf.nscale + 1, dtype=np.longdouble)
    Globals.temp = np.zeros(Globals.ntime, dtype=np.longdouble)
