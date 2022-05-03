"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM in procedural
Python version
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class Pirf:
    """
    This class represents a c-like structure of the IRF parameter set.
    """
    name: str
    nscale: int
    weight: np.ndarray
    t_scale: np.ndarray
