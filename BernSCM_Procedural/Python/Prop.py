"""
Brice Petit 000439675
Master Thesis
Transformation of BernSCM in procedural
Python version
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class Prop:
    """
    This class represents a c-like structure of propagators.
    """
    nscale: int
    x: float
    propm: np.ndarray
    propf: np.ndarray
    propfo: np.ndarray
