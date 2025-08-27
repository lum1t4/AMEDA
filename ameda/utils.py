import numpy as np
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def is_regular_grid(lon, lat):
    """Check if grid is regular"""
    dy_var = np.var(np.diff(lat, axis=0))
    dx_var = np.var(np.diff(lon, axis=1))
    return dy_var < 1e-6 and dx_var < 1e-6
