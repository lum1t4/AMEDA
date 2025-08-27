import numpy as np
from pathlib import Path
from scipy.interpolate import griddata
from matplotlib.path import Path as MPath
from typing import Tuple, List, Dict, Optional


ROOT = Path(__file__).resolve().parent

# Constants
DEG2RAD = np.pi / 180
DEG2NM = 60.0
NM2KM = 1.8520
EARTH_RADIUS = 6378.137  # km


def is_regular_grid(lon, lat):
    """Check if grid is regular"""
    dy_var = np.var(np.diff(lat, axis=0))
    dx_var = np.var(np.diff(lon, axis=1))
    return dy_var < 1e-6 and dx_var < 1e-6


def sw_dist2(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """
    Calculate distance between consecutive lat,lon coordinates in km.
    Uses the 'Plane Sailing' approximation.
    
    Parameters
    ----------
    lat : array-like
        Latitudes in decimal degrees
    lon : array-like
        Longitudes in decimal degrees
        
    Returns
    -------
    dist : ndarray
        Distances between consecutive positions in km
    """
    lat = np.asarray(lat).ravel()
    lon = np.asarray(lon).ravel()
    
    dlon = np.diff(lon)
    # Wrap longitude difference to [-180, 180]
    dlon = np.where(np.abs(dlon) > 180,
                    -np.sign(dlon) * (360 - np.abs(dlon)),
                    dlon)
    
    latrad = np.abs(lat * DEG2RAD)
    dep = np.cos((latrad[1:] + latrad[:-1]) / 2.0) * dlon
    dlat = np.diff(lat)
    
    # Distance in nautical miles then convert to km
    dist = DEG2NM * np.sqrt(dlat**2 + dep**2) * NM2KM
    
    return dist


def mean_radius(xy: np.ndarray, grid_ll: bool = True) -> Tuple[np.ndarray, float, float, np.ndarray]:
    """
    Compute the mean radius (R), perimeter (P), surface area (A),
    and barycenter (ll) of a closed contour.
    
    Parameters
    ----------
    xy : ndarray
        2xN array with x (or lon) in first row, y (or lat) in second row
    grid_ll : bool
        True if coordinates are (lon, lat), False if in km
        
    Returns
    -------
    R : ndarray
        Array of 4 different radius computations (km)
    A : float
        Surface area (km²)
    P : float
        Perimeter (km)
    ll : ndarray
        Barycenter coordinates [x, y]
    """
    lim = xy.shape[1] - 1
    
    # Barycenter computation
    ll = np.array([np.sum(xy[0, :lim]) / lim,
                   np.sum(xy[1, :lim]) / lim])
    
    # Distance from barycenter to each point
    distance = np.zeros(lim + 1)
    for point in range(lim + 1):
        xs = np.array([ll[0], xy[0, point]])
        ys = np.array([ll[1], xy[1, point]])
        
        if grid_ll:
            distance[point] = sw_dist2(ys, xs)[0] if len(xs) > 1 else 0
        else:
            distance[point] = np.sqrt((xs[1] - xs[0])**2 + (ys[1] - ys[0])**2)
    
    # Distance between consecutive points
    if grid_ll:
        distance2 = sw_dist2(xy[1, :], xy[0, :])
    else:
        distance2 = np.sqrt(np.diff(xy[0, :])**2 + np.diff(xy[1, :])**2)
    
    # Perimeter
    P = np.sum(distance2)
    
    # Surface area using triangulation
    aire = np.zeros(lim)
    for point in range(lim):
        a = distance[point]
        b = distance[point + 1]
        c = distance2[point]
        s = (a + b + c) / 2
        # Heron's formula
        aire[point] = np.sqrt(np.maximum(0, s * (s - a) * (s - b) * (s - c)))
    
    A = np.sum(np.real(aire))
    
    # Compute different radius definitions
    R = np.zeros(4)
    
    # Radius of equivalent circle
    R[0] = np.sqrt(A / np.pi) if A > 0 else 0
    
    # Radius weighted by surface area
    param = np.zeros(lim)
    for point in range(lim):
        rayonmoy = (distance[point] + distance[point + 1]) / 2
        param[point] = aire[point] * rayonmoy
    R[1] = np.sum(param) / A if A > 0 else 0
    
    # Max radius from barycenter
    R[2] = np.max(distance[:lim])
    
    # Mean radius from barycenter
    R[3] = np.mean(distance[:lim])
    
    return R, A, P, ll


def in_polygon(x: np.ndarray, y: np.ndarray, 
                xv: np.ndarray, yv: np.ndarray) -> np.ndarray:
    """
    Test if points (x, y) are inside polygon defined by vertices (xv, yv).
    
    Parameters
    ----------
    x, y : ndarray
        Coordinates to test
    xv, yv : ndarray
        Polygon vertices
        
    Returns
    -------
    inside : ndarray (bool)
        True for points inside or on the polygon
    """
    x = np.asarray(x)
    y = np.asarray(y)
    xv = np.asarray(xv)
    yv = np.asarray(yv)
    
    original_shape = x.shape
    x_flat = x.ravel()
    y_flat = y.ravel()
    
    # Create polygon path
    vertices = np.column_stack([xv, yv])
    path = MPath(vertices)
    
    # Test points
    points = np.column_stack([x_flat, y_flat])
    inside = path.contains_points(points)
    
    return inside.reshape(original_shape)


def scan_lines(CS) -> Tuple[List[Dict], np.ndarray]:
    """
    Rearrange contour data into structured format.
    
    Parameters
    ----------
    CS : contour object or array
        Contour data from matplotlib.contour or similar format
        
    Returns
    -------
    lines : list of dicts
        Each dict contains 'x', 'y', and 'l' (max y) for a contour
    lvl : ndarray
        Contour level values
    """
    lines = []
    lvl = []
    
    if hasattr(CS, 'allsegs'):  # matplotlib ContourSet
        for level_idx, level in enumerate(CS.levels):
            for seg in CS.allsegs[level_idx]:
                if len(seg) > 0:
                    lines.append({
                        'x': seg[:, 0],
                        'y': seg[:, 1],
                        'l': np.max(seg[:, 1])
                    })
                    lvl.append(level)
    else:  # Raw contour array format (like MATLAB contourc output)
        k = 0
        while k < CS.shape[1]:
            npoints = int(CS[1, k])
            level = CS[0, k]
            x_vals = CS[0, k+1:k+1+npoints]
            y_vals = CS[1, k+1:k+1+npoints]
            
            lines.append({
                'x': x_vals,
                'y': y_vals,
                'l': np.max(y_vals)
            })
            lvl.append(level)
            
            k = k + npoints + 1
    
    # Sort by maximum y coordinate
    if lines:
        sorted_indices = np.argsort([line['l'] for line in lines])
        lines = [lines[i] for i in sorted_indices]
        lvl = np.array(lvl)[sorted_indices]
    else:
        lvl = np.array([])
    
    return lines, lvl


def get_missing_val_2d(x: np.ndarray, y: np.ndarray, 
                        field: np.ndarray, missvalue: float = np.nan,
                        default: float = 0) -> np.ndarray:
    """
    Fill missing values in a 2D field using nearest neighbor interpolation.
    
    Parameters
    ----------
    x, y : ndarray
        Grid coordinates
    field : ndarray
        2D field with missing values
    missvalue : float
        Value indicating missing data
    default : float
        Default value if no valid data
        
    Returns
    -------
    field : ndarray
        Field with missing values filled
    """
    field = field.copy()
    
    if np.isnan(missvalue):
        ismask = np.isnan(field)
    else:
        ismask = (field == missvalue)
    
    isdata = ~ismask
    
    # Check if there's any valid data
    if not np.any(isdata):
        return np.full_like(field, default)
    
    if np.sum(isdata) < 6:
        default = np.min(field[isdata])
        return np.full_like(field, default)
    
    if not np.any(ismask):
        return field
    
    # Fill missing values using nearest neighbor
    if x.ndim == 1:
        x, y = np.meshgrid(x, y)
    
    valid_points = np.column_stack([x[isdata].ravel(), y[isdata].ravel()])
    valid_values = field[isdata].ravel()
    invalid_points = np.column_stack([x[ismask].ravel(), y[ismask].ravel()])
    
    filled_values = griddata(valid_points, valid_values, invalid_points,
                            method='nearest', fill_value=default)
    
    field[ismask] = filled_values
    
    return field


def cumtrapz(y: np.ndarray, x: np.ndarray = None, axis: int = -1, initial: float = 0) -> np.ndarray:
    """
    Cumulative trapezoidal integration, similar to MATLAB's cumtrapz.
    
    Parameters
    ----------
    y : ndarray
        Values to integrate
    x : ndarray, optional
        Sample points (if None, unit spacing assumed)
    axis : int
        Axis along which to integrate
    initial : float
        Initial value
        
    Returns
    -------
    result : ndarray
        Cumulative integral
    """
    from scipy.integrate import cumulative_trapezoid
    
    if x is not None:
        # Handle different input shapes
        if x.ndim == y.ndim:
            # x has same shape as y - extract along axis
            dx = np.diff(x, axis=axis)
        else:
            dx = np.diff(x)
    else:
        dx = 1
    
    result = cumulative_trapezoid(y, x, axis=axis, initial=initial)
    
    return result
