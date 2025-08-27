"""
Field computation module for AMEDA (Angular Momentum Eddy Detection Algorithm).

This module handles:
- Loading velocity fields from AVISO NetCDF data
- Computing detection fields (KE, divergence, vorticity, OW, LOW, LNAM)
- Grid interpolation and mask management
"""

import numpy as np
import xarray as xr
from scipy.interpolate import interp2d, griddata
from typing import Dict, Tuple, Optional
from dataclasses import dataclass
import numba as nb

from ameda.params import AMEDAParams, AMEDADerivedParams
from ameda.utils import get_missing_val_2d, EARTH_RADIUS


@dataclass
class DetectionFields:
    """Container for computed detection fields"""
    step: int
    ke: np.ndarray      # Kinetic energy
    div: np.ndarray     # Divergence
    vort: np.ndarray    # Vorticity
    OW: np.ndarray      # Okubo-Weiss
    LOW: np.ndarray     # Local Okubo-Weiss
    LNAM: np.ndarray    # Local Normalized Angular Momentum


def load_fields_aviso(
    dataset: xr.Dataset,
    step: int,
    params: AMEDAParams,
    resol: int = 1,
    deg: int = 1
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Load velocity fields from AVISO NetCDF data.
    
    Parameters
    ----------
    dataset : xr.Dataset
        Input NetCDF dataset
    step : int
        Time step to load
    params : AMEDAParams
        Algorithm parameters
    resol : int
        Interpolation factor (1 = no interpolation)
    deg : int
        Degradation factor (1 = no degradation)
    
    Returns
    -------
    x, y : ndarray
        Grid coordinates (lon, lat)
    mask : ndarray
        Ocean mask (1 = ocean, 0 = land)
    u, v : ndarray
        Velocity components (m/s)
    ssh : ndarray or None
        Sea surface height (m) if type_detection >= 2
    """
    print(f"Loading fields at step {step}...")
    
    # Extract data at given time step
    data_step = dataset.isel(time=step)
    
    # Get coordinates and fields
    lon0 = data_step[params.x_name].values
    lat0 = data_step[params.y_name].values
    
    # Create 2D meshgrids if needed
    if lon0.ndim == 1:
        lon0, lat0 = np.meshgrid(lon0, lat0)
    
    # Get velocity fields
    u0 = data_step[params.u_name].values
    v0 = data_step[params.v_name].values
    
    # Create mask from velocity fields
    mask0 = ~(np.isnan(u0) | np.isnan(v0))
    mask0 = mask0.astype(float)
    
    # Get SSH if needed
    ssh0 = None
    if params.type_detection >= 2:
        if params.s_name in data_step:
            ssh0 = data_step[params.s_name].values
    
    # Apply degradation if requested
    if deg != 1:
        print(f"  Degrading fields by factor {deg}")
        x = lon0[::deg, ::deg]
        y = lat0[::deg, ::deg]
        mask = mask0[::deg, ::deg]
        u = u0[::deg, ::deg]
        v = v0[::deg, ::deg]
        if ssh0 is not None:
            ssh = ssh0[::deg, ::deg]
        else:
            ssh = None
    else:
        x = lon0
        y = lat0
        mask = mask0
        u = u0.copy()
        v = v0.copy()
        ssh = ssh0.copy() if ssh0 is not None else None
    
    N, M = x.shape
    
    # Interpolation
    if resol == 1:
        print("NO INTERPOLATION")
        
        # Set NaN in land
        u[mask == 0] = np.nan
        v[mask == 0] = np.nan
        
        # Enlarge mask by 1 pixel into land
        print("Enlarging coastal mask by 1 pixel...")
        for i in range(N):
            for j in range(M):
                if mask[i, j] == 0:
                    # Check if any neighbor is ocean
                    i_min, i_max = max(i-1, 0), min(i+1, N-1)
                    j_min, j_max = max(j-1, 0), min(j+1, M-1)
                    
                    if np.sum(mask[i_min:i_max+1, j_min:j_max+1]) > 0:
                        u[i, j] = 0
                        v[i, j] = 0
                        if ssh is not None and np.isnan(ssh[i, j]):
                            ssh_neighbors = ssh[i_min:i_max+1, j_min:j_max+1]
                            ssh[i, j] = np.nanmean(ssh_neighbors)
    
    else:
        print(f"Interpolating grid by factor {resol}")
        
        # New grid size
        Ni = resol * (N - 1) + 1
        Mi = resol * (M - 1) + 1
        
        # Grid spacing
        dy = (y[1, 0] - y[0, 0]) / resol if N > 1 else 0
        dx = (x[0, 1] - x[0, 0]) / resol if M > 1 else 0
        
        # Create interpolated grid
        xi, yi = np.meshgrid(
            np.arange(Mi) * dx + x.min(),
            np.arange(Ni) * dy + y.min()
        )
        
        # Interpolate mask
        mask_interp = interp2d(x[0, :], y[:, 0], mask, kind='linear', 
                              bounds_error=False, fill_value=0)
        maski = mask_interp(xi[0, :], yi[:, 0])
        maski = np.where(maski >= 0.5, 1.0, 0.0)
        
        # Enlarge mask by 1 pixel
        maski1 = maski.copy()
        print("Enlarging coastal mask by 1 pixel...")
        for i in range(Ni):
            for j in range(Mi):
                if maski[i, j] == 0:
                    i_min, i_max = max(i-1, 0), min(i+1, Ni-1)
                    j_min, j_max = max(j-1, 0), min(j+1, Mi-1)
                    
                    if np.sum(maski[i_min:i_max+1, j_min:j_max+1]) > 0:
                        maski1[i, j] = 1
        
        # Set land velocities to 0 for interpolation
        u[mask == 0 | np.isnan(u)] = 0
        v[mask == 0 | np.isnan(v)] = 0
        
        # Interpolate velocity fields
        u_interp = interp2d(x[0, :], y[:, 0], u, kind='cubic',
                          bounds_error=False, fill_value=0)
        v_interp = interp2d(x[0, :], y[:, 0], v, kind='cubic',
                          bounds_error=False, fill_value=0)
        
        ui = u_interp(xi[0, :], yi[:, 0])
        vi = v_interp(xi[0, :], yi[:, 0])
        
        # Interpolate SSH if needed
        if ssh is not None:
            ssh1 = get_missing_val_2d(x, y, ssh)
            ssh_interp = interp2d(x[0, :], y[:, 0], ssh1, kind='cubic',
                                bounds_error=False, fill_value=np.nan)
            sshi = ssh_interp(xi[0, :], yi[:, 0])
        else:
            sshi = None
        
        # Apply enlarged mask
        ui[maski1 == 0] = np.nan
        vi[maski1 == 0] = np.nan
        if sshi is not None:
            sshi[maski1 == 0] = np.nan
        
        # Export interpolated fields
        x = xi
        y = yi
        mask = maski
        u = ui
        v = vi
        ssh = sshi
    
    return x, y, mask, u, v, ssh


@nb.jit(nopython=True, parallel=False)
def compute_lnam_low_numba(
    uu: np.ndarray,
    vv: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    okubo: np.ndarray,
    b: np.ndarray,
    f: np.ndarray,
    grid_ll: bool,
    R: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute LNAM and LOW fields using Numba for acceleration.
    
    Parameters
    ----------
    uu, vv : ndarray
        Velocity components
    x, y : ndarray
        Grid coordinates
    okubo : ndarray
        Okubo-Weiss field
    b : ndarray
        Half box size for LNAM calculation
    f : ndarray
        Coriolis parameter
    grid_ll : bool
        True if coordinates are lon/lat
    R : float
        Earth radius factor for lon/lat grids
    
    Returns
    -------
    L : ndarray
        LNAM field
    LOW : ndarray
        Local Okubo-Weiss field
    """
    N, M = uu.shape
    L = np.zeros_like(uu)
    LOW = np.full_like(uu, np.nan)
    
    borders = int(np.max(b)) + 1
    
    for i in range(borders, N - borders + 1):
        for j in range(borders, M - borders + 1):
            
            if not np.isnan(vv[i, j]):
                bi = int(b[i, j])
                
                # Calculate LOW
                OW = okubo[i-bi:i+bi+1, j-bi:j+bi+1]
                LOW[i, j] = np.nanmean(OW)
                
                # Calculate LNAM
                xlocal = x[i-bi:i+bi+1, j-bi:j+bi+1]
                ylocal = y[i-bi:i+bi+1, j-bi:j+bi+1]
                ulocal = uu[i-bi:i+bi+1, j-bi:j+bi+1]
                vlocal = vv[i-bi:i+bi+1, j-bi:j+bi+1]
                
                # Center coordinates
                coordcentre = bi
                xc = xlocal[coordcentre, coordcentre]
                yc = ylocal[coordcentre, coordcentre]
                
                if grid_ll:
                    # Convert to km distances
                    d_xcentre = (xlocal - xc) * R * np.cos(ylocal * np.pi / 180)
                    d_ycentre = (ylocal - yc) * R
                else:
                    d_xcentre = xlocal - xc
                    d_ycentre = ylocal - yc
                
                # Angular momentum calculation
                cross = d_xcentre * vlocal - d_ycentre * ulocal
                dot = ulocal * d_xcentre + vlocal * d_ycentre
                produit = np.sqrt(ulocal**2 + vlocal**2) * np.sqrt(d_xcentre**2 + d_ycentre**2)
                
                sum_cross = np.nansum(cross)
                sum_dp = np.nansum(dot) + np.nansum(produit)
                
                if sum_dp != 0:
                    L[i, j] = sum_cross / sum_dp * np.sign(f[i, j])
                else:
                    L[i, j] = 0
    
    return L, LOW


def compute_fields(
    x: np.ndarray,
    y: np.ndarray,
    mask: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    params: AMEDAParams,
    derived_params: AMEDADerivedParams,
    step: int,
    resol: int = 1
) -> DetectionFields:
    """
    Compute detection fields from velocity data.
    
    Parameters
    ----------
    x, y : ndarray
        Grid coordinates
    mask : ndarray
        Ocean mask
    u, v : ndarray
        Velocity components (m/s)
    params : AMEDAParams
        Algorithm parameters
    derived_params : AMEDADerivedParams
        Derived parameters (b, f, etc.)
    step : int
        Time step
    resol : int
        Interpolation factor
    
    Returns
    -------
    fields : DetectionFields
        Computed detection fields
    """
    print(f"Computing fields for step {step}")
    
    # Use interpolated parameters if resol > 1
    if resol > 1:
        b = derived_params.bi
        f = derived_params.fi
    else:
        b = derived_params.b
        f = derived_params.f
    
    # Ensure integer b values and handle NaN
    b = np.nan_to_num(b, nan=1.0)
    b = np.round(b).astype(int)
    b[b < 1] = 1
    
    # Calculate kinetic energy
    ke = (u**2 + v**2) / 2
    
    # Initialize derivative arrays
    dx_arr = np.zeros_like(x)
    dy_arr = np.zeros_like(y)
    dux = np.zeros_like(u)
    duy = np.zeros_like(u)
    dvx = np.zeros_like(v)
    dvy = np.zeros_like(v)
    
    # Compute spatial derivatives
    dx_arr[1:-1, 1:-1] = x[1:-1, 2:] - x[1:-1, :-2]
    dy_arr[1:-1, 1:-1] = y[2:, 1:-1] - y[:-2, 1:-1]
    
    if params.grid_ll:
        # Convert to meters
        R = EARTH_RADIUS * np.pi / 180  # km per degree
        dx_arr = dx_arr * R * np.cos(np.deg2rad(y))
        dy_arr = dy_arr * R
    
    # Convert to meters
    dx_arr = dx_arr * 1000
    dy_arr = dy_arr * 1000
    
    # Compute velocity derivatives
    dux[1:-1, 1:-1] = u[1:-1, 2:] - u[1:-1, :-2]
    duy[1:-1, 1:-1] = u[2:, 1:-1] - u[:-2, 1:-1]
    dvx[1:-1, 1:-1] = v[1:-1, 2:] - v[1:-1, :-2]
    dvy[1:-1, 1:-1] = v[2:, 1:-1] - v[:-2, 1:-1]
    
    # Avoid division by zero
    dx_arr[dx_arr == 0] = np.nan
    dy_arr[dy_arr == 0] = np.nan
    
    # Calculate Okubo-Weiss components
    sn = dux / dx_arr - dvy / dy_arr  # Shear
    ss = dvx / dx_arr + duy / dy_arr  # Strain
    om = dvx / dx_arr - duy / dy_arr  # Vorticity
    
    okubo = sn**2 + ss**2 - om**2
    
    # Calculate divergence
    div = dux / dx_arr + dvy / dy_arr
    
    # Calculate vorticity (sign-adjusted)
    vorticity = om * np.sign(f)
    
    # Calculate LNAM and LOW
    print("Computing LNAM...")
    
    # Ensure consistent float64 types for Numba
    u = u.astype(np.float64)
    v = v.astype(np.float64)
    x = x.astype(np.float64)
    y = y.astype(np.float64)
    okubo = okubo.astype(np.float64)
    f = f.astype(np.float64)
    
    # Use Numba-accelerated function
    R_km = EARTH_RADIUS * np.pi / 180 if params.grid_ll else 1.0
    L, LOW = compute_lnam_low_numba(
        u, v, x, y, okubo, b, f,
        params.grid_ll, R_km
    )
    
    # Handle NaN values
    L[np.isnan(L)] = 0
    
    # Apply mask to all fields
    fields = DetectionFields(
        step=step,
        ke=ke * mask,
        div=div * mask,
        vort=vorticity * mask,
        OW=okubo * mask,
        LOW=LOW * mask,
        LNAM=L * mask
    )
    
    print("")
    return fields


# Main processing function
def process_fields(
    nc_file: str,
    step: int,
    params: AMEDAParams,
    derived_params: AMEDADerivedParams
) -> DetectionFields:
    """
    Main function to load data and compute detection fields.
    
    Parameters
    ----------
    nc_file : str
        Path to NetCDF file
    step : int
        Time step to process
    params : AMEDAParams
        Algorithm parameters
    derived_params : AMEDADerivedParams
        Derived parameters
    
    Returns
    -------
    fields : DetectionFields
        Computed detection fields
    """
    # Load dataset
    ds = xr.open_dataset(nc_file)
    
    # Load fields
    x, y, mask, u, v, ssh = load_fields_aviso(
        ds, step, params, params.resol, params.deg
    )
    
    # Compute detection fields
    fields = compute_fields(
        x, y, mask, u, v,
        params, derived_params,
        step, params.resol
    )
    
    ds.close()
    
    return fields


if __name__ == "__main__":
    # Test with sample data
    from ameda.params import AMEDAParams, compute_derived_params
    
    DATASET = "/Users/gianlucacalo/Desktop/projects/ocean_surface_levels/data/raw/cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.0625deg_P1D.nc"
    
    # Load test data
    ds = xr.open_dataset(DATASET)
    step = 0
    
    # Get parameters with correct field names
    params = AMEDAParams(
        x_name="longitude",
        y_name="latitude",
        u_name="ugos",
        v_name="vgos",
        s_name="sla"
    )
    
    # Get first time step for parameter computation
    data_step = ds.isel(time=0)
    lon = data_step.longitude.values
    lat = data_step.latitude.values
    ugos = data_step.ugos.values
    vgos = data_step.vgos.values
    mask = ~np.isnan(ugos)
    
    # Compute derived parameters
    derived_params = compute_derived_params(lon, lat, mask, ugos, vgos, params)
    
    # Process fields
    fields = process_fields(DATASET, step, params, derived_params)
    
    print(f"Fields computed successfully:")
    print(f"  KE range: {np.nanmin(fields.ke):.6f} to {np.nanmax(fields.ke):.6f}")
    print(f"  Vorticity range: {np.nanmin(fields.vort):.6e} to {np.nanmax(fields.vort):.6e}")
    print(f"  LNAM range: {np.nanmin(fields.LNAM):.6f} to {np.nanmax(fields.LNAM):.6f}")
    
    ds.close()