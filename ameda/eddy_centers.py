"""
Eddy center detection module for AMEDA (Angular Momentum Eddy Detection Algorithm).

This module handles:
- Detection of LNAM maxima where LOW < 0
- Validation of centers with streamfunction/SSH contours
- Psi field computation from velocities
- Contour scanning and analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from typing import Tuple, List, Dict, Optional
from dataclasses import dataclass, field
from scipy.interpolate import griddata

from ameda.params import AMEDAParams, AMEDADerivedParams
from ameda.fields import DetectionFields, load_fields_aviso
from ameda.utils import (
    in_polygon, mean_radius, scan_lines, sw_dist2,
    cumtrapz, EARTH_RADIUS
)


@dataclass
class EddyCenter:
    """Container for eddy center information"""
    step: int
    type: List[int] = field(default_factory=list)  # 1 = cyclonic, -1 = anticyclonic
    x: List[float] = field(default_factory=list)   # x coordinate (lon)
    y: List[float] = field(default_factory=list)   # y coordinate (lat)
    i: List[int] = field(default_factory=list)     # column index
    j: List[int] = field(default_factory=list)     # row index


def compute_psi(
    x: np.ndarray,
    y: np.ndarray,
    mask: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    ci: int,
    cj: int,
    grid_ll: bool = True
) -> np.ndarray:
    """
    Compute streamfunction (PSI) field by spatially integrating velocities.
    
    Parameters
    ----------
    x, y : ndarray
        Grid coordinates
    mask : ndarray
        Ocean mask
    u, v : ndarray
        Velocity components (already scaled by f/g if needed)
    ci, cj : int
        Center indices in the grid
    grid_ll : bool
        True if coordinates are (lon, lat)
    
    Returns
    -------
    psi : ndarray
        Streamfunction field
    """
    # Build mask
    mask = mask.copy()
    mask[mask == 0] = np.nan
    N0, M0 = mask.shape
    
    # Set NaN to zero for integration
    u = u.copy()
    v = v.copy()
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    
    # Prepare distance matrices for the 4 domains
    km_di = np.zeros((x.shape[0], x.shape[1] - 1))
    km_dj = np.zeros((y.shape[0] - 1, y.shape[1]))
    
    # Calculate distances along i (x) direction
    for i in range(x.shape[0]):
        if grid_ll:
            # Use great circle distance for lon/lat
            lons = x[i, :]
            lats = np.full_like(lons, y[i, 0])
            km_di[i, :] = sw_dist2(lats, lons)
        else:
            km_di[i, :] = np.sqrt(np.diff(x[i, :])**2 + np.diff(y[i, :])**2)
    
    # Calculate distances along j (y) direction
    for j in range(y.shape[1]):
        if grid_ll:
            lats = y[:, j]
            lons = np.full_like(lats, x[0, j])
            km_dj[:, j] = sw_dist2(lats, lons)
        else:
            km_dj[:, j] = np.sqrt(np.diff(x[:, j])**2 + np.diff(y[:, j])**2)
    
    # Length for the 4 domains (NE, SE, NW, SW)
    lx1 = u[cj:, :].shape[0]
    lx2 = u[:cj+1, :].shape[0]
    ly1 = u[:, ci:].shape[1]
    ly2 = u[:, :ci+1].shape[1]
    
    # Adjust km_di and km_dj for integration starting at ci, cj
    di = np.zeros((lx1 + lx2 - 1, ly1 + ly2 - 1))
    dj = np.zeros((lx1 + lx2 - 1, ly1 + ly2 - 1))
    
    if ci > 0:
        di[:, :ly2-1] = np.concatenate([
            km_di[:, :ci][:, ::-1],
        ], axis=1)
    if ci < km_di.shape[1]:
        di[:, ly2:] = km_di[:, ci:]
    
    if cj > 0:
        dj[:lx2-1, :] = np.concatenate([
            km_dj[:cj, :][::-1, :],
        ], axis=0)
    if cj < km_dj.shape[0]:
        dj[lx2:, :] = km_dj[cj:, :]
    
    # Integrate first row of v along x
    if ci < v.shape[1] - 1:
        v_seg = v[cj, ci:]
        cx1_raw = cumtrapz(v_seg, initial=0)
        di_seg = di[cj, ly2:] if cj < di.shape[0] else np.zeros(len(cx1_raw))
        cx1 = cx1_raw[:min(len(cx1_raw), len(di_seg))] * di_seg[:min(len(cx1_raw), len(di_seg))]
    else:
        cx1 = np.array([])
    
    if ci > 0:
        v_seg = v[cj, ci::-1]
        cx2_raw = cumtrapz(v_seg, initial=0)
        di_seg = di[cj, ly2-1::-1] if cj < di.shape[0] else np.zeros(len(cx2_raw))
        cx2 = -cx2_raw[:min(len(cx2_raw), len(di_seg))] * di_seg[:min(len(cx2_raw), len(di_seg))]
    else:
        cx2 = np.array([])
    
    # Integrate first column of u along y
    if cj < u.shape[0] - 1:
        u_seg = u[cj:, ci]
        cy1_raw = cumtrapz(u_seg, initial=0)
        dj_seg = dj[lx2:, ci] if ci < dj.shape[1] else np.zeros(len(cy1_raw))
        cy1 = -cy1_raw[:min(len(cy1_raw), len(dj_seg))] * dj_seg[:min(len(cy1_raw), len(dj_seg))]
    else:
        cy1 = np.array([])
    
    if cj > 0:
        u_seg = u[cj::-1, ci]
        cy2_raw = cumtrapz(u_seg, initial=0)
        dj_seg = dj[lx2-1::-1, ci] if ci < dj.shape[1] else np.zeros(len(cy2_raw))
        cy2 = cy2_raw[:min(len(cy2_raw), len(dj_seg))] * dj_seg[:min(len(cy2_raw), len(dj_seg))]
    else:
        cy2 = np.array([])
    
    # Initialize PSI arrays for 4 quadrants
    psi_xy = np.zeros((N0, M0))
    psi_yx = np.zeros((N0, M0))
    
    # Northeast quadrant
    if cj < N0 - 1 and ci < M0 - 1:
        ne_u = u[cj:, ci:]
        ne_v = v[cj:, ci:]
        ne_dj = dj[lx2:, ly2:]
        ne_di = di[lx2:, ly2:]
        
        if cx1.size > 0:
            mcx11 = np.repeat(cx1[np.newaxis, :], ne_u.shape[0], axis=0)
            psi_xy11 = mcx11 - cumtrapz(ne_u, axis=0, initial=0) * ne_dj
            psi_xy[cj:, ci:] = psi_xy11[:N0-cj, :M0-ci]
        
        if cy1.size > 0:
            mcy11 = np.repeat(cy1[:, np.newaxis], ne_v.shape[1], axis=1)
            psi_yx11 = mcy11 + cumtrapz(ne_v, axis=1, initial=0) * ne_di
            psi_yx[cj:, ci:] = psi_yx11[:N0-cj, :M0-ci]
    
    # Northwest quadrant
    if cj < N0 - 1 and ci > 0:
        nw_u = u[cj:, :ci+1]
        nw_v = v[cj:, :ci+1]
        nw_dj = dj[lx2:, :ly2]
        nw_di = di[lx2:, :ly2]
        
        if cx2.size > 0:
            mcx21 = np.repeat(cx2[np.newaxis, ::-1], nw_u.shape[0], axis=0)
            psi_xy21 = mcx21 - cumtrapz(nw_u[:, ::-1], axis=0, initial=0) * nw_dj[:, ::-1]
            psi_xy[cj:, :ci+1] = psi_xy21[:N0-cj, :ci+1][:, ::-1]
        
        if cy1.size > 0:
            mcy12 = np.repeat(cy1[:, np.newaxis], nw_v.shape[1], axis=1)
            psi_yx21 = mcy12 - cumtrapz(nw_v[:, ::-1], axis=1, initial=0) * nw_di[:, ::-1]
            psi_yx[cj:, :ci+1] = psi_yx21[:N0-cj, :ci+1][:, ::-1]
    
    # Southeast quadrant
    if cj > 0 and ci < M0 - 1:
        se_u = u[:cj+1, ci:]
        se_v = v[:cj+1, ci:]
        se_dj = dj[:lx2, ly2:]
        se_di = di[:lx2, ly2:]
        
        if cx1.size > 0:
            mcx12 = np.repeat(cx1[np.newaxis, :], se_u.shape[0], axis=0)
            psi_xy12 = mcx12 + cumtrapz(se_u[::-1, :], axis=0, initial=0) * se_dj[::-1, :]
            psi_xy[:cj+1, ci:] = psi_xy12[::-1, :M0-ci][:cj+1, :]
        
        if cy2.size > 0:
            mcy21 = np.repeat(cy2[::-1, np.newaxis], se_v.shape[1], axis=1)
            psi_yx12 = mcy21 + cumtrapz(se_v[::-1, :], axis=1, initial=0) * se_di[::-1, :]
            psi_yx[:cj+1, ci:] = psi_yx12[::-1, :M0-ci][:cj+1, :]
    
    # Southwest quadrant
    if cj > 0 and ci > 0:
        sw_u = u[:cj+1, :ci+1]
        sw_v = v[:cj+1, :ci+1]
        sw_dj = dj[:lx2, :ly2]
        sw_di = di[:lx2, :ly2]
        
        if cx2.size > 0:
            mcx22 = np.repeat(cx2[np.newaxis, ::-1], sw_u.shape[0], axis=0)
            psi_xy22 = mcx22 + cumtrapz(sw_u[::-1, ::-1], axis=0, initial=0) * sw_dj[::-1, ::-1]
            psi_xy[:cj+1, :ci+1] = psi_xy22[::-1, ::-1][:cj+1, :ci+1]
        
        if cy2.size > 0:
            mcy22 = np.repeat(cy2[::-1, np.newaxis], sw_v.shape[1], axis=1)
            psi_yx22 = mcy22 - cumtrapz(sw_v[::-1, ::-1], axis=1, initial=0) * sw_di[::-1, ::-1]
            psi_yx[:cj+1, :ci+1] = psi_yx22[::-1, ::-1][:cj+1, :ci+1]
    
    # Average the two integration paths
    psi = (psi_xy + psi_yx) / 2 * mask
    
    # Enlarge PSI into land by 1 pixel
    mask[np.isnan(mask)] = 0
    
    for i in range(N0):
        for j in range(M0):
            if mask[i, j] == 0:
                i_min, i_max = max(i-1, 0), min(i+1, N0-1)
                j_min, j_max = max(j-1, 0), min(j+1, M0-1)
                
                if np.sum(mask[i_min:i_max+1, j_min:j_max+1]) > 0:
                    psi1 = psi[i_min:i_max+1, j_min:j_max+1]
                    psi[i, j] = np.nanmean(psi1)
    
    return psi


def detect_centers(
    x: np.ndarray,
    y: np.ndarray,
    mask: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    ssh: Optional[np.ndarray],
    fields: DetectionFields,
    params: AMEDAParams,
    derived_params: AMEDADerivedParams,
    step: int
) -> Tuple[EddyCenter, EddyCenter]:
    """
    Detect potential eddy centers using LNAM and LOW fields.
    
    Parameters
    ----------
    x, y : ndarray
        Grid coordinates
    mask : ndarray
        Ocean mask
    u, v : ndarray
        Velocity components
    ssh : ndarray or None
        Sea surface height
    fields : DetectionFields
        Computed detection fields
    params : AMEDAParams
        Algorithm parameters
    derived_params : AMEDADerivedParams
        Derived parameters
    step : int
        Time step
    
    Returns
    -------
    centers0 : EddyCenter
        All LNAM maxima found
    centers : EddyCenter
        Validated centers with proper streamlines
    """
    print(f"Finding potential centers at step {step}...")
    
    # Initialize centers
    centers0 = EddyCenter(step=step)
    centers = EddyCenter(step=step)
    
    # Get fields
    OW = fields.LOW
    LNAM = fields.LNAM
    LOW = np.abs(LNAM)
    LOW[OW >= 0 | np.isnan(OW)] = 0
    
    # Find contours at threshold K
    fig, ax = plt.subplots(figsize=(1, 1))
    if params.grid_reg:
        CS = ax.contour(x[0, :], y[:, 0], LOW, levels=[params.K])
    else:
        CS = ax.contour(x, y, LOW, levels=[params.K])
    plt.close(fig)
    
    # Scan each LNAM contour
    if hasattr(CS, 'allsegs'):
        # New matplotlib API - iterate through segments
        for level_idx, level in enumerate(CS.levels):
            for seg in CS.allsegs[level_idx]:
                if len(seg) < params.n_min:
                    continue
                
                xv = seg[:, 0]
                yv = seg[:, 1]
                
                # Find points inside contour
                inside = in_polygon(x, y, xv, yv)
                
                # Mask LNAM inside contour
                Lm = LNAM.copy()
                Lm[~inside] = np.nan
                
                # Find maximum absolute LNAM inside contour
                if np.any(mask * inside > 0) and np.nanmax(np.abs(Lm)) != 0:
                    max_val = np.nanmax(np.abs(Lm))
                    max_mask = (np.abs(Lm) == max_val)
                    
                    if np.any(max_mask & (mask == 1)):
                        # Get coordinates of maximum
                        j_idx, i_idx = np.where(max_mask & (mask == 1))
                        
                        if len(j_idx) > 0:
                            j_idx = j_idx[0]
                            i_idx = i_idx[0]
                            
                            xLmax = x[j_idx, i_idx]
                            yLmax = y[j_idx, i_idx]
                            
                            # Check latitude constraint
                            if not params.grid_ll or abs(yLmax) > params.lat_min:
                                # Check if not already detected
                                already_found = False
                                for k in range(len(centers0.x)):
                                    if centers0.x[k] == xLmax and centers0.y[k] == yLmax:
                                        already_found = True
                                        break
                                
                                if not already_found:
                                    centers0.type.append(np.sign(Lm[j_idx, i_idx]))
                                    centers0.x.append(xLmax)
                                    centers0.y.append(yLmax)
                                    centers0.i.append(i_idx)
                                    centers0.j.append(j_idx)
    
    print(f"  -> {len(centers0.x)} max LNAM found at step {step}")
    
    if len(centers0.x) == 0:
        print(f"!!! WARNING !!! No LNAM extrema found at step {step}")
        return centers0, centers
    
    # Validate centers with streamline criteria
    print(f"  Validating centers with streamline criteria...")
    
    # Use interpolated parameters if needed
    if params.resol > 1:
        bxi = derived_params.bxi
        Dxi = derived_params.Dxi
        Rdi = derived_params.Rdi
        f_i = derived_params.fi
        g = 9.8  # gravity
    else:
        bxi = derived_params.bx
        Dxi = derived_params.Dx
        Rdi = derived_params.Rd
        f_i = derived_params.f
        g = 9.8
    
    # Process each max LNAM
    validated_indices = []
    
    for ii in range(len(centers0.x)):
        C_I = centers0.i[ii]
        C_J = centers0.j[ii]
        xy_ci = centers0.x[ii]
        xy_cj = centers0.y[ii]
        
        # Box size around center
        bx = int(bxi[C_J, C_I])
        Dx = abs(Dxi[C_J, C_I])
        Rd = abs(Rdi[C_J, C_I])
        f = abs(f_i[C_J, C_I])
        
        # Extract local area
        j_min = max(C_J - bx, 0)
        j_max = min(C_J + bx, x.shape[0] - 1)
        i_min = max(C_I - bx, 0)
        i_max = min(C_I + bx, x.shape[1] - 1)
        
        xx = x[j_min:j_max+1, i_min:i_max+1]
        yy = y[j_min:j_max+1, i_min:i_max+1]
        mk = mask[j_min:j_max+1, i_min:i_max+1]
        uu = u[j_min:j_max+1, i_min:i_max+1]
        vv = v[j_min:j_max+1, i_min:i_max+1]
        
        # Find center index in local domain
        cj_local = C_J - j_min
        ci_local = C_I - i_min
        
        # Compute streamfunction
        streamlines_found = False
        radius = []
        
        if params.type_detection == 1 or params.type_detection == 3:
            # Compute psi from velocity
            psi1 = compute_psi(xx, yy, mk, uu * f / g * 1e3, vv * f / g * 1e3,
                              ci_local, cj_local, params.grid_ll)
            
            # Determine contours to scan
            H = np.arange(np.floor(np.nanmin(psi1)), 
                         np.ceil(np.nanmax(psi1)), 
                         params.DH)
            if len(H) > params.nH_lim:
                H = H[:params.nH_lim]
            
            # Get contours
            if len(H) > 0:
                fig, ax = plt.subplots(figsize=(1, 1))
                CS1 = ax.contour(xx[0, :], yy[:, 0], psi1, levels=H)
                plt.close(fig)
                
                # Scan streamlines
                if hasattr(CS1, 'allsegs'):
                    for level_idx, level in enumerate(CS1.levels):
                        for seg in CS1.allsegs[level_idx]:
                            if len(seg) < params.n_min:
                                continue
                            
                            xdata = seg[:, 0]
                            ydata = seg[:, 1]
                            
                            # Check if closed and contains center
                            is_closed = (xdata[0] == xdata[-1] and ydata[0] == ydata[-1])
                            contains_center = in_polygon(np.array([xy_ci]), 
                                                        np.array([xy_cj]), 
                                                        xdata, ydata)[0]
                            
                            if is_closed and contains_center:
                                # Calculate radius
                                R, _, _, _ = mean_radius(
                                    np.vstack([xdata, ydata]),
                                    params.grid_ll
                                )
                                radius.append(R[0])
                                
                                # Check if we have 2 streamlines with proper size
                                if len(radius) >= 2:
                                    if (radius[-1] >= params.nRmin * Dx and 
                                        radius[-1] <= params.nR_lim * Rd):
                                        streamlines_found = True
                                        break
                        
                        if streamlines_found:
                            break
        
        # Check SSH contours if needed
        if not streamlines_found and params.type_detection >= 2 and ssh is not None:
            sshh = ssh[j_min:j_max+1, i_min:i_max+1]
            
            if not np.all(np.isnan(sshh)):
                # Determine contours to scan
                Hs = np.arange(np.floor(np.nanmin(sshh)), 
                              np.ceil(np.nanmax(sshh)), 
                              params.DH)
                if len(Hs) > params.nH_lim:
                    Hs = Hs[:params.nH_lim]
                
                # Get contours
                if len(Hs) > 0:
                    fig, ax = plt.subplots(figsize=(1, 1))
                    CS2 = ax.contour(xx[0, :], yy[:, 0], sshh, levels=Hs)
                    plt.close(fig)
                    
                    # Scan SSH contours
                    if hasattr(CS2, 'allsegs'):
                        for level_idx, level in enumerate(CS2.levels):
                            for seg in CS2.allsegs[level_idx]:
                                if len(seg) < params.n_min:
                                    continue
                                
                                xdata = seg[:, 0]
                                ydata = seg[:, 1]
                                
                                # Check if closed and contains center
                                is_closed = (xdata[0] == xdata[-1] and ydata[0] == ydata[-1])
                                contains_center = in_polygon(np.array([xy_ci]), 
                                                            np.array([xy_cj]), 
                                                            xdata, ydata)[0]
                                
                                if is_closed and contains_center:
                                    # Calculate radius
                                    R, _, _, _ = mean_radius(
                                        np.vstack([xdata, ydata]),
                                        params.grid_ll
                                    )
                                    radius.append(R[0])
                                    
                                    # Check if we have 2 streamlines with proper size
                                    if len(radius) >= 2:
                                        if (radius[-1] >= params.nRmin * Dx and 
                                            radius[-1] <= params.nR_lim * Rd):
                                            streamlines_found = True
                                            break
                            
                            if streamlines_found:
                                break
        
        if streamlines_found:
            validated_indices.append(ii)
            print(f"   Validated max LNAM {ii} with 2 streamlines at step {step}")
    
    # Export validated centers
    for idx in validated_indices:
        centers.type.append(centers0.type[idx])
        centers.x.append(centers0.x[idx])
        centers.y.append(centers0.y[idx])
        centers.i.append(centers0.i[idx])
        centers.j.append(centers0.j[idx])
    
    print(f"  -> {len(centers.x)} potential centers found")
    print(f"    ({len(centers0.x) - len(centers.x)} max LNAM removed)")
    print("")
    
    return centers0, centers


# Main processing function
def process_centers(
    nc_file: str,
    step: int,
    fields: DetectionFields,
    params: AMEDAParams,
    derived_params: AMEDADerivedParams
) -> Tuple[EddyCenter, EddyCenter]:
    """
    Main function to detect eddy centers.
    
    Parameters
    ----------
    nc_file : str
        Path to NetCDF file
    step : int
        Time step to process
    fields : DetectionFields
        Computed detection fields
    params : AMEDAParams
        Algorithm parameters
    derived_params : AMEDADerivedParams
        Derived parameters
    
    Returns
    -------
    centers0 : EddyCenter
        All LNAM maxima found
    centers : EddyCenter
        Validated centers
    """
    import xarray as xr
    
    # Load dataset
    ds = xr.open_dataset(nc_file)
    
    # Load fields
    x, y, mask, u, v, ssh = load_fields_aviso(
        ds, step, params, params.resol, params.deg
    )
    
    # Set velocities to 0 in land for psi calculation
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    
    # Detect centers
    centers0, centers = detect_centers(
        x, y, mask, u, v, ssh,
        fields, params, derived_params,
        step
    )
    
    ds.close()
    
    return centers0, centers


if __name__ == "__main__":
    # Test with sample data
    from ameda.params import AMEDAParams, compute_derived_params
    from ameda.fields import process_fields
    
    DATASET = "/Users/gianlucacalo/Desktop/projects/ocean_surface_levels/data/raw/cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.0625deg_P1D.nc"
    
    # Load test data
    import xarray as xr
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
    
    # Process fields first
    fields = process_fields(DATASET, step, params, derived_params)
    
    # Process centers
    centers0, centers = process_centers(DATASET, step, fields, params, derived_params)
    
    print(f"Centers detected successfully:")
    print(f"  Total LNAM maxima: {len(centers0.x)}")
    print(f"  Validated centers: {len(centers.x)}")
    if len(centers.x) > 0:
        print(f"  Center types: {centers.type}")
        print(f"  First center at: ({centers.x[0]:.3f}, {centers.y[0]:.3f})")
    
    ds.close()