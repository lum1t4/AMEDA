"""
Test script for AMEDA implementation
"""

import numpy as np
import xarray as xr
from ameda.params import AMEDAParams, compute_derived_params
from ameda.fields import process_fields

# Configuration
DATASET = "/Users/gianlucacalo/Desktop/projects/ocean_surface_levels/data/raw/cmems_obs-sl_eur_phy-ssh_my_allsat-l4-duacs-0.0625deg_P1D.nc"

print("=" * 60)
print("Testing AMEDA Python Implementation")
print("=" * 60)

# Load dataset
print("\n1. Loading dataset...")
ds = xr.open_dataset(DATASET)
print(f"   Dataset shape: time={len(ds.time)}, lat={len(ds.latitude)}, lon={len(ds.longitude)}")
print(f"   Variables: {list(ds.data_vars.keys())}")

# Configure parameters
print("\n2. Setting up parameters...")
params = AMEDAParams(
    x_name="longitude",
    y_name="latitude",
    u_name="ugos",
    v_name="vgos",
    s_name="sla",
    type_detection=1,  # Use velocity only for simplicity
    resol=1,  # No interpolation
    K=0.3  # Lower threshold for more detections
)

# Get first time step
print("\n3. Extracting first time step...")
step = 0
data_step = ds.isel(time=step)
lon = data_step.longitude.values
lat = data_step.latitude.values
ugos = data_step[params.u_name].values
vgos = data_step[params.v_name].values
sla = data_step[params.s_name].values

# Check for valid data
valid_mask = ~np.isnan(ugos)
print(f"   Valid data points: {np.sum(valid_mask)} / {ugos.size}")
print(f"   Velocity range: u=[{np.nanmin(ugos):.3f}, {np.nanmax(ugos):.3f}] m/s")
print(f"                   v=[{np.nanmin(vgos):.3f}, {np.nanmax(vgos):.3f}] m/s")
print(f"   SSH range: [{np.nanmin(sla):.3f}, {np.nanmax(sla):.3f}] m")

# Compute derived parameters
print("\n4. Computing derived parameters...")
derived_params = compute_derived_params(lon, lat, valid_mask, ugos, vgos, params)
print(f"   Rossby radius range: [{np.nanmin(derived_params.Rd):.1f}, {np.nanmax(derived_params.Rd):.1f}] km")
print(f"   Grid spacing range: [{np.nanmin(derived_params.Dx):.1f}, {np.nanmax(derived_params.Dx):.1f}] km")
print(f"   LNAM box size range: [{np.nanmin(derived_params.b):.0f}, {np.nanmax(derived_params.b):.0f}] pixels")

# Process fields
print("\n5. Computing detection fields...")
try:
    fields = process_fields(DATASET, step, params, derived_params)
    print("   ✓ Fields computed successfully")
    print(f"   KE range: [{np.nanmin(fields.ke):.6f}, {np.nanmax(fields.ke):.6f}] m²/s²")
    print(f"   Vorticity range: [{np.nanmin(fields.vort):.2e}, {np.nanmax(fields.vort):.2e}] s⁻¹")
    print(f"   OW range: [{np.nanmin(fields.OW):.2e}, {np.nanmax(fields.OW):.2e}] s⁻²")
    print(f"   LNAM range: [{np.nanmin(fields.LNAM):.3f}, {np.nanmax(fields.LNAM):.3f}]")
    
    # Count potential centers
    LOW = fields.LOW
    LNAM = fields.LNAM
    potential = np.abs(LNAM) * (LOW < 0)
    high_lnam = potential > params.K
    print(f"   Potential centers (|LNAM| > {params.K} & LOW < 0): {np.sum(high_lnam)}")
    
except Exception as e:
    print(f"   ✗ Error: {e}")
    import traceback
    traceback.print_exc()

# Try center detection with simpler approach
print("\n6. Testing simplified center detection...")
try:
    from ameda.fields import load_fields_aviso
    
    # Load fields
    x, y, mask, u, v, ssh = load_fields_aviso(ds, step, params, params.resol, params.deg)
    
    # Find LNAM maxima
    LOW = fields.LOW
    LNAM = fields.LNAM
    
    # Simple threshold detection
    candidates = (np.abs(LNAM) > params.K) & (LOW < 0) & (mask > 0)
    n_candidates = np.sum(candidates)
    
    if n_candidates > 0:
        j_idx, i_idx = np.where(candidates)
        print(f"   Found {n_candidates} LNAM maxima above threshold")
        
        # Show first few centers
        for k in range(min(5, n_candidates)):
            print(f"   Center {k+1}: lon={x[j_idx[k], i_idx[k]]:.3f}°, "
                  f"lat={y[j_idx[k], i_idx[k]]:.3f}°, "
                  f"LNAM={LNAM[j_idx[k], i_idx[k]]:.3f}")
    else:
        print("   No candidates found")
        
except Exception as e:
    print(f"   ✗ Error: {e}")

ds.close()

print("\n" + "=" * 60)
print("Test complete!")
print("=" * 60)