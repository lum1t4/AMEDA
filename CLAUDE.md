# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview
AMEDA, **Angular Momentum Eddy Detection and tracking Algorithm**, is an algorithm that follows Nencioli algorithm in Matlab (Nencioli *et al.* 2010) to detect and to track moving eddy structures in the ocean with a strong geostrophic signature. The method records time series of center positions of the eddies, along with their size and intensity; it gives also an history of the merging and splitting events.  

It is a hybrid detection method based on the computation of Local Normalised Angular Momentum developed at the LMD by Nadia Mkhinnin and Alexander Stegner (Mkhinnin *et al.* 2014) which allows to locate a center of an eddy of any intensity and propose a contour based on the maximum radial velocity.  

The advantage of the algorithm resides in its capacity to analyse various gridded velocity fields of oceanic currents (velocimetry imagery, high-frequency radar, satellite, numerical model – see figure 1) without fine tuning of the parameters. It has been calibrated on SSALTO/DUACS products, velocities from PIV laboratory imagery and ROMS model output (Le Vu *et al.* 2018).  


see @README.md to understand the overall structure and usage of the AMEDA algorithm.

### Workflow Pattern
1. **Parameter Initialization**: `mod_eddy_params` loads configuration from `keys_sources_*.m` and computes:
   - Grid parameters (Dx, coordinates, Coriolis parameter)
   - Rossby radius (Rd) from pre-computed files
   - LNAM computation parameters (b, bx)
   - Detection thresholds (K for LNAM threshold)
   - Tracking parameters (V_eddy, Dt, cut_off)
2. **Field Computation**: `mod_fields` computes detection fields for each time step:
   - Kinetic Energy (KE): (u²+v²)/2
   - Divergence: ∂u/∂x + ∂v/∂y
   - Vorticity: ∂v/∂x - ∂u/∂y
   - Okubo-Weiss (OW): strain² + shear² - vorticity²
   - Local Okubo-Weiss (LOW): averaged OW
   - Local Normalized Angular Momentum (LNAM): key detection parameter
3. **Center Detection**: `mod_eddy_centers` finds LNAM maxima where:
   - LOW < 0 (rotation dominates strain)
   - |LNAM| > K (above threshold)
   - At least 2 closed streamfunction/SSH contours exist
4. **Shape Analysis**: `mod_eddy_shapes` determines eddy boundaries:
   - Scans closed streamlines around each center
   - Finds "speed radius" (contour with maximum mean velocity)
   - Computes radius, velocity, and ellipse parameters
5. **Tracking**: `mod_eddy_tracks` connects eddies across time:
   - Uses cost function based on distance, radius change, and time gap
   - Resolves conflicts with optimal assignment algorithm
   - Allows gaps up to Dt days
6. **Post-processing**: `mod_merging_splitting` filters and analyzes:
   - Identifies merging and splitting events
   - Removes short-lived eddies (< cut_off days)

### Required MATLAB To
olboxes
- Curve Fitting Toolbox (for shape coefficient computation)
- Parallel Computing Toolbox (for parfor loops)
- m_map (external, for map projections)

### Performance Considerations
- Use parallel computing with `cpus` parameter in configuration
- Compile MEX files for critical functions (InPolygon, assignmentoptimal)
- Test with degraded resolution using `deg` parameter before full runs
- Memory-intensive for large domains and long time series
- Parallelized versions available for fields, centers, and shapes modules
- Multi-file processing option splits long time series into yearly chunks

### Data Format Requirements
Input NetCDF files must contain:
- Velocity components: `u`, `v` (or `U`, `V`)
- Coordinates: `x`/`lon`, `y`/`lat`
- Time dimension
- Optional: SSH field for visualization

Output MAT files:
- `param_eddy_tracking.mat`: All computed parameters
- `fields.mat` / `fields_inter.mat`: Computed detection fields (KE, vort, OW, LNAM)
- `eddy_centers.mat`: Detected centers at each processing stage
- `eddy_shapes.mat`: Eddy boundaries and characteristics
- `eddy_tracks.mat`: Raw connected tracks
- `eddy_tracks2.mat`: Final filtered tracks with merging/splitting events