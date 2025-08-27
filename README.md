# AMEDA: Angular Momentum Eddy Detection and tracking Algorithm
AMEDA, **Angular Momentum Eddy Detection and tracking Algorithm**, is an algorithm that follows Nencioli algorithm in Matlab (Nencioli *et al.* 2010) to detect and to track moving eddy structures in the ocean with a strong geostrophic signature. The method records time series of center positions of the eddies, along with their size and intensity; it gives also an history of the merging and splitting events.  

It is a hybrid detection method based on the computation of Local Normalised Angular Momentum developed at the LMD by Nadia Mkhinnin and Alexander Stegner (Mkhinnin *et al.* 2014) which allows to locate a center of an eddy of any intensity and propose a contour based on the maximum radial velocity.  

The advantage of the algorithm resides in its capacity to analyse various gridded velocity fields of oceanic currents (velocimetry imagery, high-frequency radar, satellite, numerical model – see figure 1) without fine tuning of the parameters. It has been calibrated on SSALTO/DUACS products, velocities from PIV laboratory imagery and ROMS model output (Le Vu *et al.* 2018).  

The code is freely available and regularly updated on github to improve its efficiency and to correct bugs. Feel free to test it on your own field and ask for any advice if necessary.

# First steps with AMEDA
## Requirements
To run AMEDA at its best performance you need a Matlab License with at least version 7.3 and some of Matlab’s toolboxes:

- Matlab tested versions (R009a, R2014a and 2017a).

- Toolboxes:
  - **Curve Fitting Toolbox** to specifically compute the shape coefficient of the V-R profile. You must deactivate the `streamlines` key if you don’t want to use this toolbox.
  - **Parallel Computing Toolbox** to run the parallelised version of the code on a single node (tested up to 32 CPUs for CROCO outputs fields on datarmor cluster).
  - **Distributed Computing Server** should be necessary to run the parallelised version of the code on the un cluster. This has not been tested yet, but here is the documentation: [Mathworks documentation](https://fr.mathworks.com/support/product/DM/installation/ver_current.html).
  - **Matlab Compiler.** This tool allows to compile the AMEDA code and produce an executable which does not need a Matlab licence to run. This toolbox helps particularly in analyzing a long simulation when there are a limited number of available licences.
  - **m_map** (Pawlowicz, 2020) or any tools for projecting your results on a map.

## How to deploy AMEDA?

In the main directory called **AMEDA** you can find the main routines which have to be modified by user to run the code. You may place this directory and all its content on your usual Matlab directory.  

From the downloaded files you also have to untar **Rossby_radius.tar** which contains matrixes of the Rossby radius, an essential part of the parameters computation automatically computed by AMEDA during the early steps of the algorithm. Figure 3 shows the detailed deployment of the code.

Routines of the code developed for AMEDA are located in the **sources** subdirectory and routines taken and adapted from [Mathworks](https://fr.mathworks.com/) or from model post-processing tools are located in the **tools** subdirectory, along with a list of associated licences in the `licences.txt` file.  

👉 To increase AMEDA speed performance, we suggest you to build **MEX files** according to your Matlab environment from the few files provided in C or Fortran language.
### Routine and components of the downloaded code AMEDA
```
AMEDA
├── Contents.txt
├── keys_sources_AVISO_DYNED_MED_adt.m
├── keys_sources_AVISO_tuto_dt.m
├── keys_sources_AVISO_tuto_nrt.m
├── keys_sources_AVISO.m
├── keys_sources_CROCO_2D_tuto_CRE.m
├── keys_sources_CROCO_test_deg2.m
├── keys_sources_NEMO_test.m
├── MAIN_AMEDA_AVISO_tuto_dt.m
├── MAIN_AMEDA_AVISO_tuto_nrt.m
├── MAIN_AMEDA_CROCO_2D_tuto_CRE.m
├── MAIN_AMEDA_CROCO_tuto_test.m
├── MAIN_AMEDA_DYNED_MED.m
├── MAIN_AMEDA_MED_multi.m
├── MAIN_AMEDA_MED_nopool.m
├── MAIN_AMEDA_NEMO_multi.m
├── MAIN_AMEDA_NEMO_test.m
├── plot # Visualization routines
│   ├── plot_shapes_example1.m
│   ├── plot_shapes_example2.m
│   ├── plot_tracking_example3.m
│   └── plot_tracking_example4.m
├── README.md
├── Rossby_radius # First baroclinic Rossby radius data files
│   ├── global_KNTN.mat
│   ├── read_global_Rossby_radius.m
│   ├── Rossby_radius.mat
│   ├── global_Rossby_radius.mat
│   ├── Rossby_Radius_WOA_Barocl1.mat
│   └── rossrad.dat
├── Rossby_radius.tar
├── sources # Core algorithm implementations
│   ├── compute_best_fit.m
│   ├── compute_curve.m
│   ├── compute_ellip.m
│   ├── compute_psi.m
│   ├── concat_eddy.m
│   ├── eddy_dim.m
│   ├── integrate_vel.m
│   ├── load_fields_AVISO.m
│   ├── load_fields_CROCO_2D.m
│   ├── load_fields_CROCO.m
│   ├── load_fields_HFR.m
│   ├── load_fields_HYCOM.m
│   ├── load_fields_NEMO.m
│   ├── load_fields_PIV.m
│   ├── load_fields_ROMS.m
│   ├── make_netcdf_from_tracks.m
│   ├── max_curve.m
│   ├── mean_radius.m
│   ├── min_dist_shapes.m
│   ├── mod_eddy_centers.m # Detects eddy centers using LNAM maxima
│   ├── mod_eddy_params.m # Loads configuration and computes detection parameters
│   ├── mod_eddy_shapes.m # Computes eddy dimensions and contours
│   ├── mod_eddy_tracks_nopool.m # Links eddies across time steps
│   ├── mod_eddy_tracks_pool.m
│   ├── mod_fields.m # Computes detection fields (KE, vorticity, OW, LNAM)
│   ├── mod_init.m  # Initializes data structures and output files
│   ├── mod_merging_splitting.m # Detects merging/splitting events
│   └── scan_lines.m
├── start.m
├── tools # Utility functions and third-party tools
   ├── assignmentoptimal.c
   ├── assignmentoptimal.m
   ├── csf.m
   ├── fitellipse.m
   ├── get_Dx_from_ll.m
   ├── get_missing_val_2d.m
   ├── get_z_croco.m
   ├── InPolygon.c
   ├── InPolygon.m
   ├── license.txt
   ├── nanmean.m
   ├── nanstd.m
   ├── nansum.m
   ├── nanvar.m
   ├── sw_dist.m
   ├── sw_dist2.m
   ├── TaubinNTN.m
   └── vinterp.m
```

Additional details:
- `MAIN_AMEDA_*.m`: Entry points for different data sources
- `keys_sources_*.m`: Configuration files for each data source

# AMEDA algorithm workflow in details

## 1. User routines

In the AMEDA folder itself, you will find 2 kinds of Matlab files and also the file `Contents.txt` which describes the workflow followed by the keys routines during a full analysis by AMEDA. The 2 others kind of routines that user needs to modify for running AMEDA are:

- **MAIN_AMEDA\*.m** is the Matlab routine you launch to run AMEDA. You need a good understanding of this routine to be able to modify `MAIN_AMEDA.m` files for your own usage. The postfix `multi` is used when there is a need to time split a very long input file into more than one analysis in order to prevent a possible crash during a long computation time.  

  For a multi-years run for example, `AMEDA_MAIN_*_multi.m` will run, for each year, the detection analysis and will record the results in files with the year as a postname accordingly. Then `AMEDA_MAIN_*_multi.m` will run the tracking by loading the concatenated output files containing the yearly results.  

  The postfixes `nopool` and `pool` define two versions of the same routine in which a `parfor loop` or a simple `for loop` is used, respectively.

- **keys_sources\*.m** contain definitions of paths and fields names, as well as setting of keys and options depending on the simulation you want to execute. These distinct routines have been created at different moments of the code development. Thus, some of them may have missing keys or parameters that you will have to add if you want to use these AMEDA sources.

Details on modifying and using these routines are described extensively in the next section, **running AMEDA**.

---

## 2. Parameters computation

Running one of the `MAIN_AMEDA*.m` with its associated `keys_sources*.m` begins by the parameters computation (`mod_eddy_params.m`). Many parameters like the Coriolis parameter or deformation radius are computed in this routine. You can consult the different parameters imposed or computed in Table 1.  

It is possible to modify these parameters by changing the value directly in the routine or you can also remove some from the routine and define them in your user routine `keys_sources*.m`. We will describe in detail what user parameters need to be defined in the section **running AMEDA**.  

All the 2D-fields and parameters necessary for running AMEDA as well as the different paths of data and results defined in `keys_sources*.m` are recorded in the file `param_eddy_tracking.mat` at the end of this stage. This file is recorded in the folder defined by `path_out` in `keys_sources*.m`.

⚠️ **Warning**: Defining correctly the `path_out` variable is very important because it will be added to the Matlab path as a ‘global variable’ and scanned by the code when AMEDA loads a file. For this reason, if the value of `path_out` must change due to changes in the configuration or the simulation, please start a new Matlab session.

---

## 3. Main routines

After the parameter have been computed as described above, and that the Matlab structures pre-allocation or update (`mod_init.m`) is successfully completed, the computational part of AMEDA algorithm will be executed. The execution of the AMEDA algorithm utilizes the four following routines:

1. **mod_fields.m**  
   Computes, through a finite spatial elements method, the 2D fields of the following variables: kinetic energy, divergence, vorticity, Okubo-Weiss, Local Okubo-Weiss (LOW) and Local Normalised Angular Momentum (LNAM). Among them, the two latter (LOW and LNAM) will be used to detect eddy centers by the next script. LOW and LNAM and the others fields are included in a structure array called `{detection_fields(t)}` and save in the `fields.mat` file located in `path_out`.

2. **mod_eddy_centers.m**  
   Detects the potential eddies centers present in the domain using the LNAM and LOW fields. In fact, when the precision of the grid is lower than the deformation radius like for AVISO fields, AMEDA interpolates the LNAM and LOW on a thinner grid by a factor `resol` (cf. Table 1). This interpolation is made to locate centers at a precision higher than the size of the typical deformation radius of the domain.  

   - Firstly, the position of the maxima `max[LNAM(LOW<0)>K]` are resolved, with `K` referring to the LNAM threshold defined in `mod_eddy_params.m`.  
   - Then, the potential centers saved are these max LNAM surrounded by at least two closed contours of the streamfunction (`psi`) or `ssh`.  

   Potential eddy centers are saved/updated as the structure array `{center(t)}` in the file `eddy_centers.mat` in `path_out`.

3. **mod_eddy_shapes.m**  
   Computes the shapes (if any) of eddies identified by their potential centers in the previous stage. Characteristic shapes are closed contours defined by the maximum mean tangential velocity along the closed streamlines around the center, hereby called as **speed radius**.  

   These contours have a maximum (i.e. `nR_lim`, cf. Table 1) and a minimum (i.e. `nRmin`, cf. Table 2b) size. The corresponding identified eddy centers are saved as `{centers2(t)}` in the file `eddy_centers.mat` in `path_out` and speed radius together with other features as `{shapes1(t)}` in `eddy_shapes.mat` in `path_out`.  

   The routine also records the outermost contours in `{shapes1(t)}` and defines contours which include two different eddy centers, called **dual eddies**, as `{shapes2(t)}` in the same `eddy_shapes.mat` file.

4. **mod_eddy_tracks.m**  
   Performs eddy tracking through the eddy centers positions and shapes features computed as described above. Eddy tracks are inter-connected by comparing the detected eddy fields at successive time steps.  

   An eddy at time *t* is assumed to be the new position of an eddy of the same type detected at time *t-dt*, only if the distance *d* between the two centers is smaller than the distance D:  
   D = V_eddy * (1+dt)/2 + r_maxN + r_maxm

   where `r_maxN` is the mean radius of each eddy time-averaged on the last `D_stp` tracked time steps (cf. Table 1) and `r_maxm` the new eddy radius to be tested at `t`. If an eddy does not connect any new detection after a number of time steps `Dt` (cf. Table 1), it is considered dissipated and the track *n* ends. In case two or more eddies are found within the same area defined by D, the track is connected to the centers which minimize the cost function of an assignment matrix considering all the past centers N at `t-dt` which are not already dissipated and the new ones M at `t`. The cost function is a NxM matrix of the `C_nm` elements:

   C_nm = sqrt( d/D² + δR/R_max² + δRo/Ro² + dt/Dt/2² )
   Unfiltered tracked eddies are saved/updated as the structure array `{traks(n)}` in `eddy_tracks.mat` in `path_out` where *n* is the index number of the recorded track.

On the first three stages, the routines are not time step dependent, and thus the parallelization is applied in the `MAIN_AMEDA*.m` routine to compute many time steps at the same time thanks to the `parfor loop`.
On the last stage, tracking is dependent on the previous steps and thus the `parfor loop` cannot be applied with the time step. A `parfor loop` is included inside `mod_eddy_tracks_pool.m` at the cost function calculation. Unfortunately, the gain in terms of computation time is not important when the number of eddies tracked increases too much (i.e. in a big domain or a long run) due to the increase of the allocation time to the pools with the size of the tracks structure.   While we have not modified the code to limit the size of this matrix we advice you to only use `mod_eddy_tracks_nopool.m`.


## 4. Flags and filtering

At the end of the AMEDA process, the unfiltered **tracks** structure is analysed to resolve merging and splitting events from eddy tracks. A time filtering is also applied, as two times the turnover time if `cut_off` is 0 (by default) or the value of `cut_off` otherwise.  

In order to avoid time filtering and to save every eddy tracks whatever their life length you must set `cut_off` to the value `dps` which is settled in your `keys_sources*.m` file (cf. Table 2b).  

The result of this flagging (merging or splitting) and filtering process is recorded as `{tracks2(n)}` in `eddy_tracks2.mat` in `path_out`.


# Running the AMEDA walkthrough

## 1. Build input files

Apart from the two user routines (`MAIN_AMEDA*.m` and `keys_sources*.m`) that you need to adapt to your purpose, you first need to prepare your input files in agreement with the AMEDA specifications. The input interface is taken in charge by the routines `load_fields*.m`. Many specifications are allowed, but only some are absolutely necessary:

- Every processed file must be of reasonable size for your computer buffer memory (typically few Go <10 Go). Keep this in mind when producing your input files.
- Time is given in terms of unitless and equally spaced steps as the time step duration is fixed in the `keys_sources*.m` script.
- The dimensions of the fields output of `load_fields*.m` needs to be `[field] = (x | longitude, y | latitude)`.
- The velocity units as output of `load_fields*.m` must be in m.s⁻¹.

You could change these specifications. In that case you will need to produce your own `load_fields*.m` with your own requirement. For example, create a `load_fields_model.m` to read all the different fields from the same file like a 3D output model containing all the fields.  

You might refer to the proper `load_fields_source.m` routine belonging to your source (AVISO, NEMO, …) and you can even modify it to properly build your input files.  

---

## 2. Choose or prepare the `rossby_radius` file

The first assumption when using AMEDA is that you want to detect eddies of size close to the first baroclinic Rossby radius of deformation (**Rd**). Indeed, Rd is the first needed parameter from which the other parameters are calculated in `mod_eddy_params.m`.  

- AMEDA code includes a global Rd 2D-field in the file `Rossby_radius/global_Rossby_radius.mat`, built from Chelton et al. 1998, that you can use for any area in the ocean.
- This global field is somewhat coarse (1°x1°). If you work on the Mediterranean sea, you should use `Rossby_radius/Rossby_radius.mat` (resolution 1/8°x1/8°).  
- You might also create your own Rd `.mat` file by following the script `Rossby_radius/read_global_Rossby_radius.m` which creates a `.mat` file from a `.dat` file.

⚠️ **Warning**: You should pay attention to use a continuous longitude grid in your area of interest. For instance, avoid longitude ranges of 180°E to -180°W or 359°E to 0°.

---

## 3. Choose and set `MAIN_AMEDA`

`MAIN_AMEDA*.m` are the master routines executing all the stages described in the previous section **AMEDA workflow**. Each `MAIN_AMEDA*.m` depicts relevant schema depending on the pre-processed input files available and the specific configuration you want to apply.  

- Once you have your input netcdf files in the right format and units you can run `MAIN_AMEDA`.
- You need to add the paths `AMEDA/`, `sources/`, `tools/`, `Rossby_radius/` and `m_map/` in the Matlab path with the command `addpath`.
- This can be done through the `start.m` script, as indicated in the `MAIN_AMEDA*.m` routines.

In practice, to run AMEDA, you must choose and copy a `MAIN_AMEDA*.m` and a `keys_sources*.m` file, based on your problem setup and hardware availability, and rename them. For example:  

- Copy and create `MAIN_AMEDA_AVISO_test.m` and `keys_sources_AVISO_test.m` from `MAIN_AMEDA_MED.m` and `keys_sources_AVISO_MED.m`.

### The choice of the `MAIN_AMEDA*.m` script depends on:
- **Source availability**: is it already known by AMEDA?  
- **Parallelization**: does your hardware setup allow for using the `parpool` function?  
- **Steps/file size**: do you have more than 1000 steps or a big input file more than 10Go?  

### Parameters to modify in `MAIN_AMEDA*.m`:
- `source`: drives the `load_fields_source.m` used to read fields (e.g. `AVISO`, `NEMO`, `HFR`, `CROCO`).
- `keys`: specifies the configuration/domain (e.g. `MED`, `test`).
- `update`: integer to follow an existing analysis or complete one.
- `stepF`: optional integer to stop after N steps (for quick tests).
- `deg`: degrade input field by subsampling (default is 1, no degradation).
- `cpus`: number of CPUs (uses **Parallel Computing Toolbox**).
- `list`: array of starting steps (used in `MAIN_AMEDA_*_multi.m`).
- `name`: string used as postname in output files.

---

## 4. Set `keys_sources`

`keys_sources*.m` is the other user routine interface where you set your paths and your option keys belonging to the configuration. It is called at the beginning of the algorithm to describe the configuration and to add the results directory, `path_out`, to the Matlab path.  

It is recommended to name `keys_sources_source_keys.m` with exactly the same names as the variables `source` and `keys` set in your `MAIN_AMEDA*.m`.

---

## Table 2a: Paths, files and names definitions in `keys_sources*.m`

| Paths       | Definition                          | Usage or value | Remarks |
|-------------|--------------------------------------|----------------|---------|
| source      | Source of the input | Calls proper `load_fields_source.m` | Same as in `MAIN_AMEDA*.m` |
| config      | Configuration | Defines path and input file name | Different from `keys` or `domain` |
| sshype      | Type of ssh (adt, sla, model, …) | Distinguish different sources | Optional, flexible |
| runname     | Sensibility test name | For test runs | Optional |
| postname    | Specific analysis | For cropped common input file | |
| path_in     | Input files directory | Used to open input files | Use previous names |
| path_result | Output files common directory | Define path for results | Optional but convenient |
| path_out    | Output files specific directory | Save outputs + param file | Important: added to Matlab path |
| path_tracks | Satellite tracks directory | For plotting | Useful for satellite comparisons |
| path_data   | In situ data directory | For drifters or ARGO | |
| path_rossby | Rossby radius directory | Open First Baroclinic Rossby radius | Must exist or be created |
| nc_dim      | Input x,y,mask field file | Netcdf with grid data I x J | |
| nc_u, nc_v  | Input u, v velocity fields | Netcdf file with u,v values I x J x L | |
| nc_ssh      | Input ssh field file | Netcdf file with ssh data I x J x L | |
| x,y,m_name, … | Field names in nc files | Match standard input | Adaptable |
| mat_Rd      | Rossby radius file absolute name | e.g. `global_Rossby_radius.mat` | Or Med sea file |
| name_Rd     | Rossby radius name in `mat_Rd` | `lon_Rd`, `lat_Rd` | Must be 2D |

---

## Table 2b: Parameters definition and keys activation in `keys_sources*.m`

| Parameter       | Definition | Value | Usage | Remarks |
|-----------------|------------|-------|-------|---------|
| Rd_typ | Rossby radius typical value | 12 km (Med sea) | Used to limit Rd to 1/3 of Rd_typ | Default AVISO |
| nRmin | Minimal size for Rmax | 0.5 Dx | Min radius of recorded eddies | Default AVISO |
| T | Period | 24\*3600 sec | For Coriolis param | Default ocean |
| dps | Turnover per step | 1 / step | For tracking | Daily step default |
| level | Vertical level | 1 | For >1 level netcdf | |
| grid_ll | Grid type | 0 or 1 | Cartesian (km) or Earth coords (deg) | |
| grid_reg | Regular grid? | 0 or 1 | 0 irregular, 1 regular | |
| type_detection | Detection type | 1, 2, 3 | Velocity fields, SSH, or both | |
| extended_diags | Extra diagnostics | 1 | 0 not saved, 1 saved | Post-processing |
| streamlines | Compute V-R profiles | 1 | Needs Curve Fitting Toolbox | |
| daystreamfunction | Steps of V-R | 1:stepF | Time steps for streamlines key | |
| periodic | Periodic fields | 0 or 1 | Repeat field at edges | |
| nrt | Near real time | 1 | 0 not activated, 1 activated | |

---

## 5. Launch AMEDA

You can run the `MAIN_AMEDA*.m` routine all at once or partially by means of a Matlab interface. Otherwise, you could launch the script by detaching it from your connection while printing the output log in a specific file, if necessary.  

To run your `MAIN_AMEDA*.m` in such a detached job manner, use for example:

```bash
~/MATLAB$ nohup matlab -nodestkop -nodisplay < ./AMEDA/MAIN_AMEDA_*.m > job_out.txt &
