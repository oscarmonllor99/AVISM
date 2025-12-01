# AVISM: Algorithm for Void Identification in coSMology

## About

Developed at the Departament d'Astronomia i Astrofísica of Universitat de València by Óscar Monllor-Berbegal in collaboration with David Vallés-Pérez, Susana Planelles and Vicent Quilis. This work has been supported by the European Union NextGenerationEU (PRTR-C17.I1), the Spanish Ministerio de Ciencia e Innovación(ASFAE/2022/001 and PID2022-138855NB-C33), the Generalitat Valenciana (CIPROM/2022/49), and Óscar Monllor-Berbegal acknowledges support from Universitat de València through an Atracció de Talent fellowship.

Citation: _Monllor-Berbegal et al. 2025_ ([A&A](https://doi.org/10.1051/0004-6361/202554513)). 

## Index of contents 

1. [Brief description](#brief-description)

2. [Repository organisation](#repository-organisation)

3. [Installation](#installation)

4. [Running the code](#running-the-code)

5. [Input](#input)

6. [Code configuration](#code-configuration)

7. [Output](#output)


## Brief description

**AVISM** is a void finder algorithm combining both geometrical and dynamical information in order to identify voids within a specific input region in which a given set of particles (or a grid) have matter (density) and velocities well defined. It is able to handle cosmological simulation outputs, dark matter halo catalogues and galaxy surveys as inputs. To do so, it interpolates the data onto a uniform grid by means of an efficient parallel k-d tree implementation, allowing the code to calculate each cell particle neighbours quickly. With the density and velocity fields defined on the grid, the code performs three main steps to identify voids:

1. All cells satysfing $\delta > \delta_1$ and $\nabla \cdot \mathbf{v} > 0$ are marked as potential void centres.
2. These cells are expanded, starting from the cell with maximum divergence in descending order, until one of the following criteria is fulfilled:
    * $|\nabla\delta| > |\nabla\delta|_\text{th}$
    * $\nabla \cdot \mathbf{v} < \nabla \cdot \mathbf{v}_\text{th}$
    * $\delta > \delta_\text{2}$
3. After the previous step, the code obtains a set of overlapping cubes covering all regions potentially belonging to a void. The next step consists of merging these cubes in a volume-ordered way, starting with the biggest, such that voids are built all synchronously, starting by the core regions (biggest cubes) until reaching the boundaries (smaller cubes).

This process can be repeated at different levels of resolution using finer grids, leading to a list of voids-in-voids. The void finder characteristics, several applications (mock test, cosmological simulation data and galaxy surveys), and a thorough comparison with other void finders can be found in _Monllor-Berbegal et al. 2025 accepted, A&A_ ([ArXiv](https://arxiv.org/abs/2509.25329)). 

## Repository organisation

The source code can be found inside the `src` folder with all `.f90` files. Several `Python` tools are provided inside the `tools` folder. These are examples to handle simulated and galaxy survey data in order to prepare a proper input for the void finder. Inside `test1` one can find the source code used to create the mock voids test described in _Monllor-Berbegal et al. 2025 accepted, A&A_ ([ArXiv](https://arxiv.org/abs/2509.25329)). Finally, the `config` folder contains `voids.dat`, which is the main configuration file provided to the code. We also give an example `Makefile` and `run.sh` in order to properly compile and execute the code.

## Installation

### Download

Clone the repository into the desired directory by running:

```
git clone https://github.com/oscarmonllor99/AVISM.git
```

Then access the root directory:

```
cd AVISM-main
```

### Dependencies

Our code needs few dependencies: a Fortran compiler and OpenMP for shared-memory parallelisation.

* [gfortran](https://gcc.gnu.org/wiki/GFortran)
* [OpenMP](https://www.openmp.org/)
* (optional) [HDF5](https://www.hdfgroup.org/solutions/hdf5/)

HDF5 is only used for some kind of input data, such as reading directly from Arepo outputs.

### Make

We already provide the user with an example Makefile to compile the code. Inside the root directory, where we can find `src` and `bin`, one can compile the code by simply running

```
make COMP=1
```

We also provide a debugging option:

```
make COMP=2
```

Both options will create the `avism.x` executable inside the root directory.

Furthermore, two more compilation options are available: `HDF5` (defaults to 0) and `PERIODIC` (defaults to 0). The first specifies if the [HDF5](https://www.hdfgroup.org/) library is needed (`HDF5=1` ) to read the input data. The second tells the code to use periodic boundary conditions (`PERIODIC=1`), if needed.

## Running the code

In order to run the code, the user needs to execute `avism.x` with a suitable memory and parallelisation (OpenMP) configuration. The `run.sh` file exemplifies this, with a proper setup for most **AVISM** application cases. If a large grid is used (say $1024^3$ cells) or a huge number of particles is provided to the code (more than $10^9$), the user may need to modify some parameters like `stacksize` or `memoryuse`.

## Input data

As of today, the code allows for 4 different types of input:

1. **MASCLET input:** `Option 0`
   
   Since the code was originally conceived to find voids within the MASCLET simulation outputs, this input option is still supported. Furthermore, extending it to other adaptive-mesh refinement (AMR) codes should be straightforward. If this option is chosen, the void finder transforms the AMR grid into a uniform representation with the required resolution. Also, if this input is chosen, input files should be inside `path_to_AVISM/simu_masclet`.

2. **Particle input:** `Option 1`
   
   In order to extend its applicability, we implement the `particle.f90` and `kdtree.f90` modules to process particle data. The code transforms particle fields (masses and velocities) onto a uniform grid representation by leveraging an efficient k-d tree implementation, enabling fast neighbour searches for each cell. Inside the `tools` and `test1` directories the user can find examples for preparing particle data as input for the void finder. The particle input must have the following structure:

   * `int64: N, float32: ZETA` 
   * `float32(1:N): X` (Mpc)
   * `float32(1:N): Y` (Mpc)
   * `float32(1:N): Z` (Mpc)
   * `float32(1:N): V_x` (km/s)
   * `float32(1:N): V_y` (km/s)
   * `float32(1:N): V_z` (km/s)
   * `float32(1:N): M` ($M_{\odot}$)

   Being `N` the total number of particles, `ZETA` the redshift corresponding to the snapshot and `X, Y, Z`, `V_x, V_y, V_z`, `M` the position, velocity and mass of each particle. If points do not belong to the same time hypersurface (e.g a galaxy survey), `ZETA` can be understood as the mean redshift.

   In the particular case when particle velocities are unknown (for instance, in a galaxy survey), the user can specify this situation in `config/voids.dat` and linear theory will be applied to compute the velocity divergence field. In that case, the input must contain the following information:

   * `int64: N, float32: ZETA` 
   * `float32(1:N): X` (Mpc)
   * `float32(1:N): Y` (Mpc)
   * `float32(1:N): Z` (Mpc)
   * `float32(1:N): M` ($M_{\odot}$)

   When this type of input is selected, the code expects a `bin_file_partXXXXX` binary file. The user must specify an `iteration` or `snapshot` number `XXXXX`, as this allow the code to be run on several iterations without stopping. If this feature is not needed (for example analysing a galaxy survey), one can    simply provide **AVISM** with a bin_file_part00001 file and tell the code to find voids just in iteration `1`. Moreover, `bin_file_partXXXXX` files must be inside the `path_to_AVISM/input_data` directory.

   The python script `tools/uchuu2avism.py` serves as an example to properly prepare a particle input from a simulation output (in this case, a halo catalogue from [Mini-Uchuu](https://www.skiesanduniverses.org/Simulations/Uchuu/). Similarly, `tools/galaxy_survey.py` shows how to prepare a simple galaxy survey input (2MRS [John P. Huchra et al 2012 ApJS 199 26](https://iopscience.iop.org/article/10.1088/0067-0049/199/2/26) in that case), although we strongly recommend preprocessing galaxy surveys with external tools such as [CORAS](https://github.com/rlilow/CORAS), [Neural Networks](https://github.com/rlilow/2MRS-NeuralNet), or utilising constrained simulations of the Local Universe (e.g., see [Manticore-Local](https://arxiv.org/abs/2505.10682)) to obtain full reconstructions/descriptions of the density and velocity fields (non-linear in the last two examples), thus allowing to fully leverage the void finder capabilities. 

3. **Grid input:** `Option 2`
   
   If an AMR simulation snapshot is previosly processed and transformed into a uniform grid or if, for instance, the reconstruction of the density and velocity fields from a galaxy survey has been carried out by an external tool on a uniform grid, the user may need to apply **AVISM** directly on this data structure. In this case, the grid input must provide the following information:
   
   * `float32: ZETA` 
   * `float32(1:NX,1:NY,1:NZ): DELTA`
   * `float32(1:NX,1:NY,1:NZ): V_x` (km/s)
   * `float32(1:NX,1:NY,1:NZ): V_y` (km/s)
   * `float32(1:NX,1:NY,1:NZ): V_z` (km/s)

   With `ZETA` the redshift and `DELTA` and `(V_x, V_y, V_z)` the density contrast ($\delta = \rho / \rho_B - 1$) and velocity fields defined on the grid.

   When this type of input is selected, the code expects a `bin_file_gridXXXXX` binary file inside the `path_to_AVISM/input_data` directory. In analogy to the particle input, the user must specify an `iteration` or `snapshot` number `XXXXX` but, if this feature is not needed, one can simply provide **AVISM** with a bin_file_grid00001 and tell the code to find voids just in iteration `1`.

   The python scripts `tools/linear2M++_2_avism.py` and `tools/manticore2avism.py` serve as examples to properly prepare a grid input from a full linear reconstruction ([2M++_linear](https://cosmicflows.iap.fr/)) and a non-linear constrained simulation ([Manticore-Local](https://arxiv.org/abs/2505.10682))), respectively, of the Local Universe consisting of a uniform grid with the density and velocity fields defined.

4. **Arepo input:** `Option 3`

   A reader for [Arepo](https://arepo-code.org/) cosmological simulations (particularly the [IllustrisTNG](https://www.tng-project.org/) suite) is provided. This way, the user can give as input to the void finder all gas or dark matter particles from a snapshot. Care must be taken, however, as the [HDF5](https://www.hdfgroup.org/) library has to be properly installed and linked inside the Makefile. If this input is chosen, input files should be inside `path_to_AVISM/simu_arepo`.

   Note that, in order to use all dark matter or gas particles from Arepo's simulation, the internal AVISM reader (`Option 3`) has to be used. The TNG readers inside `tools` are just examples to handle [IllustrisTNG](https://www.tng-project.org/) halo catalogues to prepare a particle (`Option 2`) input.

## Run configuration

**AVISM** runs are highly customizable via `config/voids.dat`, where the user can specify different values for the physical parameters governing the void-finding procedure. Nevertheless, first, the user must introduce the right configuration for the code to run properly:

```
*******************************************************************************
*       General parameters block                                      
*******************************************************************************
files: first, last, every ---------------------------------------------------->
1300,1300,1
cells for the coarser grid --------------------------------------------------->
128,128,128
min and max levels of the hierarchies to be used by voidfind   --------------->
0,0
Hubble constant (h), omega matter -------------------------------------------->
0.678,0.31
comoving box size ------------------------------------------------------------>
147.5
periodic boundary conditions -> 0: no, 1: yes -------------------------------->
1
```

In this first block, the user must specify the iteration range, where `first` is the first iteration to be processed and `last` is the last one, in steps of `every`. Next, if voids-in-voids are required, $\ell_{min}$ and $\ell_{max}$ levels must be introduced, with the grid size $N_x, N_y, N_z$ corresponding to the coarser grid ($\ell_{min}$), gaining a factor of 2 in resolution as we go up in levels. If the user only needs a void population, a single level $\ell_{min} = \ell_{max}$ can be introduced, with a proper resolution. Last, $h$, $\Omega_m$ and $L$ are specified, with $L$ the size of the bounding box, and periodic boundary conditions can be activated.

```
*******************************************************************************
*       Void finder input parameters                                      
*******************************************************************************
maximum number of voids for dimensioning ------------------------------------->
30000000
density contrast threshold for centers --------------------------------------->
-0.8
density contrast threshold for the edges ------------------------------------->
10.
density gradient threshold --------------------------------------------------->
0.25
velocity divergence threshold ------------------------------------------------>
0.
min void radius (in Mpc) ----------------------------------------------------->
3.
min void radius to look for subvoids (in Mpc) -------------------------------->
3.
```

In this block, the physical thresholds for performing the void-finding algorithm are defined. The default values have been rigorously tested, but the user is free to change them. In descending order, and keeping the notation of _Monllor-Berbegal et al. 2025 accepted, A&A_ ([ArXiv](https://arxiv.org/abs/2509.25329)), we have $\delta_1$, $\delta_2$, $\nabla \delta_{th}$ and $\nabla \cdot \mathbf{v}_\text{th}$.

```
*******************************************************************************
*       Type of data to process                                   
*******************************************************************************
type of data -> 0:MASCLET, 1:BinPart, 2:BinGrid, 3:AREPO(snap) --------------->
0
IF MASCLET or Grid data: NX, NY, NZ (INPUT grid size) ------------------------>
256,256,256
IF AREPO data: files per snapshot, PartType (1:gas, 2:dm), DM mass (Msun) ---->
100,2,470000000.0
```

Here, the user must specify the input data format and, in the case of using `Option 0` or `2`, the input grid size should be supplied, with the previous $N_x, N_y, N_z$ values in `General parameters block` being ignored. On the other hand, if 'Option 3' is chosen (that is, Arepo data), the user must specify the number of files per snapshot, the particle type used as matter tracer and the dark matter particle mass if 'PartType=2'.

```
*******************************************************************************
*       Particle data handling parameters                                   
*******************************************************************************
maximum number of particles to read (DIMENSIONING) --------------------------->
2000000000
particle velocity is available -> 0: no, 1: yes ------------------------------>
1
interpolation scheme for particle velocity field (0:TSC, 1:SPH) -------------->
1
interpolation scheme for particle density field (0:TSC, 1:SPH) --------------->
1
Number of nearest neighbours to use to transform into grid ------------------->
32
```

Particle data can be handled in different ways, as the code allows for different velocity and density interpolations. We recommend, however, the SPH option, as it is well-optimized and tested, yielding better results than the TSC kernel. Here, the user must also specify if the particle velocity field is available and the number $N_\text{ngh}$ of nearest neighbours to use for the SPH interpolation. Depending on the grid size and the amount of particles provided as input, this quantity could be changed accordingly, although $N_\text{ngh} = 32$ is a proper value for a balanced case in which cells and particles are approximately in the same number.


## Output

**AVISM** provides two main output files: `voidsXXXXX` and `mapXXXXX`. 

### `voidsXXXXX`

`voidsXXXXX` corresponds to the void catalogue and contains information about the run (number of levels used for the grid, coarse grid size, box size, ...), information about each level (number of voids, filling fraction, ...), as well as each void main properties (ID, centre, volume, mass, ...). The structure is as follows:

----------------------------------------------------------------------------------

* $N_\ell$ /  $\ell_{min}$ /  $\ell_{max}$  /  $N_x^0$ /  $N_Y^0$ /  $N_z^0$  /  $L$
  
   - $\ell$ / $N_{cubes}$ / $N_{voids}$ / $N_{\ell-1}$ / FF / $\langle  \rho \rangle$
     
      - ID / $X$ / $Y$ / $Z$ / $X_G$ / $Y_G$ / $Z_G$ / Vol / $R$ / $\overline{\rho}$ / $\epsilon$ / IP / ID($\ell-1$) / $R(\ell-1)$ / Mass
        
        .
        .
        .
        
        (for all voids at this level)

      .
      .
      .
     
      (for all levels)

---------------------------------------------------------------------------------- 


Below, we provide three tables (one for each type of information given in `voidsXXXXX`) describing all variables listed before:

| Run variable  | Description |
| ------------- | ------------- |
| $N_\ell$  |  Number of grid levels |
| $\ell_{min}$ and $\ell_{max}$ | Minimum and maximum grid levels  |
| $N_x^0$, $N_y^0$, $N_z^0$ | Coarse (minimum) grid size in each cartesian direction|
| $L$  | Size (in Mpc) of the box in which the particles or grid are defined |


| Level variable  | Description |
| ------------- | ------------- |
| $\ell$  | Which level  |
| $N_{cubes}$ | Number of cubes found by the first void-finding step  |
| $N_{voids}$ | Final number of voids after merging and post-processing |
| $N_{\ell-1}$  | Number of voids in the previous level (parent voids) |
| FF  | Volume filling fraction of voids at this grid level |
| $\langle  \rho \rangle$ | Mean density used to define the density contrast



| Void property  | Description |
| ------------- | ------------- |
| ID | Void ID, corresponding to the ID of the biggest cube belonging to it |
| $X$, $Y$, $Z$| Void centre coordinates, defined as the centre of the biggest cube belonging to it |
| $X_G$, $Y_G$, $Z_G$| Void volume-weighed (geometrical) centre |
| Vol | Void total volume (in $\text{Mpc}^3$) |
| $R$ | Void effective radius (in Mpc)|
| $\overline{\rho}$ | Mean density (in matter background density units) inside the void |
| $\epsilon$ | Void ellipticity |
| IP | Void inverse porosity |
| ID($\ell-1$) | If applicable, ID of parent void at level $\ell-1$ |
| $R(\ell-1)$ | If applicable, radius of parent void at level $\ell-1$|
| Mass | Mass inside the void (in $M_{\odot}$)|



### `mapXXXXX`

The `mapXXXXX` file contains, for each level, the density and divergence velocity fields defined on the 3D grid used to identify voids, as well as the integer field indicating to which void each cell belongs. The basic structure for this file is as follows:

``` 
int32(1:NX,1:NY,1:NZ): MARCA
float32(1:NX,1:NY,1:NZ): DELTA
float32(1:NX,1:NY,1:NZ): DIVERGENCE

      .
      .
      .
     
 (for all levels, without spacing between blocks)
```

Below, we provide a brief description of the fields:

| Field | Description |
| ------------- | ------------- |
| MARCA  | Indicates to which void does a cell belong |
| DELTA | Density contrast field ($\delta$) utilised by the void-finding algorithm at this level |
| DIVERGENCE | Velocity divergence field ($\nabla \cdot \mathbf{v}$) utilised by the void-finding algorithm at this level |


### read_voids.py

Inside `tools` the user can find `read_voids.py`. This reader helps load **AVISM**'s output into Python for analysis or further post-processing. 
