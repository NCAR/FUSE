# FUSE

## Description and credits

This is a source code repository for the **Framework for Understanding Structural Errors** or **FUSE**. FUSE initial implementation is described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735).

This implementation involves four main additional features:

1. a snow module described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736),
2. a calibration mode using [Duan et al. (WRR, 1992)](http://dx.doi.org/10.1029/91WR02985) shuffled complex evolution alogrithm (SCE),
3. a distributed mode enabling to run FUSE on a grid, and
4. all the input, output and parameter files are now NetCDF files.

## FUSE modes

FUSE can be run in four complentary modes:

1. `run_def` runs FUSE with default parameter values,
2. `run_pre` runs FUSE with with a pre-defined set of parameter values,
3. `calib_sce` runs FUSE in a SCE-calibration mode,
4. `run_best` runs FUSE using the best (highest RMSE) parameter set found by SCE.

To get FUSE running, follow the following steps.

## A. Fork this repository and compile FUSE
1. Fork this repository to a directory `$(MASTER)` on your machine.
2. Change directory to `$(MASTER)/build/` and edit the `Makefile`, by:
   1. defining the name of the master directory (line 10),
   2. defining the fortran compiler (line 196),
   3. defining the path to the NetCDF libraries (lines 198-219, note that the NetCDF libraries must be compiled using the same compiler that you are using to run the program ).
 4. Compile the code (i.e. type `make` or `make -f Makefile`).
 
## B. Populate the setup directory
The setup directory must contain the following files:

   1. The file `M_DECISIONS` (can be called anything, and in the example is called `fuse_zDecisions_snow.txt`) describes the different options available in the FUSE modeling framework. Each of these modeling decisions is described in detail by [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735), except decision 9 described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736).
   2. The file `CONSTRAINTS` (can be called anything, and in the example is called `fuse_zConstraints_snow.txt`) defines the default parameter values and lower and upper parameter bounds. The list of parameters corresponds to those described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735) and [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736). There is a lot in this file, the important columns are the default parameter values and lower and upper parameter bounds.
   3. The file `MOD_NUMERIX` (can be called anything, and in the example is called `fuse_zNumerix.txt`) defines decisions regarding the numerical solution technique. Examples of the impact of these decisions are described by [Clark and Kavetski (WRR 2010)](http://dx.doi.org/10.1029/2009WR008894) and [Kavetski and Clark (WRR 2010)](http://dx.doi.org/10.1029/2009WR008896).
   4. The file `FORCINGINFO` (can be called anything, and in the example is called `us_05585000_input_info.txt`) provides metadata for the NetCDF input file. It defines the name of the input file and of the variables,  and also defines the start of the similation, the end of the warm-up period, and the end of the simulation.
   5. The file `FILEMANAGER` (can be called anything, and in the example is called `us_09210500_902_fuse_file_manager.txt`) defines the name of the files listed above and the directories in which FUSE settings, FUSE input, and FUSE output are kept.

## C. Populate the input directory
1. The setup directory must contain the following files:

   1. The file `M_DECISIONS` (can be called anything, and in the example is called `fuse_zDecisions_snow.txt`) describes the different options available in the FUSE modeling framework. Each of these modeling decisions is described in detail by [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735), except decision 9 described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736).
   2. The file `CONSTRAINTS` (can be called anything, and in the example is called `fuse_zConstraints_snow.txt`) defines the default parameter values and lower and upper parameter bounds. The list of parameters corresponds to those described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735) and [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736). There is a lot in this file, the important columns are the default parameter values and lower and upper parameter bounds.
   
## D. Run the puppy
Background: different driver programs were written to fulfill different objectives. In each driver program the options are defined through command-line arguments. This facilitates running multiple instances of the executable on a cluster. The filenames are quite long (and descriptive) to avoid different model runs writing to the same file.

Forcing and streamflow data for the catchment [USGS 08013000 Calcasieu River near Glenmora, LA](http://waterdata.usgs.gov/nwis/inventory/?site_no=08013000&agency_cd=USGS&amp;) are provided as an example in `$(MASTER)/input`. FUSE can be run using an uniform random sampling of feasible parameter space using the Sobol sequence. The NetCDF output in `$(MASTER)/output` was created using the following command line argument:
```
./fuse_URS.exe fuse_direktor_08013000.txt 08013000 070 2 0 1.e-2 1.e-2 1.0000000000 10
```
where
`$1` is the muster file,
`$2` if the ID of the basin,
`$3` is the ID of the FUSE model,
`$4` is the method used to temporally integrate model equations (2 is implicit Euler),
`$5` is a switch between fixed and adaptive sub-steps,
`$6` is the absolute tolerance for defining the length of sub-steps (adaptive sub-stepping),
`$7` is the relative tolerance for defining the length of sub-steps (adaptive sub-stepping),
`$8` is the maximum length of the time step (days) and
`$9` is the number of random samples desired.

These arguments override the information provided in the control files, specifically:
* `$2` is used to define the name of the `FORCINGINFO` file, and overwrites the information provided in `fuse_fileManager.txt`.
* `$3` is used to define the FUSE model used, and overwrites the information provided in `M_DECISIONS` (UNLESS the ID is negative, in which case the model decisions are read from the file. The list of model indices is defined in `$(MASTER)/settings/fuse_rModelList.txt`.
* `$4` through `$8` overwrites the information provided in the `MOD_NUMERIX` file.

## D. Content of the output directory
1. The setup directory must contain the following files:

   1. The file `M_DECISIONS` (can be called anything, and in the example is called `fuse_zDecisions_snow.txt`) describes the different options available in the FUSE modeling framework. Each of these modeling decisions is described in detail by [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735), except decision 9 described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736).
   2. The file `CONSTRAINTS` (can be called anything, and in the example is called `fuse_zConstraints_snow.txt`) defines the default parameter values and lower and upper parameter bounds. The list of parameters corresponds to those described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735) and [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736). There is a lot in this file, the important columns are the default parameter values and lower and upper parameter bounds.

## F. Compile SCE
1. The code of the shuffled complex evolution method (SCE, in file `$(MASTER)/build/FUSE_SRC/FUSE_SCE/sce.f`, [Duan et al., 1992](http://dx.doi.org/10.1029/91WR02985)) was written in F77, so it must be compiled separately. We compile it using `ifort` and the following flags:
  ```
  ifort -c -fixed -O3 -r8 sce.f  
  ```

2. If necessary, rename the compiled file, so that it can be found by `$(MASTER)/build/FUSE_SRC/FUSE_SCE/URS_driver_sce.f90`, which by default will be looking for a file named `sce.o`.

3. Adapt and compile `$(MASTER)/build/Makefile_sce` following the steps A1 to A3.

### License
FUSE is distributed under the GNU Public License Version 3. For details see the file `LICENSE` in the FUSE root directory or visit the [online version](http://www.gnu.org/licenses/gpl-3.0.html).
