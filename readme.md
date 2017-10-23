# FUSE

## Description and credits

This is a source code repository for the **Framework for Understanding Structural Errors** or **FUSE**. The initial implementation is described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735). This implementation involves four main additional features:

1. a snow module described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736),
2. a calibration mode using [Duan et al. (WRR, 1992)](http://dx.doi.org/10.1029/91WR02985) shuffled complex evolution alogrithm (SCE),
3. a distributed mode enabling to run FUSE on a grid, and
4. all the input, output and parameter files are now NetCDF files.

## FUSE modes and case studies

FUSE can be run in four complementary modes:

1. `run_def` runs FUSE with default parameter values,
2. `run_pre` runs FUSE with with a pre-defined set of parameter values,
3. `calib_sce` runs FUSE in an SCE-calibration mode,
4. `run_best` runs FUSE using the best (lowest RMSE) parameter set found by SCE.

To get you started with FUSE, we provide files for two case studies involing modeling at different scales:

* catchment scale: forcing and streamflow data for the [USGS 09066300 MIDDLE CREEK NEAR MINTURN, CO.](https://waterdata.usgs.gov/nwis/inventory/?site_no=09066300&agency_cd=USGS&) catchment - available [here](
https://www.dropbox.com/s/e8doyoju5fvsrvr/fuse_case_study_1.zip?dl=0) for download,  
* grid scale: forcing on a 1/8th degree grid for a small 10x10 domain - soon available for download.

Follow the following steps to run FUSE.

## A. Fork this repository and compile FUSE
1. Fork this repository to a directory `$(MASTER)` on your machine.
1. Change directory to `$(MASTER)/build/` and edit the `Makefile`, by:
   1. defining the name of the master directory (line 10),
   2. defining the fortran compiler (line 196),
   3. defining the path to the NetCDF libraries (lines 198-219, note that the NetCDF libraries must be compiled using the same compiler that you are using to run the program ).
 1. Compile the SCE code (see Section G below).
 1. Compile the FUSE code (type `make`).
 1.  Change to `$(MASTER)/bin/` and try running FUSE by typing `./fuse.exe`. If the output is `1st command-line argument is missing (DatString)`, you probably have compiled FUSE correctly. 
 
## B. Populate the bin directory
1. Move the file `fm_902_us.txt` provided for the catchment case study to `$(MASTER)/bin/`. The file is a `FILEMANAGER`. It defines the paths to the FUSE `settings`, `input`, `output` directories and the name of files essential to run FUSE (described  in Sections C and D). 
1. Edit the lines 4 to 6 of `fm_902_us.txt` and enter the path of the `fuse_case_study_1` directory on your machine.
   
## C. Populate the setup directory
The `setup` directory must contain the following files (provided for the catchment case study):

   1. The file `M_DECISIONS` (called `fuse_zDecisions_902.txt` in the case studies) describes the different options available in the FUSE modeling framework. These modeling decisions are described in detail by [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735), except decision 9 described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736).
   2. The file `CONSTRAINTS` (called `fuse_zConstraints_snow.txt` in the case studies) defines in particular the default parameter values and lower and upper parameter bounds. The list of parameters corresponds to those described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735) and [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736). 
   3. The file `MOD_NUMERIX` (called `fuse_zNumerix.txt` in the case studies) defines decisions regarding the numerical solution technique. Examples of the impact of these decisions are described by [Clark and Kavetski (WRR 2010)](http://dx.doi.org/10.1029/2009WR008894) and [Kavetski and Clark (WRR 2010)](http://dx.doi.org/10.1029/2009WR008896).
   4. The file `FORCINGINFO` (called `us_09066300_input_info.txt` in the catchment case study) provides metadata for the NetCDF input file. It defines the name of the input file and of the variables,  and also defines the start of the similation, the end of the warm-up period, and the end of the simulation.

## D. Populate the input directory
The `input` directory must contain the following files (provided for the catchment case study):

   1. The file `forcefile` (called `us_09066300_input.nc` in the catchment case study) contains the input data either in 2D (3D) arrays for modeling at the catchment (grid) scale. The name of this file is prescribed by `FORCINGINFO` (see C4).
   2. The file `BFILE` (called `us_09066300_elev_bands.nc` in the catchment case study) describe the elevation bands required when the snow module is on. The dimensions of this file must match that of `forcefile`.
   
Note that the dimension of the NetCDF files will determine if FUSE is run at the catchment or grid-scale. FUSE will look for the variables `lat` and `lon` and if they are arrays, it will run on the grid they define. This means that NetCDF input files for a single catchment must also include the variables `lat` and `lon`.
   
## E. Run the puppy

Run FUSE unsing default parameter values at the catchment scale:
```
./fuse.exe fm_902_us.txt us_09066300 902 run_def
```

then calibrate it:

```
./fuse.exe fm_902_us.txt us_09066300 902 calib_sce
```

then run it with the best SCE parameter set:

```
./fuse.exe fm_902_us.txt us_09066300 902 run_best
```

where
`$1` is the muster file,
`$2` is the ID of the basin,
`$3` is the ID of the FUSE model,
`$4` is the FUSE mode.

Note that:
* the string `$2`_ `$3` will be included in name of most of the input and output files,
* `$4` will define the name of `M_DECISIONS`.

## F. Content of the output directory
Running FUSE in its different modes will create the following files in the `output` directory (provided for the catchment case study for comparison purposes):
* the files whose name contains `runs` contain the simulations, 
* the files whose name contains `para` contain the parameter values,
* the last element of the file name indicates which FUSE mode was used.
   
## G. Compile SCE
The code of the shuffled complex evolution method (`$(MASTER)/build/FUSE_SRC/FUSE_SCE/sce.f`) was written in F77, so it must be compiled separately. If you use `ifort`, try the following flags:
  ```
  ifort -O2 -c -fixed sce_16plus.f
  ```

If necessary, rename the compiled file, so that it can be found by the `Makefile`, which by default will be looking for a file named `sce_16plus.o`.

### License
FUSE is distributed under the GNU Public License Version 3. For details see the file `LICENSE` in the FUSE root directory or visit the [online version](http://www.gnu.org/licenses/gpl-3.0.html).
