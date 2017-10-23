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

* Catchment scale: forcing and streamflow data for the [USGS 08013000 Calcasieu River near Glenmora, LA](http://waterdata.usgs.gov/nwis/inventory/?site_no=08013000&agency_cd=USGS&amp;) catchment available [here](
https://www.dropbox.com/s/ht4hqegcvu60x2m/fuse_catch_ex.zip?dl=0)  
* Grid scale: forcing on a 1/8th degree grid for a small 10x10 domain in Colorado available [here](
https://www.dropbox.com/s/c3a23549rp57sen/fuse_grid_ex.zip?dl=0)  . 

Follow the following steps to run FUSE.

## A. Fork this repository and compile FUSE
1. Fork this repository to a directory `$(MASTER)` on your machine.
1. Change directory to `$(MASTER)/build/` and edit the `Makefile`, by:
   1. defining the name of the master directory (line 10),
   2. defining the fortran compiler (line 196),
   3. defining the path to the NetCDF libraries (lines 198-219, note that the NetCDF libraries must be compiled using the same compiler that you are using to run the program ).
 1. Compile SCE code (see Section F below)
 1. Compile FUSE code (type `make` or `make -f Makefile`).
 1.  Change to `$(MASTER)/bin/` and try running FUSE by typing `./fuse_dist`. If the output is `1st command-line argument is missing (DatString)`, you probably have compiled FUSE correctly. 

## B. Populate the setup directory
The `setup` directory must contain the following files (provided for the two case studies):

   1. The file `M_DECISIONS` (called `fuse_zDecisions_snow.txt` in the case studies) describes the different options available in the FUSE modeling framework. These modeling decisions are described in detail by [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735), except decision 9 described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736).
   2. The file `CONSTRAINTS` (called `fuse_zConstraints_snow.txt` in the case studies) defines the default parameter values and lower and upper parameter bounds. The list of parameters corresponds to those described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735) and [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736). There is a lot in this file, the important columns are the default parameter values and lower and upper parameter bounds.
   3. The file `MOD_NUMERIX` (called `fuse_zNumerix.txt` in the case studies) defines decisions regarding the numerical solution technique. Examples of the impact of these decisions are described by [Clark and Kavetski (WRR 2010)](http://dx.doi.org/10.1029/2009WR008894) and [Kavetski and Clark (WRR 2010)](http://dx.doi.org/10.1029/2009WR008896).
   4. The file `FORCINGINFO` (called `us_05585000_input_info.txt` in the case studies) provides metadata for the NetCDF input file. It defines the name of the input file and of the variables,  and also defines the start of the similation, the end of the warm-up period, and the end of the simulation.
   5. The file `FILEMANAGER` (called `fm_901_us.txt` in the case studies) defines the name of the files listed above and the directories in which FUSE settings, FUSE input, and FUSE output are kept.

## C. Populate the input directory
The `input` directory must contain the following files (provided for the two case studies):

   1. The file `NAMEHERE` (can be called anything, and in the example is called `us_03237500_elev_bands.nc`) contains the input data either in 2D (3D) arrays for modeling at the catchment (grid) scale.
   2. The file `NAMEHERE` (can be called anything, and in the example is called `us_06037500_elev_bands.nc`) describe the elevation bands required when the snow module is on as 1D (2D) arrays for modeling at the catchment (grid) scale.
   
Note that the dimension of the NetCDF files will determine if FUSE is run at the catchment or grid-scale. FUSE will look for the variables `lat` and `lon` and if they are arrays, it will run on the grid they define. This means that NetCDF input files for a single catchment must also include the variables `lat` and `lon` 
   
   
## D. Run the puppy

Run FUSE unsing default parameter values at the catchment scale:
```
./fuse_snow_dist_catch.exe fm_901_us.txt us_11264500 run_def
```

or at the grid scale:

```
./fuse_snow_dist_catch.exe fm_901_us.txt us_11264500 run_def
```

where
`$1` is the muster file,
`$2` is the ID of the basin,
`$3` is the FUSE mode.

## E. Content of the output directory
Running FUSE in its different modes will create the following files in the `output` directory (provided for the two case studies for comparison purposes):
   * `run_def`: the file `NAMEHERE` (called `us_06906800_694_para.nc` and `us_06784000_694_runs.nc`) contains model parameters.
   * `run_pre`: the file `NAMEHERE` (called `us_06906800_694_para.nc` and `us_06784000_694_runs.nc`) contains model parameters.
   * `calib_sce`: the file `NAMEHERE` (called `us_10336660_902_para_best.nc`) contains model parameters determined by SCE.
   * `run_sce`: the file `NAMEHERE` (called `us_06746095_902_runs_best.nc`) contains model parameters determined by SCE.
   
## F. Compile SCE
The code of the shuffled complex evolution method (`$(MASTER)/build/FUSE_SRC/FUSE_SCE/sce.f`) was written in F77, so it must be compiled separately. We compile it using `ifort` and the following flags:
  ```
  ifort -c -fixed -O3 -r8 sce.f  
  ```

If necessary, rename the compiled file, so that it can be found by `$(MASTER)/build/FUSE_SRC/FUSE_SCE/URS_driver_sce.f90`, which by default will be looking for a file named `sce.o`.

### License
FUSE is distributed under the GNU Public License Version 3. For details see the file `LICENSE` in the FUSE root directory or visit the [online version](http://www.gnu.org/licenses/gpl-3.0.html).
