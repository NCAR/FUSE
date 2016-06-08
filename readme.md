# FUSE 

## Description and credits

This is a source code repository for the Framework for Understanding Structural Errors or FUSE. FUSE initial implementation is described in the following paper:

 * Clark, M. P., Slater, A. G., Rupp, D. E., Woods, R. A., Vrugt, J. A., Gupta, H. V., Wagener, T. and Hay, L. E.: Framework for Understanding Structural Errors (FUSE): A modular framework to diagnose differences between hydrological models, Water Resour. Res., 44(12), [doi:10.1029/2007WR006735](http://dx.doi.org/10.1029/2007WR006735), 2008.

If you use FUSE, please credit this publication.

##A. Compile
1. Change directory to `$(MASTER)/build/`, where `$(MASTER)` is the directory where you un-tar'd the package.
1. Edit the Makefile, by
  1. Defining the name of the master directory (line 12);
  1. Defining the path to the NetCDF libraries (lines 162-163); and
  1. Define the fortran compiler - note that the NetCDF libraries must be compiled using the same compiler that you are using to run the program.
1. Compile the code (i.e., type `make` or `make -f Makefile`)

##B. Define the files to be used
1. In the directory $(MASTER)/bin, edit the file `fuse_MusterDirektor.txt` to point to the file that defines all of the files required to run FUSE (i.e., change the path to wherever you want to put the file manager file).
1. Now, in the file defined in B1 above (which by default is `$(MASTER)/settings/fuse_fileManager.txt`), modify the path to define the directories where you would like to keep FUSE settings, FUSE input, and FUSE output. In the example given these directories are simply
  * `$(MASTER)/settings/`
  * `$(MASTER)/input/`
  * `$(MASTER)/output/`
but it may be necessary to store input and output data on a different disk partition that is not backed up (hence, the flexibility). Note that there is also flexibility to change the name of the control files.

##C. Assemble control and input files
1. The file `M_DECISIONS` (can be called anything, and by default is called `fuse_zDecisions.txt`) describes the different options available in the FUSE modeling framework. Each of these modeling decisions is described in detail by Clark et al. (WRR, 2008).
2. The file `CONSTRAINTS` (can be called anything, and by default is called `fuse_zConstraints.txt`) defines the default parameter values and lower and upper parameter bounds. The list of parameters corresponds to those described in Clark et al. (WRR, 2008). There is alot in this file, mostly used for hierarchacal Bayesian modeling, but the important columns are the default parameter values and lower and upper parameter bounds (everything else is used for research that is still underway).
3. The file `MOD_NUMERIX` (can be called anything, and by default is called `fuse_zNumerix.txt`) defines decisions regarding the numerical solution technique. Examples of the impact of these decisions are described by Clark and Kavetski (WRR 2010) and Kavetski and Clark (WRR 2010).
4. The file `FORCINGINFO` (can be called anything, and by default is called `forcinginfo.txt`) provides metadata for the model input files. It defines the name of the data file, the number of columns in the data file, the column numbers for precip, potential ET, and runoff, the number of header lines, and the row numbers for the start of the similation, the end of the warm-up period, and the end of the simulation.
55. The model input file resides in the location defined by `fuse_fileManager.txt` and has the name defined in `FORCINGINFO`. The only real restriction is that it is an ASCII file. The time step of the forcing data (in days) is defined on the second line of the forcing data file.

##D. Run the puppy
Background: My modis operandi is to write different driver programs to fulfill different objectives. The example provided is for the case of uniform random sampling through the feasible parameter space using the Sobol sequence.

In each driver program that I write I define the options used through command-line arguments. This facilitates running multiple instances of the executable on a cluster. The filenames are quite long (and descriptive) to avoid different model runs writing to the same file.

The NetCDF output in $(MASTER)/output was created using the following command line argument:
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
`$8` is the maximum length of the time step (days), 
`$9` is the number of random samples desired.

These arguments override the information provided in the control files, specifically: 
`$2` is used to define the name of the "FORCINGINFO" file, and overwrites the information
   provided in fuse_fileManager.txt.
`$3` is used to define the FUSE model used, and overwrites the information provided in `M_DECISIONS` (UNLESS the ID is negative, in which case the model decisions are read from the file. The list of model indices is defined in `$(MASTER)/settings/fuse_rModelList.txt`. 
`$4` through `$8` overwrites the information provided in the `MOD_NUMERIX` file
`$9` simply defines the number of random samples desired.

##E. Plot the results
Plot the content of the input and output files, for instance using the script r_scripts/plot_fuse_input_output.R. This will make basic consistency tests, e.g. check that the length of the input and output time series and the indices in the settings file are consistent. Since parameter values were obtained from a uniformrandom sampling, do not expect a good fit of the observed discharge at this stage.

##F. Calibrate FUSE using SCE
1. The code of the shuffled complex evolution method (file sce.f, Duan et al., 1992, 1993) was written in F77. It must be compiled separately. Compile it with ifort or whichever F77 compiler you like, e.g. using: 

```
ifort -O2 -c -fixed sce.f
```

2. If necessary, rename the compiled file, so that it can be found by `FUSE_SCE/URS_driver_sce.f90`, which by default will be looking for a file nammed `sce.o`.

3. Adapt and compile `Makefile_sce` following the steps A1 to A3.

4. Calibrate FUSE by adapating and executing one following command lines - the second will take 
a while to run:

```
./fuse_URS_sce.exe fuse_direktor_08013000.txt 08013000 070 2 0 1.e-2 1.e-2 1.0000000000 3 50 3 1.e-3
./fuse_URS_sce.exe fuse_direktor_08013000.txt 08013000 070 2 0 1.e-2 1.e-2 1.0000000000 3 5000 3 1.e-3
```

where `$1` to $8 are as above, `$9` is the number of calibrated parameter sets desired, `$10` is the maximum number of trials before optimization is terminated, `$11` is the number of shuffling loops the objective function must change by `PCENTO` (max=9), `$12` is the percentage `PCENTO` (1 is 1%).

5. Note that the objective function is RMSE, defined in the file `FUSE_SCE/fuse_rsme.f90` and called by the wrapper `FUSE_SCE/functn.f90`.

6. Note that in `FUSE_SCE/URS_driver_sce.f90` the following line turns off the production 
of time series outputs to save space: 

```
OUTPUT_FLAG = .FALSE.    ! .TRUE. if desire time series output
```

However, the value of the objective function of each model run is backed up together with parameter values, so that the parameter set associated with the best results can then be run (see Section G below).

##G. Run FUSE with calibrated parameter values
1. Adapt and compile `Makefile_sce_merge` following the steps A1 to A3.

2. Create a file containing the name of the netCDF produced by the SCE calibration (e.g. `sce_params_08013000_FUSE_070__NMETH-2_SSTEP-0_2_0_1.e-2_1.e-2_1.0000000000_SCE_5000_3_1.e-3.nc`)

3. Run FUSE for the parameter set leading to the lowest RMSE by adapting and executing the following command line:

```
./fuse_URS_sce_merge.exe fuse_direktor_08013000.txt 08013000 sce_ncfiles_08013000.txt
```

where
`$1` is the muster file as above, `$2` is the ID of the basin run as above and `$3` is the name of the file containing the name of the netcdf file produced by SCE

##H. Plot the calibrated model runs
This will check that the RMSE values returned by SCE are consistent with those computed in R using the observed discharge from the input file and the simulated discharge from the output file.

##I. Final remarks
Martyn has a number of other driver programs written to produce the results presented in Clark and Kavetski (WRR 2010) and Kavetski and Clark (WRR 2010). Have a look through the papers and if you see something there that you think would be useful for your project, then just let him know and he'll provide some guidance on how to get it done.

### License

FUSE is distributed under the GNU Public License Version 3. For details see the file `COPYING` in the FUSE root directory or visit the [online version](http://www.gnu.org/licenses/gpl-3.0.html).

