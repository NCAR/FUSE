MODULE metaparams
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Describe all parameters used in the model (used to define NetCDF output files, etc.)
! ---------------------------------------------------------------------------------------
! variable definitions
USE nrtype
USE multibands
USE model_defn,ONLY:SMODL
USE model_defnames
IMPLICIT NONE
CHARACTER(LEN=11), DIMENSION(200)      :: PNAME       ! parameter names
CHARACTER(LEN=52), DIMENSION(200)      :: PDESC       ! parameter long names (description of variable)
CHARACTER(LEN= 8), DIMENSION(200)      :: PUNIT       ! parameter units
INTEGER(I4B)                           :: I           ! loop through parameter sets
INTEGER(I4B)                           :: IBAND       ! loop through bands
CHARACTER(LEN=2)                       :: TXT_IBAND   ! band index as a character
INTEGER(I4B)                           :: NOUTPAR     ! number of model parameters for output
CONTAINS
! ---------------------------------------------------------------------------------------
SUBROUTINE PARDESCRIBE()
I=0  ! initialize counter
! adjustable model parameters
I=I+1; PNAME(I)='RFERR_ADD  '; PDESC(I)='additive rainfall error                            '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='RFERR_MLT  '; PDESC(I)='multiplicative rainfall error                      '; PUNIT(I)='-       '
I=I+1; PNAME(I)='MAXWATR_1  '; PDESC(I)='maximum total storage in the upper layer           '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='MAXWATR_2  '; PDESC(I)='maximum total storage in the lower layer           '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='FRACTEN    '; PDESC(I)='fraction total storage as tension storage          '; PUNIT(I)='-       '
I=I+1; PNAME(I)='FRCHZNE    '; PDESC(I)='fraction tension storage in recharge zone          '; PUNIT(I)='-       '
I=I+1; PNAME(I)='FPRIMQB    '; PDESC(I)='fraction of baseflow in primary reservoir          '; PUNIT(I)='-       '
I=I+1; PNAME(I)='RTFRAC1    '; PDESC(I)='fraction of roots in the upper layer               '; PUNIT(I)='-       '
I=I+1; PNAME(I)='PERCRTE    '; PDESC(I)='percolation rate                                   '; PUNIT(I)='mm day-1'
I=I+1; PNAME(I)='PERCEXP    '; PDESC(I)='percolation exponent                               '; PUNIT(I)='-       '
I=I+1; PNAME(I)='SACPMLT    '; PDESC(I)='percolation multiplier in the SAC model            '; PUNIT(I)='-       '
I=I+1; PNAME(I)='SACPEXP    '; PDESC(I)='percolation exponent in the SAC model              '; PUNIT(I)='-       '
I=I+1; PNAME(I)='PERCFRAC   '; PDESC(I)='fraction of percolation to tension storage         '; PUNIT(I)='-       '
I=I+1; PNAME(I)='FRACLOWZ   '; PDESC(I)='fraction of soil excess to lower zone              '; PUNIT(I)='-       '
I=I+1; PNAME(I)='IFLWRTE    '; PDESC(I)='interflow rate                                     '; PUNIT(I)='mm day-1'
I=I+1; PNAME(I)='BASERTE    '; PDESC(I)='baseflow rate                                      '; PUNIT(I)='mm day-1'
I=I+1; PNAME(I)='QB_POWR    '; PDESC(I)='baseflow exponent                                  '; PUNIT(I)='-       '
I=I+1; PNAME(I)='QB_PRMS    '; PDESC(I)='baseflow depletion rate                            '; PUNIT(I)='-       '
I=I+1; PNAME(I)='QBRATE_2A  '; PDESC(I)='baseflow depletion rate for primary reservoir      '; PUNIT(I)='day-1   '
I=I+1; PNAME(I)='QBRATE_2B  '; PDESC(I)='baseflow depletion rate for secondary reservoir    '; PUNIT(I)='day-1   '
I=I+1; PNAME(I)='SAREAMAX   '; PDESC(I)='maximum saturated area                             '; PUNIT(I)='-       '
I=I+1; PNAME(I)='AXV_BEXP   '; PDESC(I)='ARNO/VIC b exponent                                '; PUNIT(I)='-       '
I=I+1; PNAME(I)='LOGLAMB    '; PDESC(I)='mean value of the log-transformed topographic index'; PUNIT(I)='log m   '
I=I+1; PNAME(I)='TISHAPE    '; PDESC(I)='shape parameter for the topo index Gamma distribtn '; PUNIT(I)='-       '
I=I+1; PNAME(I)='TIMEDELAY  '; PDESC(I)='time delay in runoff (routing)                     '; PUNIT(I)='day     '
I=I+1; PNAME(I)='MBASE      '; PDESC(I)='snow model base melt temperature                   '; PUNIT(I)='deg.C   '
I=I+1; PNAME(I)='MFMAX      '; PDESC(I)='snow model maximum melt factor                     '; PUNIT(I)='mm/(C-d)'
I=I+1; PNAME(I)='MFMIN      '; PDESC(I)='snow model minimum melt factor                     '; PUNIT(I)='mm/(C-d)'
I=I+1; PNAME(I)='PXTEMP     '; PDESC(I)='rain-snow partition temperature                    '; PUNIT(I)='deg.C   '
I=I+1; PNAME(I)='OPG        '; PDESC(I)='maximum relative precip difference across the bands'; PUNIT(I)='-       '
I=I+1; PNAME(I)='LAPSE      '; PDESC(I)='maximum temperature difference across the bands    '; PUNIT(I)='deg.C   '
! derived model parameters
I=I+1; PNAME(I)='MAXTENS_1  '; PDESC(I)='maximum tension storage in the upper layer         '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='MAXTENS_1A '; PDESC(I)='maximum storage in the recharge zone               '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='MAXTENS_1B '; PDESC(I)='maximum storage in the lower zone                  '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='MAXFREE_1  '; PDESC(I)='maximum free storage in the upper layer            '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='MAXTENS_2  '; PDESC(I)='maximum tension storage in the lower layer         '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='MAXFREE_2  '; PDESC(I)='maximum free storage in the lower layer            '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='MAXFREE_2A '; PDESC(I)='maximum storage in the primary baseflow reservoir  '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='MAXFREE_2B '; PDESC(I)='maximum storage in the secondary baseflow reservoir'; PUNIT(I)='mm      '
I=I+1; PNAME(I)='RTFRAC2    '; PDESC(I)='fraction of roots in the lower layer               '; PUNIT(I)='-       '
I=I+1; PNAME(I)='QBSAT      '; PDESC(I)='baseflow at saturation (derived parameter)         '; PUNIT(I)='mm day-1'
I=I+1; PNAME(I)='POWLAMB    '; PDESC(I)='mean value of power-transformed topographic index  '; PUNIT(I)='m**(1/n)'
I=I+1; PNAME(I)='MAXPOW     '; PDESC(I)='max value of power-transformed topographic index   '; PUNIT(I)='m**(1/n)'
! model bands parameters
IF(SMODL%iSNOWM.EQ.iopt_temp_index) THEN !loop through snow model bands
 I=I+1; PNAME(I)='N_BANDS    '; PDESC(I)='number of basin bands in model                     '; PUNIT(I)='=       '
 I=I+1; PNAME(I)='Z_FORCING  '; PDESC(I)='elevation of model forcing data                    '; PUNIT(I)='m       '
 DO IBAND=1,N_BANDS
  WRITE(TXT_IBAND,'(I2)') IBAND              ! convert band no. to text
  IF (IBAND.LT.10) TXT_IBAND(1:1) = '0'      ! pad with zeros 
  I=I+1; PNAME(I)='Z_MID'//TXT_IBAND//'    '; PDESC(I)='basin band mid-point elevation                     '; PUNIT(I)='m       '
  I=I+1; PNAME(I)='AF'//TXT_IBAND//'       '; PDESC(I)='basin band area fraction                           '; PUNIT(I)='-       '
 END DO
ENDIF
! numerical solution parameters
I=I+1; PNAME(I)='SOLUTION   '; PDESC(I)='0=explicit euler; 1=implicit euler                 '; PUNIT(I)='-       '
I=I+1; PNAME(I)='TIMSTEP_TYP'; PDESC(I)='0=fixed time steps; 1=adaptive time steps          '; PUNIT(I)='-       '
I=I+1; PNAME(I)='INITL_GUESS'; PDESC(I)='0=old state; 1=explicit half-step; 2=expl full-step'; PUNIT(I)='-       '
I=I+1; PNAME(I)='JAC_RECOMPT'; PDESC(I)='0=variable; 1=constant sub-step; 2=const full step '; PUNIT(I)='-       '
I=I+1; PNAME(I)='CK_OVRSHOOT'; PDESC(I)='0=always take full newton step; 1=line search      '; PUNIT(I)='-       '
I=I+1; PNAME(I)='SMALL_ESTEP'; PDESC(I)='0=step truncation; 1=look-ahead; 2=step absorption '; PUNIT(I)='-       '
I=I+1; PNAME(I)='ERRTRUNCABS'; PDESC(I)='absolute temporal truncation error tolerance       '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='ERRTRUNCREL'; PDESC(I)='relative temporal truncation error tolerance       '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='ERRITERFUNC'; PDESC(I)='iteration convergence tolerance for function values'; PUNIT(I)='mm      '
I=I+1; PNAME(I)='ERR_ITER_DX'; PDESC(I)='iteration convergence tolerance for dx             '; PUNIT(I)='-       '
I=I+1; PNAME(I)='THRESH_FRZE'; PDESC(I)='threshold for freezing the Jacobian                '; PUNIT(I)='mm      '
I=I+1; PNAME(I)='FSTATE_MIN '; PDESC(I)='fractional minimum value of state                  '; PUNIT(I)='-       '
I=I+1; PNAME(I)='STEP_SAFETY'; PDESC(I)='safety factor in step-size equation                '; PUNIT(I)='-       '
I=I+1; PNAME(I)='RMIN       '; PDESC(I)='minimum step size multiplier                       '; PUNIT(I)='-       '
I=I+1; PNAME(I)='RMAX       '; PDESC(I)='maximum step size multiplier                       '; PUNIT(I)='-       '
I=I+1; PNAME(I)='NITER_TOTAL'; PDESC(I)='maximum number of iterations in the implicit scheme'; PUNIT(I)='-       '
I=I+1; PNAME(I)='MIN_TSTEP  '; PDESC(I)='minimum time step length                           '; PUNIT(I)='day     '
I=I+1; PNAME(I)='MAX_TSTEP  '; PDESC(I)='maximum time step length                           '; PUNIT(I)='day     '
! parameter identifier
I=I+1; PNAME(I)='SOBOL_INDX '; PDESC(I)='indentifier for Sobol parameter set                '; PUNIT(I)='-       '
NOUTPAR=I
END SUBROUTINE PARDESCRIBE
END MODULE metaparams
