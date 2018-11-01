MODULE metaoutput
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Describe all variables used in the model (used to define NetCDF output files, etc.)
! ---------------------------------------------------------------------------------------
! variable definitions
USE nrtype
USE multibands,ONLY:N_BANDS
USE model_defn,ONLY:SMODL
USE model_defnames
IMPLICIT NONE
CHARACTER(LEN=11), DIMENSION(200)      :: VNAME       ! variable names
CHARACTER(LEN=52), DIMENSION(200)      :: LNAME       ! variable long names (descrition of variable)
CHARACTER(LEN=13), DIMENSION(200)      :: VUNIT       ! variable units
INTEGER(I4B)                           :: I           ! loop through variables
INTEGER(I4B)                           :: NOUTVAR     ! number of output variables
INTEGER(I4B)                           :: ISNW        ! loop through SWE states
CHARACTER(LEN=2)                       :: TXT_ISNW    ! band index as a character
CONTAINS
! ---------------------------------------------------------------------------------------
SUBROUTINE VARDESCRIBE()
I=0  ! initialize counter
! model forcing
I=I+1; VNAME(I)='ppt        '; LNAME(I)='precipitation rate                                 '; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='pet        '; LNAME(I)='potential evapotranspiration rate                  '; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='temp       '; LNAME(I)='mean air temperature                               '; VUNIT(I)='deg.C        '
I=I+1; VNAME(I)='obsq       '; LNAME(I)='observed runoff                                    '; VUNIT(I)='mm timestep-1'
! model states
I=I+1; VNAME(I)='tens_1     '; LNAME(I)='tension storage in the upper layer                 '; VUNIT(I)='mm           '
I=I+1; VNAME(I)='tens_1a    '; LNAME(I)='tension storage in the soil excess zone            '; VUNIT(I)='mm           '
I=I+1; VNAME(I)='tens_1b    '; LNAME(I)='tension storage in the soil recharge zone          '; VUNIT(I)='mm           '
I=I+1; VNAME(I)='free_1     '; LNAME(I)='free storage in the upper layer                    '; VUNIT(I)='mm           '
I=I+1; VNAME(I)='watr_1     '; LNAME(I)='total storage in the upper layer                   '; VUNIT(I)='mm           '
I=I+1; VNAME(I)='tens_2     '; LNAME(I)='tension storage in the lower layer                 '; VUNIT(I)='mm           '
I=I+1; VNAME(I)='free_2     '; LNAME(I)='free storage in the lower layer                    '; VUNIT(I)='mm           '
I=I+1; VNAME(I)='free_2a    '; LNAME(I)='free storage in the primary baseflow reservoir     '; VUNIT(I)='mm           '
I=I+1; VNAME(I)='free_2b    '; LNAME(I)='free storage in the secondary baseflow reservoir   '; VUNIT(I)='mm           '
I=I+1; VNAME(I)='watr_2     '; LNAME(I)='total storage in the lower layer                   '; VUNIT(I)='mm           '

IF(SMODL%iSNOWM.EQ.iopt_temp_index) THEN !loop through snow model bands

 print *, 'Creating variables for the snow model for ', N_BANDS ,'elevation bands'

 I=I+1; VNAME(I)='swe_tot    '; LNAME(I)='total storage as snow                              '; VUNIT(I)='mm           '

 DO ISNW=1,N_BANDS ! output each for each snow model band
  WRITE(TXT_ISNW,'(I2)') ISNW              ! convert band no. to text
  IF (ISNW.LT.10) TXT_ISNW(1:1) = '0'      ! pad with zeros
  I=I+1; VNAME(I)='swe_z'//TXT_ISNW//' '! first create SWE band series
  LNAME(I)='elevation band snow water equivalent               '; VUNIT(I)='mm           '
  I=I+1; VNAME(I)='snwacml_z'//TXT_ISNW  ! then the accumulation series
  LNAME(I)='new band snowpack accumulation, in water equivalent'; VUNIT(I)='mm timestep-1'
  I=I+1; VNAME(I)='snwmelt_z'//TXT_ISNW  ! then the melt series
  LNAME(I)='band snowpack melt, in water equivalent            '; VUNIT(I)='mm timestep-1'
 END DO

ENDIF

! model fluxes
I=I+1; VNAME(I)='eff_ppt    '; LNAME(I)='effective precipitation rate                       '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='satarea    '; LNAME(I)='saturated area                                     '; VUNIT(I)='-            '
I=I+1; VNAME(I)='qsurf      '; LNAME(I)='surface runoff                                     '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='evap_1a    '; LNAME(I)='evaporation from soil excess zone                  '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='evap_1b    '; LNAME(I)='evaporation from soil recharge zone                '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='evap_1     '; LNAME(I)='evaporation from the upper soil layer              '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='evap_2     '; LNAME(I)='evaporation from the lower soil layer              '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='rchr2excs  '; LNAME(I)='flow from recharge zone to excess zone             '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='tens2free_1'; LNAME(I)='flow from tension to free storage, lower layer     '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='oflow_1    '; LNAME(I)='bucket overflow from upper soil layer              '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='tens2free_2'; LNAME(I)='flow from tension to free storage, lower layer     '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='qintf_1    '; LNAME(I)='interflow                                          '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='qperc_12   '; LNAME(I)='percolation from upper to lower soil layer         '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='qbase_2    '; LNAME(I)='baseflow                                           '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='qbase_2a   '; LNAME(I)='baseflow from primary baseflow reservoir           '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='qbase_2b   '; LNAME(I)='baseflow from secondary baseflow reservoir         '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='oflow_2    '; LNAME(I)='bucket overflow from lower soil layer              '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='oflow_2a   '; LNAME(I)='bucket overflow from primary baseflow reservoir    '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='oflow_2b   '; LNAME(I)='bucket overflow from secondary baseflow reservoir  '; VUNIT(I)='mm timestep-1'
! errors in model states (due to excessive extrapolation)
I=I+1; VNAME(I)='err_tens_1 '; LNAME(I)='excessive extrapolation: upper tension storage     '; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='err_tens_1a'; LNAME(I)='excessive extrapolation: upper excs tension storage'; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='err_tens_1b'; LNAME(I)='excessive extrapolation: upper rech tension storage'; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='err_free_1 '; LNAME(I)='excessive extrapolation: upper free storage        '; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='err_watr_1 '; LNAME(I)='excessive extrapolation: upper total storage       '; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='err_tens_2 '; LNAME(I)='excessive extrapolation: lower tension storage     '; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='err_free_2 '; LNAME(I)='excessive extrapolation: lower free storage        '; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='err_free_2a'; LNAME(I)='excessive extrapolation: 1st baseflow reservoir    '; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='err_free_2b'; LNAME(I)='excessive extrapolation: 2nd baseflow reservoir    '; VUNIT(I)='mm day-1     '
I=I+1; VNAME(I)='err_watr_2 '; LNAME(I)='excessive extrapolation: lower total storage       '; VUNIT(I)='mm day-1     '
! time check
I=I+1; VNAME(I)='chk_time   '; LNAME(I)='length of time step included in weighted average   '; VUNIT(I)='days         '
! model numerix
I=I+1; VNAME(I)='num_funcs  '; LNAME(I)='number of function calls                           '; VUNIT(I)='-            '
I=I+1; VNAME(I)='numjacobian'; LNAME(I)='number of times the Jacobian is calculated         '; VUNIT(I)='-            '
I=I+1; VNAME(I)='sub_accept' ; LNAME(I)='number of sub-steps accepted (taken)               '; VUNIT(I)='-            '
I=I+1; VNAME(I)='sub_reject' ; LNAME(I)='number of sub-steps tried but rejected             '; VUNIT(I)='-            '
I=I+1; VNAME(I)='sub_noconv' ; LNAME(I)='number of sub-steps tried that did not converge    '; VUNIT(I)='-            '
I=I+1; VNAME(I)='max_iterns' ; LNAME(I)='maximum number of iterations in implicit euler     '; VUNIT(I)='-            '
! model runoff (for BATEA, assumed to be last)
I=I+1; VNAME(I)='q_instnt   '; LNAME(I)='instantaneous runoff                               '; VUNIT(I)='mm timestep-1'
I=I+1; VNAME(I)='q_routed   '; LNAME(I)='routed runoff                                      '; VUNIT(I)='mm timestep-1'

print *, 'Setting NOUTVAR (number of forcing, state and flux variables) to', I
NOUTVAR=I

END SUBROUTINE VARDESCRIBE
END MODULE metaoutput
