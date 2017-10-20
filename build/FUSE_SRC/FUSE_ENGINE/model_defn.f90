! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark
! Modified by Brian Henn to include snow model, 6/2013
! ---------------------------------------------------------------------------------------
MODULE model_defn
 USE nrtype
 ! FUSE version
 character(*),parameter::FUSE_version="FUSE 1.0"
 logical,parameter::FUSE_enabled=.true.
 ! list of combinations in each model component
 INTEGER, PARAMETER :: NDEC = 9                           ! number of model decisions
 TYPE DESC
  CHARACTER(LEN=16)                    :: MCOMPONENT      ! description of model component
 END TYPE DESC
 TYPE(DESC), DIMENSION(2)              :: LIST_RFERR      ! rainfall error
 TYPE(DESC), DIMENSION(3)              :: LIST_ARCH1      ! upper-layer architecture
 TYPE(DESC), DIMENSION(4)              :: LIST_ARCH2      ! lower-layer architecture
 TYPE(DESC), DIMENSION(3)              :: LIST_QSURF      ! surface runoff
 TYPE(DESC), DIMENSION(3)              :: LIST_QPERC      ! percolation
 TYPE(DESC), DIMENSION(2)              :: LIST_ESOIL      ! evaporation
 TYPE(DESC), DIMENSION(2)              :: LIST_QINTF      ! interflow
 TYPE(DESC), DIMENSION(2)              :: LIST_Q_TDH      ! time delay in runoff
 TYPE(DESC), DIMENSION(2)              :: LIST_SNOWM      ! snow model
 ! structure that holds (x) unique combinations
 TYPE UMODEL
  INTEGER(I4B)                         :: MODIX           ! model index
  CHARACTER(LEN=256)                   :: MNAME           ! model name
!  CHARACTER(LEN=16)                    :: RFERR           ! rainfall error
  INTEGER(I4B)                         :: iRFERR
!  CHARACTER(LEN=16)                    :: ARCH1           ! upper-layer architecture
  INTEGER(I4B)                         :: iARCH1
!  CHARACTER(LEN=16)                    :: ARCH2           ! lower-layer architecture
  INTEGER(I4B)                         :: iARCH2
!  CHARACTER(LEN=16)                    :: QSURF           ! surface runoff
  INTEGER(I4B)                         :: iQSURF
!  CHARACTER(LEN=16)                    :: QPERC           ! percolation
  INTEGER(I4B)                         :: iQPERC
!  CHARACTER(LEN=16)                    :: ESOIL           ! evaporation
  INTEGER(I4B)                         :: iESOIL
!  CHARACTER(LEN=16)                    :: QINTF           ! interflow
  INTEGER(I4B)                         :: iQINTF
!  CHARACTER(LEN=16)                    :: Q_TDH           ! time delay in runoff
  INTEGER(I4B)                         :: iQ_TDH
  INTEGER(I4B)                         :: iSNOWM           ! snow
   END TYPE UMODEL
 ! structure to hold model state names
 TYPE SNAMES
!  CHARACTER(LEN=8)                     :: SNAME           ! state name
  INTEGER(I4B)                         :: iSNAME          ! integer value of state name
 END TYPE SNAMES
 ! structure to hold model flux names
 TYPE FNAMES
  CHARACTER(LEN=16)                    :: FNAME           ! state name
 END TYPE FNAMES
! max steps in routing function
  INTEGER(I4B),PARAMETER::NTDH_MAX=500
! model definitions
 CHARACTER(LEN=256)                    :: FNAME_NETCDF_RUNS    ! NETCDF output filename for model runs
 CHARACTER(LEN=256)                    :: FNAME_NETCDF_PARA    ! NETCDF output filename for model parameters
 CHARACTER(LEN=256)                    :: FNAME_NETCDF_PARA_SCE   ! NETCDF output filename for model parameters produced by SCE
 CHARACTER(LEN=256)                    :: FNAME_NETCDF_PARA_PRE   ! NETCDF filename for pre-defined model parameters set
 CHARACTER(LEN=256)                    :: FNAME_PREFIX    ! prefix for desired output files
 CHARACTER(LEN=256)                    :: FNAME_TEMPRY    ! prefix for temporary output files
 CHARACTER(LEN=256)                    :: FNAME_ASCII     ! ASCII output filename
 TYPE(UMODEL),DIMENSION(5000)          :: AMODL           ! (model definition -- all)
 TYPE(UMODEL)                          :: SMODL           ! (model definition -- single model)
 TYPE(SNAMES),DIMENSION(7)             :: CSTATE          ! (list of model states for SMODL)
 INTEGER(I4B)                          :: NSTATE=0        ! number of model states
 TYPE(FNAMES),DIMENSION(50)            :: C_FLUX          ! (list of model fluxes for SMODL)
 INTEGER(I4B)                          :: N_FLUX=0        ! number of model fluxes
 ! --------------------------------------------------------------------------------------
END MODULE model_defn
