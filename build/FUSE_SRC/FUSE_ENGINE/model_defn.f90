MODULE model_defn
 USE nrtype
 ! FUSE version
 character(*),parameter::FUSE_version="FUSE 1.0"
 logical,parameter::FUSE_enabled=.true.
 ! list of combinations in each model component
 TYPE DESC
  CHARACTER(LEN=10)                    :: MCOMPONENT      ! description of model compopnent
 END TYPE DESC
 TYPE(DESC), DIMENSION(2)              :: LIST_RFERR      ! rainfall error
 TYPE(DESC), DIMENSION(3)              :: LIST_ARCH1      ! upper-layer architecture
 TYPE(DESC), DIMENSION(4)              :: LIST_ARCH2      ! lower-layer architecture
 TYPE(DESC), DIMENSION(3)              :: LIST_QSURF      ! surface runoff
 TYPE(DESC), DIMENSION(3)              :: LIST_QPERC      ! percolation
 TYPE(DESC), DIMENSION(2)              :: LIST_ESOIL      ! evaporation
 TYPE(DESC), DIMENSION(2)              :: LIST_QINTF      ! interflow
 TYPE(DESC), DIMENSION(2)              :: LIST_Q_TDH      ! time delay in runoff
 ! structure that holds (x) unique combinations
 TYPE UMODEL
  INTEGER(I4B)                         :: MODIX           ! model index
  CHARACTER(LEN=256)                   :: MNAME           ! model name
!  CHARACTER(LEN=10)                    :: RFERR           ! rainfall error
  INTEGER(I4B)                         :: iRFERR
!  CHARACTER(LEN=10)                    :: ARCH1           ! upper-layer architecture
  INTEGER(I4B)                         :: iARCH1  
!  CHARACTER(LEN=10)                    :: ARCH2           ! lower-layer architecture
  INTEGER(I4B)                         :: iARCH2  
!  CHARACTER(LEN=10)                    :: QSURF           ! surface runoff
  INTEGER(I4B)                         :: iQSURF  
!  CHARACTER(LEN=10)                    :: QPERC           ! percolation
  INTEGER(I4B)                         :: iQPERC  
!  CHARACTER(LEN=10)                    :: ESOIL           ! evaporation
  INTEGER(I4B)                         :: iESOIL  
!  CHARACTER(LEN=10)                    :: QINTF           ! interflow
  INTEGER(I4B)                         :: iQINTF  
!  CHARACTER(LEN=10)                    :: Q_TDH           ! time delay in runoff
  INTEGER(I4B)                         :: iQ_TDH
   END TYPE UMODEL
 ! structure to hold model state names
 TYPE SNAMES
!  CHARACTER(LEN=6)                     :: SNAME           ! state name
  INTEGER(I4B)                         :: iSNAME          ! integer value of state name
 END TYPE SNAMES
 ! structure to hold model flux names
 TYPE FNAMES
  CHARACTER(LEN=11)                    :: FNAME           ! state name
 END TYPE FNAMES
! max steps in routing function
  INTEGER(I4B),PARAMETER::NTDH_MAX=500
! model definitions
 CHARACTER(LEN=256)                    :: FNAME_NETCDF    ! NETCDF output filename
 CHARACTER(LEN=256)                    :: FNAME_PREFIX    ! prefix for desired output files
 CHARACTER(LEN=256)                    :: FNAME_TEMPRY    ! prefix for temporary output files
 CHARACTER(LEN=256)                    :: FNAME_ASCII     ! ASCII output filename
 INTEGER(I4B),PARAMETER                :: OUTFILE_UNIT=21 ! unit for output file
 TYPE(UMODEL),DIMENSION(5000)          :: AMODL           ! (model definition -- all)
 TYPE(UMODEL)                          :: SMODL           ! (model definition -- single model)
 TYPE(SNAMES),DIMENSION(6)             :: CSTATE          ! (list of model states for SMODL)
 INTEGER(I4B)                          :: NSTATE=0        ! number of model states 
 TYPE(FNAMES),DIMENSION(50)            :: C_FLUX          ! (list of model fluxes for SMODL)
 INTEGER(I4B)                          :: N_FLUX=0        ! number of model fluxes 
 ! --------------------------------------------------------------------------------------
END MODULE model_defn
