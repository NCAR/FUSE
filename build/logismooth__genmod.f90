        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun  7 22:07:59 2016
        MODULE LOGISMOOTH__genmod
          INTERFACE 
            PURE FUNCTION LOGISMOOTH(STATE,STATE_MAX,PSMOOTH)
              REAL(KIND=8), INTENT(IN) :: STATE
              REAL(KIND=8), INTENT(IN) :: STATE_MAX
              REAL(KIND=8), INTENT(IN) :: PSMOOTH
              REAL(KIND=8) :: LOGISMOOTH
            END FUNCTION LOGISMOOTH
          END INTERFACE 
        END MODULE LOGISMOOTH__genmod