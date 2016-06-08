        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun  7 22:07:59 2016
        MODULE LNSRCH__genmod
          INTERFACE 
            SUBROUTINE LNSRCH(XOLD,FOLD,G,P,X,F,STPMAX,CHECK,FUNC)
              REAL(KIND=8), INTENT(IN) :: XOLD(:)
              REAL(KIND=8), INTENT(IN) :: FOLD
              REAL(KIND=8), INTENT(IN) :: G(:)
              REAL(KIND=8), INTENT(INOUT) :: P(:)
              REAL(KIND=8), INTENT(OUT) :: X(:)
              REAL(KIND=8), INTENT(OUT) :: F
              REAL(KIND=8), INTENT(IN) :: STPMAX
              LOGICAL(KIND=4), INTENT(OUT) :: CHECK
              INTERFACE 
                FUNCTION FUNC(X)
                  REAL(KIND=8), INTENT(IN) :: X(:)
                  REAL(KIND=8) :: FUNC
                END FUNCTION FUNC
              END INTERFACE 
            END SUBROUTINE LNSRCH
          END INTERFACE 
        END MODULE LNSRCH__genmod
