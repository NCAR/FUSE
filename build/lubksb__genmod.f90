        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  1 11:18:02 2016
        MODULE LUBKSB__genmod
          INTERFACE 
            SUBROUTINE LUBKSB(A,INDX,B)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              INTEGER(KIND=4), INTENT(IN) :: INDX(:)
              REAL(KIND=8), INTENT(INOUT) :: B(:)
            END SUBROUTINE LUBKSB
          END INTERFACE 
        END MODULE LUBKSB__genmod
