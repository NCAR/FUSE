        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep  1 11:18:02 2016
        MODULE LUDCMP__genmod
          INTERFACE 
            SUBROUTINE LUDCMP(A,INDX,D)
              REAL(KIND=8), INTENT(INOUT) :: A(:,:)
              INTEGER(KIND=4), INTENT(OUT) :: INDX(:)
              REAL(KIND=8), INTENT(OUT) :: D
            END SUBROUTINE LUDCMP
          END INTERFACE 
        END MODULE LUDCMP__genmod
