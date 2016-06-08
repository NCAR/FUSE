        !COMPILER-GENERATED INTERFACE MODULE: Fri Oct 16 17:25:20 2015
        MODULE GSER_V__genmod
          INTERFACE 
            FUNCTION GSER_V(A,X,GLN)
              REAL(KIND=8), INTENT(IN) :: A(:)
              REAL(KIND=8), INTENT(IN) :: X(:)
              REAL(KIND=8) ,OPTIONAL, INTENT(OUT) :: GLN(:)
              REAL(KIND=8) :: GSER_V(SIZE(A))
            END FUNCTION GSER_V
          END INTERFACE 
        END MODULE GSER_V__genmod
