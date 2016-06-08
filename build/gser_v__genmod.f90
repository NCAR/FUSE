        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun  7 22:07:59 2016
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
