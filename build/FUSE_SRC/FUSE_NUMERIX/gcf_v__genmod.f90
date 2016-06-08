        !COMPILER-GENERATED INTERFACE MODULE: Fri Oct 16 17:25:20 2015
        MODULE GCF_V__genmod
          INTERFACE 
            FUNCTION GCF_V(A,X,GLN)
              REAL(KIND=8), INTENT(IN) :: A(:)
              REAL(KIND=8), INTENT(IN) :: X(:)
              REAL(KIND=8) ,OPTIONAL, INTENT(OUT) :: GLN(:)
              REAL(KIND=8) :: GCF_V(SIZE(A))
            END FUNCTION GCF_V
          END INTERFACE 
        END MODULE GCF_V__genmod
