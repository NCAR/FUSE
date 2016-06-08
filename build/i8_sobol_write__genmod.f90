        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun  7 22:07:53 2016
        MODULE I8_SOBOL_WRITE__genmod
          INTERFACE 
            SUBROUTINE I8_SOBOL_WRITE(M,N,SKIP,R,FILE_OUT_NAME)
              INTEGER(KIND=8) :: N
              INTEGER(KIND=8) :: M
              INTEGER(KIND=8) :: SKIP
              REAL(KIND=8) :: R(M,N)
              CHARACTER(*) :: FILE_OUT_NAME
            END SUBROUTINE I8_SOBOL_WRITE
          END INTERFACE 
        END MODULE I8_SOBOL_WRITE__genmod
