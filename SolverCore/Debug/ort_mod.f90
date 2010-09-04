        !COMPILER-GENERATED INTERFACE MODULE: Sun Sep 05 01:26:30 2010
        MODULE ORT_mod
          INTERFACE 
            SUBROUTINE ORT(FI,M,L,EPS)
              INTEGER(KIND=4) :: L
              INTEGER(KIND=4) :: M
              COMPLEX(KIND=8) :: FI(M,L)
              REAL(KIND=8) :: EPS
            END SUBROUTINE ORT
          END INTERFACE 
        END MODULE ORT_mod
