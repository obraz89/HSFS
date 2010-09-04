        !COMPILER-GENERATED INTERFACE MODULE: Sun Sep 05 01:26:30 2010
        MODULE SYSR_mod
          INTERFACE 
            SUBROUTINE SYSR(A,B,X,N)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(N,N)
              REAL(KIND=8) :: B(N)
              REAL(KIND=8) :: X(N)
            END SUBROUTINE SYSR
          END INTERFACE 
        END MODULE SYSR_mod
