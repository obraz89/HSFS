        !COMPILER-GENERATED INTERFACE MODULE: Sun Sep 05 01:26:28 2010
        MODULE SYS_mod
          INTERFACE 
            SUBROUTINE SYS(A,B,X,N)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(N,N)
              COMPLEX(KIND=8) :: B(N)
              COMPLEX(KIND=8) :: X(N)
            END SUBROUTINE SYS
          END INTERFACE 
        END MODULE SYS_mod
