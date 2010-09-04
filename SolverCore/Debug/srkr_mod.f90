        !COMPILER-GENERATED INTERFACE MODULE: Sun Sep 05 01:26:28 2010
        MODULE SRKR_mod
          INTERFACE 
            SUBROUTINE SRKR(F,Y,H,M,OP)
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: F(M)
              REAL(KIND=8) :: Y
              REAL(KIND=8) :: H
              EXTERNAL OP
            END SUBROUTINE SRKR
          END INTERFACE 
        END MODULE SRKR_mod
