
:: copy this file to your case dir and launch !!!

:: Change %HSFSRoot% to your HSFS project root folder !!!

@set HSFSRoot=E:/science/devel/HSFS

@set HSFSLib=%HSFSRoot%/lib
@set HSFSBin=%HSFSRoot%/bin_x64_Debug

@if %HSSTAB_CASE_DIR%.==. set HSSTAB_CASE_DIR=%CD% 

:: set paths to required libs

@set WX_DIR=%HSFSRoot%/lib/wxMSW-2.8.11
@set WX_ARCH=vc_x64_dll

@set MPICH2_DIR=%HSFSRoot%/lib/MPICH2

@set CGNS_DIR=%HSFSRoot%/lib/cgns-3.2/vc-rlz-win64

@set MKL_DIR=%HSFSRoot%/lib/mkl-xe-2015

@set PATH=%PATH%;%INTEL_COMPILERS%/bin/intel64;%WX_DIR%/lib/%WX_ARCH%/

@start %HSFSRoot%/StabSolverUnited.sln

===============================
