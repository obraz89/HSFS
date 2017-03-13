
:: copy this file to your case dir and launch !!!

:: Change %HSFSRoot% to your HSFS project root folder !!!

@set HSFSRoot=E:/science/devel/HSFS

@set HSFSLib=%HSFSRoot%/lib
@set HSFSBin=%HSFSRoot%/bin_x64_Debug

@if %HSSTAB_CASE_DIR%.==. set HSSTAB_CASE_DIR=%CD% 

:: set paths to required libs

::@set WX_DIR=%HSFSRoot%/lib/wxMSW-2.8.11
::@set WX_ARCH=vc_x64_dll

@set WX_DIR=%HSFSRoot%\lib_vc14\wx-vc14-rlz-win64
@set WX_ARCH=vc140_x64_dll

@set MPICH2_DIR=%HSFSRoot%/lib_vc14/msmpi-win64

@set HDF5_DIR=E:\science\devel\HSFS\lib_vc14\hdf5_1.1
@set CGNS_DIR=E:\science\devel\HSFS\lib_vc14\cgns-vc14-rlz-win64

@set MKL_DIR=%HSFSRoot%/lib/mkl-xe-2015

@set PATH=%PATH%;%INTEL_COMPILERS%/bin/intel64;%WX_DIR%/lib/%WX_ARCH%/

@start %HSFSRoot%/StabSolverUnited.sln

===============================
