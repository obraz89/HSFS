
:: @if not "%WindowsSdkDir%" == "" (
::	set "PATH=%WindowsSdkDir%bin;%PATH%"
::	set "INCLUDE=%WindowsSdkDir%include;%INCLUDE%"
::	set "LIB=%WindowsSdkDir%lib;%LIB%"
::)
@set HSFLOWStabBin=%CD%/bin
@set HSFSRoot=E:/science/devel/HSFS

@set WX_DIR=%HSFSRoot%\lib_vc14\wx-vc14-rlz-win64
@set WX_ARCH=vc140_x64_dll

@set HDF5_BIN=E:\lib\hdf5\bin

@set PATH=%PATH%;%INTEL_COMPILERS%/bin/intel64/;%WX_DIR%/lib/%WX_ARCH%;%HDF5_BIN%

::@if %HSSTAB_CASE_DIR%.==. set HSSTAB_CASE_DIR=%CD% 

::@start %HSFlowStabRoot%/StabSolverUnited.sln
::call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\vcvarsall.bat amd64"

::start "M= 4 wedge aoa=0 Re1=7.381x10^6 " /LOW cmd /K "%HSFLOWStabBin%\HSFlowStab-main.exe -l log.txt ."
start "M= 3.5 plate aoa=0" /LOW cmd /K "mpiexec -n 1 %HSFLOWStabBin%\HSFlowStab-main.exe -l log.txt ."
