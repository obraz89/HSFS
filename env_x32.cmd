
:: @if not "%WindowsSdkDir%" == "" (
::	set "PATH=%WindowsSdkDir%bin;%PATH%"
::	set "INCLUDE=%WindowsSdkDir%include;%INCLUDE%"
::	set "LIB=%WindowsSdkDir%lib;%LIB%"
::)

@set PATH=%PATH%;%INTEL_COMPILERS%/bin/ia32/
@start StabSolverUnited.sln
