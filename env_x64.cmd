
:: @if not "%WindowsSdkDir%" == "" (
::	set "PATH=%WindowsSdkDir%bin;%PATH%"
::	set "INCLUDE=%WindowsSdkDir%include;%INCLUDE%"
::	set "LIB=%WindowsSdkDir%lib;%LIB%"
::)

@set PATH=%PATH%;%INTEL_COMPILERS%/bin/intel64/;%WX_DIR%/lib/%WX_ARCH%/
@start StabSolverUnited.sln

