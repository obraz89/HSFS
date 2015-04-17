#!/bin/bash
#
# Make all projects
#

PROJ='
../shared
../components/SmallMat
../components/PhysCommon
../components/Profile
../components/MF.CGNS2D
../components/MF.CGNS3D
../components/PF.GlobSearch_Lapack
../components/PF.LocSearch
../components/WPTrack
../HSFLowStab-main
'

for p in $PROJ
do
#   echo $p
    pushd $p
    make || exit 1
    popd
done
