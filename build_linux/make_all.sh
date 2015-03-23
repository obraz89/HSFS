#!/bin/bash
#
# Make all projects
#

PROJ='
../shared
'

for p in $PROJ
do
#   echo $p
    pushd $p
    make || exit 1
    popd
done
