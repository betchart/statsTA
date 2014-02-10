#!/bin/bash

export SCRAM_ARCH=slc5_amd64_gcc472
source /vols/cms/grid/setup.sh
cd /home/hep/bbetchar/work/CMSSW_6_2_0_patch1/src && eval `scram runtime -sh` && cd - >& /dev/null
