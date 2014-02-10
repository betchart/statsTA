#!/bin/bash

source /uscmst1/prod/sw/cms/shrc prod
export SCRAM_ARCH=slc5_amd64_gcc472
cd /uscms/home/bbetchar/work/CMSSW_6_2_6/src && eval `scram runtime -sh` && cd - >& /dev/null
