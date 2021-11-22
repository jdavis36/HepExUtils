#!/bin/bash
set -euo pipefail
cd /afs/cern.ch/work/j/jejeffre/public/Higgs/CMSSW_10_2_5/src 
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/j/jejeffre/public/HEP_Ex_Tools/HepExUtils/OnShell_EFT/Categorization
TREEFILENAME=$(sed $1!d alltrees.txt)
./gammaH_Categorization.py -i $TREEFILENAME -o Test_Tag_gammaH -b branchlist.txt 2>&1
