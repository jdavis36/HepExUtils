#!/bin/bash
set -euo pipefail
cd /afs/cern.ch/work/l/lkang/CJLST/TreeTagger/
TREEFILENAME=$(sed $1!d alltrees.txt)
./batchTreeTagger.py -i $TREEFILENAME -o /eos/user/l/lkang/Active_Research/Discriminants/TaggedTrees/ -b branches.txt 2>&1
