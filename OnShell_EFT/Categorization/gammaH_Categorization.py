#!/usr/bin/env python

import sys, getopt
sys.path.insert(0, '/afs/cern.ch/work/j/jejeffre/public/HEP_Ex_Tools/HepExUtils/OnShell_Utils/')
import os
import glob
import ROOT
from math import sqrt
import time
from os import path
#from pathlib import Path
import re
from tqdm import trange, tqdm
import numpy as np
from array import *
import random
from collections import Counter
from decimal import *
from OnShell_Category import *
#from rootpy.tree import Tree, BoolArrayCol, ShortArrayCol
#from rootpy.io import root_open
import root_numpy
from root_numpy import array2tree, tree2array
# Load Root Macro#

def main(argv):
    inputfile = ''
    outputfile = ''
    branchfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:b:",["ifile=","ofile=","bfile="])
    except getopt.GetoptError:
        print('batchTreeTagger.py -i <inputfile> -o <outputfile> -b <branchfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('batchTreeTagger.py -i <inputfile> -o <outputfile> -b <branchfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-b", "--bfile"):
            branchfile = arg
    
    if not all([inputfile, outputfile, branchfile]):
        print('batchTreeTagger.py -i <inputfile> -o <outputfile> -b <branchfile>')
        sys.exit(2)

    print("\n================ Reading user input ================\n")

    print('Input CJLST TTree is', inputfile)
    print('Output directory is', outputfile)
    print('Branch list file is', branchfile)

    print("\n================ Processing user input ================\n")


    filename = inputfile #treelist[numfile]
    branchlistpath = branchfile
    tagtreepath = outputfile

    ind = filename.split("/").index("200205_CutBased")

    tagtreefile = "/".join(filename.split("/")[ind:])
    tagtreefilename = tagtreepath+"/"+tagtreefile

    print(filename, "\n")
    print(tagtreefilename, "\n")

    if os.path.exists(tagtreefilename):
        print("ERROR: " + tagtreefilename + " already exists!\n")
    #elif os.path.exists("temp/"+'/'.join(tagtreefilename.split("/")[:-1])):
    #	print("ERROR: " + "temp/"+'/'.join(tagtreefilename.split("/")[:-1]) + " already exists!\n")
    else:
        print("Pre-existing output TTree not found --- safe to proceed")
        branchlist = []
        branchdict = {}
        datatypedict = {}
        
        with open(branchlistpath) as f:
            blist = [line.rstrip() for line in f]

        for branch in blist:
            if branch: branchlist.append(branch)
        
        
        branchdict["EventTag"] = [] #np.array([])
        branchdict["Bin40"] = [] #np.array([])
        for branch in branchlist:
            branchdict[branch] = [] #np.array([])
            datatypedict[branch] = []
	#======================  load splines ========================#

	spline_list = init_spline()
	
        #===================== load output File ======================
	path = '/'.join(tagtreefilename.split("/")[:-1])
        os.makedirs(path)
        #tagtree = ROOT.TTree("eventTree", "eventTree")
        #for tree in ["candTree", "candTree_failed"]:
        for tree in ["candTree"]:
            print("\n================ Reading events from " + tree + " ================\n")

            f = ROOT.TFile(filename)
            t = f.Get("ZZTree/"+tree)
            t.Print()
            treebranches = [ x.GetName() for x in t.GetListOfBranches() ]
            
            print(t.GetEntries())
	    
	    #### create branches to be calculated ####
	    eventTag=[]
	    Bin40=[]
            for i, entry in enumerate(t):
        	#=========================== Tagging event by category ===========================
	          if i % 1000 == 0:
        	     print (i)

	          M = t.ZZMass
		  Bin40.append(f.Get("ZZTree/Counters").GetBinContent(40))
		  #Checks to make sure values wont crash tagger if they do crash we will return -999 on the tagger output array#

		  # Signal VBF variables that will fail
		  if t.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal * t.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal * t.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal * t.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal * t.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal * t.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal == 0:
		    eventTag.append(-999)
                  # Bkg VBF variables that will fail 
		  elif t.pConst_JJVBF_BKG_MCFM_JECNominal * t.pConst_HadZH_BKG_MCFM_JECNominal * t.pConst_HadWH_BKG_MCFM_JECNominal * t.pConst_JJQCD_BKG_MCFM_JECNominal * t.pConst_JJVBF_BKG_MCFM_JECNominal * t.pConst_HadZH_BKG_MCFM_JECNominal * t.pConst_HadWH_BKG_MCFM_JECNominal * t.pConst_JJQCD_BKG_MCFM_JECNominal == 0:
		    eventTag.append(-999)
		  # QCD Scale variables that will fail #
		  elif t.p_HadZH_mavjj_true_JECNominal * t.p_HadWH_mavjj_true_JECNominal == 0:
		    eventTag.append(-999)
		  # VBF 1 Jet and VBF 2 jet variables that will fail #
		  elif t.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal * t.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal == 0:
		    eventTag.append(-999)
		  # WH and Vh variables that will fail #
		  elif t.p_HadWH_mavjj_JECNominal * t.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal * t.p_HadZH_mavjj_JECNominal * t.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal == 0:
		    eventTag.append(-999)
		  else:
        	    tag =  Tag(t.nExtraLep, t.nExtraZ, t.nCleanedJetsPt30, t.nCleanedJetsPt30BTagged_bTagSF, t.JetQGLikelihood, t.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, t.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, t.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, t.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, t.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, t.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, t.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, t.p_HadWH_mavjj_JECNominal, t.p_HadWH_mavjj_true_JECNominal, t.p_HadZH_mavjj_JECNominal, t.p_HadZH_mavjj_true_JECNominal, t.JetPhi, t.ZZMass, t.ZZPt, t.PFMET, t.PhotonIsCutBasedLooseID, t.PhotonPt, False, False, spline_list)
	            eventTag.append(tag)
	#========================== Cropping the tree to remove branches you dont want ====================	
	for branchname in treebranches:
          if branchname not in branchlist:
	    t.SetBranchStatus(branchname,0)
	ftagtree = ROOT.TFile.Open(tagtreefilename, "CREATE")
	newtree = t.CloneTree()
	#========================= Load up the branches you want to add =====================
	tag=np.array(eventTag,dtype=[('tag',np.int32)])
	Bin40=np.array(Bin40,dtype=[('Bin40',np.float)])
        array2tree(tag,tree=newtree)
	array2tree(Bin40,tree=newtree)
	newtree.Write()   
	# clean up TTags #
	#ftagtree.Delete("*;1")
	#ftagtree.Delete("*;2")
	print(ftagtree.ls())
	

if __name__ == "__main__":
    main(sys.argv[1:])
