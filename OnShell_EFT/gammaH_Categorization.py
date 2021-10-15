#!/usr/bin/env python

import sys, getopt
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

    datatype = {'vector<float>':'vector<float>','vector<short>':'vector<short>','vector<long>':'vector<long>','vector<bool>':'vector<bool>','Char_t':'B','UChar_t':'b','Bool_t':'O','Short_t':'S','UShort_t':'s','Counter':'I','Int_t':'I','UInt_t':'i','Double_t':'D','Double32_t':'d','Float_t':'F','Float16_t':'f','Long_t':'G','ULong_t':'g','Long64_t':'L','ULong64_t':'l','CharStar':'C'}

    if os.path.exists(tagtreefilename):
        print("ERROR: " + tagtreefilename + " already exists!\n")

    else:
        print("Pre-existing output TTree not found --- safe to proceed")
        branchlist = []
        branchdict = {}

#         with open(branchlistpath) as file:
#             while (line := file.readline().rstrip()):
#                 branchlist.append(line)
                

        with open(branchlistpath) as f:
            blist = [line.rstrip() for line in f]

        for branch in blist:
            if branch: branchlist.append(branch)
        
#         branchlist = ["xsec", "GenHMass"]
        
        branchdict["EventTag"] = [] #np.array([])
        branchdict["Bin40"] = [] #np.array([])
        for branch in branchlist:
            branchdict[branch] = [] #np.array([])
        
	#======================  load splines ========================#

	spline_list = init_spline()
	
        for tree in ["candTree", "candTree_failed"]:
        #for tree in ["candTree"]:
            print("\n================ Reading events from " + tree + " ================\n")

            f = ROOT.TFile(filename)
            t = f.Get("ZZTree/"+tree)
            t.Print()
            treebranches = [ x.GetName() for x in t.GetListOfBranches() ]
            print t.GetEntries()
	    for i, entry in enumerate(t, start=1):
		    if i % 1000 == 0:
			print (i)
		 #=========================== Fill failed events with dummy and skip to loop over branches ===========================

		    branchdict["Bin40"].append(f.Get("ZZTree/Counters").GetBinContent(40)) #= np.append(branchdict["Bin40"], f.Get("ZZTree/Counters").GetBinContent(40))
                    if tree == "candTree_failed":
                        branchdict["EventTag"].append(-999) #= np.append(branchdict["EventTag"], -999)
                        break
		#=========================== Tagging event by category ===========================

                    M = t.ZZMass
    #               M = t.GenHMass
                    tag =  Tag(t.nExtraLep, t.nExtraZ, t.nCleanedJetsPt30, t.nCleanedJetsPt30BTagged_bTagSF, t.JetQGLikelihood, t.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, t.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal, t.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, t.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal, t.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal, t.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, t.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, t.p_HadWH_mavjj_JECNominal, t.p_HadWH_mavjj_true_JECNominal, t.p_HadZH_mavjj_JECNominal, t.p_HadZH_mavjj_true_JECNominal, t.JetPhi, t.ZZMass, t.ZZPt, t.PFMET, t.PhotonIsCutBasedLooseID,t.PhotonPt, False, False, spline_list)
		#=========================== Saving category tag ===========================
                    
                    branchdict["EventTag"].append(tag) 
		
                #=========================== Loop over branches to copy values ===========================
            
                    for branchcand in branchlist:
                        if branchcand in treebranches:
                            exec("global datatypecode; datatypecode = datatype[t.GetLeaf('"+branchcand+"').GetTypeName()]")
                            exec("global branchval; branchval = t." + branchcand)
                            branchdict[branchcand].append(branchval) # = np.append(branchdict[branchcand], branchval)
                        else:
                            branchdict[branchcand].append(-999) # = np.append(branchdict[branchcand], -999)
                    	    break

        
	#================= Make output directory ====================
	out_path = '/'.join(tagtreefilename.split("/")[:-1])
	os.makedirs(out_path)
        #Path("/".join(tagtreefilename.split("/")[:-1])).mkdir( 0o755, True, True )
        ftagtree = ROOT.TFile(tagtreefilename, "CREATE")
        tagtree = ROOT.TTree("eventTree", "eventTree")
        
        print("\n================ Create new TTree and fill with branches ================\n")
        
        for key in branchdict.keys():
            print(len(branchdict[key]), key)
            exec("value"+key.replace('-', 'm')+" = array('"+datatypecode.replace("F","f")+"', [0])")
            exec("tagtree.Branch('"+key.replace('-', 'm')+"', value"+key.replace('-', 'm')+", '"+key.replace('-', 'm')+"/"+datatypecode+"')")
        
        print("\n================ Fill with branches in event loop ================\n")
        
        for elem in trange(len(branchdict[key])):
            #for key in branchdict.keys():
	        #exec("value"+key.replace('-', 'm') +"[0] = branchdict[key][elem]")
            tagtree.Fill()

        print("\n================ Check new TTree branches ================\n")
            
        tagtree.Print()
        tagtree.Write("", ROOT.TObject.kOverwrite)
        ftagtree.Close()

        print("\n", Counter(branchdict["EventTag"]))

if __name__ == "__main__":
    main(sys.argv[1:])
