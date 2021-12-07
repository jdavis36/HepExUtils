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
from root_numpy import array2tree, tree2array
sys.path.insert(0, '/afs/cern.ch/work/j/jejeffre/public/HEP_Ex_Tools/HepExUtils/OnShell_Utils/')
from Calc_Weight import *
import numpy as np

### 

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print('CountCategories.py -i <inputdir> -o <outputfile> ')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('CountCategories.py -i <inputdir> -o <outputfile> ')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
  	elif opt in ("-o", "--ofile"):
            outputfile = arg    
    if not all([inputfile,outputfile]):
        print('CountCategories.py -i <inputfile> -o <outputfile> ')
        sys.exit(2)

    print("\n================ Reading user input ================\n")

    print('Input Directory is', inputfile)
    
    print("\n================ Processing user input ================\n")


    filename = inputfile #treelist[numfile]
        
    print(filename, "\n")
    
    # List the tags definitions #	
    Untagged      = 0
    VBF1jTagged   = 1
    VBF2jTagged   = 2
    VHLeptTagged  = 3
    VHHadrTagged  = 4
    ttHLeptTagged = 5
    ttHHadrTagged = 6
    VHMETTagged   = 7
    Boosted       = 8
    gammaHTagged  = 9
    

    print("\n================ Reading events from " + filename + " ================\n")
    
    # Open the output file #

    out = open(outputfile, 'a')

    # Load the array of total count of events #

    count=np.array([0,0,0,0,0,0,0,0,0,0])
    
    # Load luminosity variable #

    lumi_2016=35.9
    lumi_2017=41.5
    lumi_2018=59.7
    out.write("Untagged	   VBF1jTagged	  VBF2jTagged	VHLeptTagged	VHHadrTagged	ttHLeptTagged	ttHHadrTagged	VHMETTagged	Boosted		gammaG"+'\n')
    
    # Recursively Loop through directories
    for root, subdirs, files in os.walk(filename):
      #print('--\nroot = ' + root)
      list_file_path = os.path.join(root, 'ZZ4lAnalysis.root')
      if os.path.isfile(list_file_path):
	nUntagged      = 0.
        nVBF1jTagged   = 0.
        nVBF2jTagged   = 0.
        nVHLeptTagged  = 0.
        nVHHadrTagged  = 0.
        nttHLeptTagged = 0.
        nttHHadrTagged = 0.
        nVHMETTagged   = 0.
        nBoosted       = 0.
	ngammaHTagged  = 0.
        print('list_file_path = ' + list_file_path)
        f = ROOT.TFile.Open(list_file_path)
        t = f.Get('candTree')
	t.Print()
        tag = tree2array(tree=t,branches=["tag"])
	tag = tag.astype(int)
        print(tag)
	# Calculate Weights according to custom #
	#eventweight = Calc_Event_Weight_200205(t,list_file_path)
        eventweight = Calc_Event_Weight_2021_gammaH(t,list_file_path)

	# ======== Choose Luminosity per year ==========#
	lumi = 0 
	if "2016" in list_file_path:
	  lumi = lumi_2016
	elif "2017" in list_file_path:
	  lumi = lumi_2017
	elif "2018" in list_file_path:
	  lumi = lumi_2018
	# Sum up the event weights for each category #
	for i in range(len(tag)):
	  if tag[i] == Untagged:
	    nUntagged = nUntagged + eventweight[i] * lumi #/ weightsum
	  if tag[i] == VBF1jTagged:
	    nVBF1jTagged = nVBF1jTagged + eventweight[i] * lumi #/ weightsum 
	  if tag[i] == VBF2jTagged:
	    nVBF2jTagged = nVBF2jTagged + eventweight[i]  * lumi #/ weightsum
	  if tag[i] == VHLeptTagged:
	    nVHLeptTagged = nVHLeptTagged + eventweight[i] * lumi #/ weightsum
	  if tag[i] == VHHadrTagged:
	    nVHHadrTagged = nVHHadrTagged + eventweight[i] * lumi #/ weightsum
	  if tag[i] == ttHLeptTagged:
	    nttHLeptTagged = nttHLeptTagged + eventweight[i] * lumi #/ weightsum
	  if tag[i] == ttHHadrTagged:
	    nttHHadrTagged = nttHHadrTagged + eventweight[i] * lumi #/ weightsum
	  if tag[i] == VHMETTagged:
	    nVHMETTagged = nVHMETTagged + eventweight[i] * lumi #/ weightsum
	  if tag[i] == Boosted:
	    nBoosted = nBoosted + eventweight[i] * lumi #/ weightsum
	  if tag[i] == gammaHTagged:
	    ngammaHTagged = ngammaHTagged + eventweight[i] * lumi #/ weightsum
	
	count=count+[nUntagged,nVBF1jTagged,nVBF2jTagged,nVHLeptTagged,nVHHadrTagged,nttHLeptTagged,nttHHadrTagged,nVHMETTagged,nBoosted,ngammaHTagged]
	
	trimmed_file=list_file_path.split('/')[-2]+"_"+list_file_path.split('/')[-3]
        
	out.write(trimmed_file+':	'+str.format('{0:.4f}',nUntagged)+'     '+str.format('{0:.4f}',nVBF1jTagged)+'     '+str.format('{0:.4f}',nVBF2jTagged)+'      '+str.format('{0:.4f}',nVHLeptTagged)+'     '+str.format('{0:.4f}',nVHHadrTagged)+'      '+str.format('{0:.4f}',nttHLeptTagged)+'      '+str.format('{0:.4f}',nttHHadrTagged)+'     '+str.format('{0:.4f}',nVHMETTagged)+'      '+str.format('{0:.4f}',nBoosted)+'     '+str.format('{0:.4f}',ngammaHTagged)+'\n')
    

    # Write total event tags on full data set #
    out.write("Total:"+' 	'+str.format('{0:.4f}',count[0])+' 	'+str.format('{0:.4f}',count[1])+' 	'+str.format('{0:.4f}',count[2])+' 	'+str.format('{0:.4f}',count[3])+' 	'+str.format('{0:.4f}',count[4])+' 	'+str.format('{0:.4f}',count[5])+' 	'+str.format('{0:.4f}',count[6])+' 	'+str.format('{0:.4f}',count[7])+' 	'+str.format('{0:.4f}',count[8])+' 	'+str.format('{0:.4f}',count[9]))


if __name__ == "__main__":
    main(sys.argv[1:])
