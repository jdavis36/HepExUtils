#!/usr/bin/python3

import sys, getopt
import os
import glob
import ROOT
from math import sqrt
import time
#from os import path
from pathlib import Path
import re
from tqdm import trange, tqdm
import numpy as np
from array import *
import random
from collections import Counter
from decimal import *

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

    #=========================== Set input file path and output file path ===========================
    
    filename = inputfile #treelist[numfile]
    branchlistpath = branchfile
    tagtreepath = outputfile

    ind = filename.split("/").index("200205_CutBased")

    tagtreefile = "/".join(filename.split("/")[ind:])
    tagtreefilename = tagtreepath+tagtreefile

    print(filename, "\n")
    print(tagtreefilename, "\n")

    datatype = {'vector<float>':'vector<float>','vector<short>':'vector<short>','vector<long>':'vector<long>','Char_t':'B','UChar_t':'b','Bool_t':'O','Short_t':'S','UShort_t':'s','Counter':'I','Int_t':'I','UInt_t':'i','Double_t':'D','Double32_t':'d','Float_t':'F','Float16_t':'f','Long_t':'G','ULong_t':'g','Long64_t':'L','ULong64_t':'l','CharStar':'C'}

    D2jetZHSpline = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_DjjZH_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
    D2jetWHSpline = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_DjjWH_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
    D2jetVBFSpline = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_DjjVBF_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")

    DbkgkinSpline4e = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_Dbkgkin_4e_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
    DbkgkinSpline4mu = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_Dbkgkin_4mu_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
    DbkgkinSpline2e2mu = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_Dbkgkin_2e2mu_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")

    DggbkgkinSpline4e = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_Dggbkgkin_4e_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
    DggbkgkinSpline4mu = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_Dggbkgkin_4mu_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
    DggbkgkinSpline2e2mu = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_Dggbkgkin_2e2mu_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")

    DbkgjjEWQCDSpline4lHadVH = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_HadVHTagged_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
    DbkgjjEWQCDSpline2l2lHadVH = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_HadVHTagged_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")

    DbkgjjEWQCDSpline4lJJVBF = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_JJVBFTagged_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
    DbkgjjEWQCDSpline2l2lJJVBF = ROOT.TFile("/eos/user/l/lkang/Active_Research/Discriminants/RecoMEConstants/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_JJVBFTagged_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")

    #WPCshift = 1
    WPCshift2jv = 0.46386/(1. - 0.46386)
    WPCshift2jz = 0.91315/(1. - 0.91315)
    WPCshift2jw = 0.88384/(1. - 0.88384)

    #=========================== Check existence of output and set up target branches ===========================
   
    print("================ Check output location and set up branches ================\n")
 
    if Path(tagtreefilename).exists():
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
        branchdict["Dbkg"] = [] #np.array([])
        branchdict["Dbsi"] = [] #np.array([])
        branchdict["Bin40"] = [] #np.array([])
        for branch in branchlist:
            branchdict[branch] = [] #np.array([])
        
        for tree in ["candTree", "candTree_failed"]:
        #for tree in ["candTree"]:
            print("\n================ Reading events from " + tree + " ================\n")

            f = ROOT.TFile(filename)
            t = f.Get("ZZTree/"+tree)

            treebranches = [ x.GetName() for x in t.GetListOfBranches() ]
            
            for ent in trange(t.GetEntries()):
                
                #=========================== Loop over events ===========================
                
                while t.GetEntry(ent):
                    
                    #=========================== Fill failed events with dummy and skip to loop over branches ===========================
                    
                    branchdict["Bin40"].append(f.Get("ZZTree/Counters").GetBinContent(40)) #= np.append(branchdict["Bin40"], f.Get("ZZTree/Counters").GetBinContent(40))
                    if tree == "candTree_failed":
                        branchdict["EventTag"].append(-999) #= np.append(branchdict["EventTag"], -999)
                        branchdict["Dbkg"].append(-999) # = np.append(branchdict["Dbkg"], -999)
                        branchdict["Dbsi"].append(-999) # = np.append(branchdict["Dbsi"], -999)
                        break
                        
                    #=========================== Tagging event by category ===========================

                    M = t.ZZMass
    #                 M = t.GenHMass

                    tag = "none"

                    const2jv = WPCshift2jv * D2jetVBFSpline.Eval(M)
                    vars2jv = [t.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal, t.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal]
                    if not any(var < 0 for var in vars2jv):
                        if (vars2jv[0]+const2jv*vars2jv[1]) != 0:
                            D2jv = vars2jv[0]/(vars2jv[0]+const2jv*vars2jv[1])
                            if D2jv > 0.5:
                                if (((t.nCleanedJetsPt30>1) and (t.nCleanedJetsPt30<4)) and (t.nCleanedJetsPt30BTagged<2)) or ((t.nCleanedJetsPt30>3) and (t.nCleanedJetsPt30BTagged==0)):
                                    tag = "VBF"
                        else: 
                            D2jv = -999


                    const2jz = WPCshift2jz * D2jetZHSpline.Eval(M)
                    vars2jz = [t.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal, t.p_HadZH_mavjj_JECNominal, t.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, t.p_HadZH_mavjj_true_JECNominal]
                    if not any(var < 0 for var in vars2jz):
                        if (vars2jz[0]*vars2jz[1]+const2jz*vars2jz[2]*vars2jz[3]) != 0:
                            D2jz = vars2jz[0]*vars2jz[1]/(vars2jz[0]*vars2jz[1]+const2jz*vars2jz[2]*vars2jz[3])
                            if D2jz > 0.5:
                                if ((t.nCleanedJetsPt30>1) and (t.nCleanedJetsPt30<4)) or ((t.nCleanedJetsPt30>3) and (t.nCleanedJetsPt30BTagged==0)):
                                    tag = "VH"
                        else: 
                            D2jz = -999


                    const2jw = WPCshift2jw * D2jetWHSpline.Eval(M)
                    vars2jw = [t.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal, t.p_HadWH_mavjj_JECNominal, t.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, t.p_HadWH_mavjj_true_JECNominal]
                    if not any(var < 0 for var in vars2jw):
                        if (vars2jw[0]*vars2jw[1]+const2jw*vars2jw[2]*vars2jw[3]) != 0:
                            D2jw = vars2jw[0]*vars2jw[1]/(vars2jw[0]*vars2jw[1]+const2jw*vars2jw[2]*vars2jw[3])
                            if D2jw > 0.5:
                                if ((t.nCleanedJetsPt30>1) and (t.nCleanedJetsPt30<4)) or ((t.nCleanedJetsPt30>3) and (t.nCleanedJetsPt30BTagged==0)):
                                    tag = "VH"
                        else: 
                            D2jw = -999
                    
                    #=========================== Saving category tag ===========================
                    
                    if tag == "VBF": branchdict["EventTag"].append(1) #= np.append(branchdict["EventTag"], 1)
                    elif tag == "VH": branchdict["EventTag"].append(2) # = np.append(branchdict["EventTag"], 2)
                    else: branchdict["EventTag"].append(0) # = np.append(branchdict["EventTag"], 0)
                                           
                    #=========================== Calculating EW discriminants ===========================

                    ZZflav = t.Z1Flav * t.Z2Flav
#                     ZZflav = t.GenZ1Flav * t.GenZ2Flav
                    
                    if tag == "VBF" or tag == "VH":

                        if tag == "VH":
                            if (abs(ZZflav)==11*11*11*11 or abs(ZZflav)==2*11*11*11*11 or abs(ZZflav)==2*11*11*2*11*11): const1 = DbkgjjEWQCDSpline4lHadVH.Eval(M)
                            elif (abs(ZZflav)==11*11*13*13 or abs(ZZflav)==2*11*11*13*13 or abs(ZZflav)==2*11*11*2*13*13): const1 = DbkgjjEWQCDSpline2l2lHadVH.Eval(M)
                            elif (abs(ZZflav)==13*13*13*13 or abs(ZZflav)==2*13*13*13*13 or abs(ZZflav)==2*13*13*2*13*13): const1 = DbkgjjEWQCDSpline4lHadVH.Eval(M)
                            else:
                                branchdict["Dbkg"].append(-999) # = np.append(branchdict["Dbkg"], -999)
                                branchdict["Dbsi"].append(-999) # = np.append(branchdict["Dbsi"], -999)
                                break
                                
                        elif tag == "VBF":
                            if (abs(ZZflav)==11*11*11*11 or abs(ZZflav)==2*11*11*11*11 or abs(ZZflav)==2*11*11*2*11*11): const1 = DbkgjjEWQCDSpline4lJJVBF.Eval(M)
                            elif (abs(ZZflav)==11*11*13*13 or abs(ZZflav)==2*11*11*13*13 or abs(ZZflav)==2*11*11*2*13*13): const1 = DbkgjjEWQCDSpline2l2lJJVBF.Eval(M)
                            elif (abs(ZZflav)==13*13*13*13 or abs(ZZflav)==2*13*13*13*13 or abs(ZZflav)==2*13*13*2*13*13): const1 = DbkgjjEWQCDSpline4lJJVBF.Eval(M)
                            else:
                                branchdict["Dbkg"].append(-999) # = np.append(branchdict["Dbkg"], -999)
                                branchdict["Dbsi"].append(-999) # = np.append(branchdict["Dbsi"], -999)
                                break

                        var = [t.p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal, t.p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal, t.p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal, t.p_JJVBF_BKG_MCFM_JECNominal, t.p_HadZH_BKG_MCFM_JECNominal, t.p_HadWH_BKG_MCFM_JECNominal, t.p_JJQCD_BKG_MCFM_JECNominal, t.p_HadZH_mavjj_JECNominal, t.p_HadZH_mavjj_true_JECNominal, t.p_HadWH_mavjj_JECNominal, t.p_HadWH_mavjj_true_JECNominal, t.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal, t.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal, t.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal, t.pConst_JJVBF_BKG_MCFM_JECNominal, t.pConst_HadZH_BKG_MCFM_JECNominal, t.pConst_HadWH_BKG_MCFM_JECNominal, t.pConst_JJQCD_BKG_MCFM_JECNominal]
                        
                        if any(v == 0 for v in [var[11], var[12], var[13], var[14], var[15], var[16], var[17], var[8], var[10]]) or any(v < 0 for v in var): 
                            D1 = -999

                        else:
                            vbf = var[0]/var[11]
                            zh = var[1]/var[12]
                            wh = var[2]/var[13]
                            constA = 1./(1./var[11]+1./var[12]+1./var[13])

                            vbs = var[3]/var[14]
                            zzz = var[4]/var[15]
                            wzz = var[5]/var[16]
                            qcdzz = var[6]/var[17]
                            constB = 1./(1./var[14]+1./var[15]+1./var[16]+1./var[17])

                            scale_Pmjj_vb=1
                            scale_Pmjj_z = var[7]/var[8]
                            scale_Pmjj_w = var[9]/var[10]

                            vbf *= scale_Pmjj_vb
                            vbs *= scale_Pmjj_vb

                            zh *= scale_Pmjj_z
                            zzz *= scale_Pmjj_z

                            wh *= scale_Pmjj_w
                            wzz *= scale_Pmjj_w

                            PA = (vbf + zh + wh)*constA
                            PB = (vbs + zzz + wzz + qcdzz)*constB

                            d = (PA+const1*PB)
                            if d != 0:
                                D1 = PA/d
                            else:
                                D1 = -999

                        if tag == "VH":
                            if (abs(ZZflav)==11*11*11*11 or abs(ZZflav)==2*11*11*11*11 or abs(ZZflav)==2*11*11*2*11*11): const2 = DbkgjjEWQCDSpline4lHadVH.Eval(M)
                            elif (abs(ZZflav)==11*11*13*13 or abs(ZZflav)==2*11*11*13*13 or abs(ZZflav)==2*11*11*2*13*13): const2 = DbkgjjEWQCDSpline2l2lHadVH.Eval(M)
                            elif (abs(ZZflav)==13*13*13*13 or abs(ZZflav)==2*13*13*13*13 or abs(ZZflav)==2*13*13*2*13*13): const2 = DbkgjjEWQCDSpline4lHadVH.Eval(M)
                            else: 
                                branchdict["Dbkg"].append(-999) # = np.append(branchdict["Dbkg"], -999)
                                branchdict["Dbsi"].append(-999) # = np.append(branchdict["Dbsi"], -999)
                                break

                        elif tag == "VBF":
                            if (abs(ZZflav)==11*11*11*11 or abs(ZZflav)==2*11*11*11*11 or abs(ZZflav)==2*11*11*2*11*11): const2 = DbkgjjEWQCDSpline4lJJVBF.Eval(M)
                            elif (abs(ZZflav)==11*11*13*13 or abs(ZZflav)==2*11*11*13*13 or abs(ZZflav)==2*11*11*2*13*13): const2 = DbkgjjEWQCDSpline2l2lJJVBF.Eval(M)
                            elif (abs(ZZflav)==13*13*13*13 or abs(ZZflav)==2*13*13*13*13 or abs(ZZflav)==2*13*13*2*13*13): const2 = DbkgjjEWQCDSpline4lJJVBF.Eval(M)
                            else: 
                                branchdict["Dbkg"].append(-999) # = np.append(branchdict["Dbkg"], -999)
                                branchdict["Dbsi"].append(-999) # = np.append(branchdict["Dbsi"], -999)
                                break


                        var = [t.p_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal, t.p_HadZH_S_SIG_ghz1_1_MCFM_JECNominal, t.p_HadWH_S_SIG_ghw1_1_MCFM_JECNominal, t.p_JJVBF_BKG_MCFM_JECNominal, t.p_HadZH_BKG_MCFM_JECNominal, t.p_HadWH_BKG_MCFM_JECNominal, t.p_JJQCD_BKG_MCFM_JECNominal, t.p_HadZH_mavjj_JECNominal, t.p_HadZH_mavjj_true_JECNominal, t.p_HadWH_mavjj_JECNominal, t.p_HadWH_mavjj_true_JECNominal, t.pConst_JJVBF_S_SIG_ghv1_1_MCFM_JECNominal, t.pConst_HadZH_S_SIG_ghz1_1_MCFM_JECNominal, t.pConst_HadWH_S_SIG_ghw1_1_MCFM_JECNominal, t.pConst_JJVBF_BKG_MCFM_JECNominal, t.pConst_HadZH_BKG_MCFM_JECNominal, t.pConst_HadWH_BKG_MCFM_JECNominal, t.pConst_JJQCD_BKG_MCFM_JECNominal, t.p_JJVBF_S_BSI_ghv1_1_MCFM_JECNominal, t.p_HadZH_S_BSI_ghz1_1_MCFM_JECNominal, t.p_HadWH_S_BSI_ghw1_1_MCFM_JECNominal]
                        
                        if any(v == 0 for v in [var[11], var[12], var[13], var[14], var[15], var[16], var[17], var[8], var[10]]) or any(v < 0 for v in var): 
                            D2 = -999

                        else:
                            vbf = var[0]/var[11]
                            zh = var[1]/var[12]
                            wh = var[2]/var[13]
                            constA = 1./(1./var[11]+1./var[12]+1./var[13])

                            vbs = var[3]/var[14]
                            zzz = var[4]/var[15]
                            wzz = var[5]/var[16]
                            qcdzz = var[6]/var[17]
                            constB = 1./(1./var[14]+1./var[15]+1./var[16]+1./var[17])

                            vbf_vbs_int = var[18]*(1./var[11]+1./var[14]) - vbf - vbs
                            zh_zzz_int = var[19]*(1./var[12]+1./var[15]) - zh - zzz
                            wh_wzz_int = var[20]*(1./var[13]+1./var[16]) - wh - wzz

                            scale_Pmjj_vb=1
                            scale_Pmjj_z = var[7]/var[8]
                            scale_Pmjj_w = var[9]/var[10]

                            vbf *= scale_Pmjj_vb
                            vbs *= scale_Pmjj_vb
                            vbf_vbs_int *= scale_Pmjj_vb

                            zh *= scale_Pmjj_z
                            zzz *= scale_Pmjj_z
                            zh_zzz_int *= scale_Pmjj_z

                            wh *= scale_Pmjj_w
                            wzz *= scale_Pmjj_w
                            wh_wzz_int *= scale_Pmjj_w

                            PA = (vbf + zh + wh)*constA
                            PB = (vbs + zzz + wzz + qcdzz)*constB
                            d = (PA+const2*PB)
                            if (d != 0) and (constA*constB >= 0) and (const2 >= 0):
                                Pint = (vbf_vbs_int + zh_zzz_int + wh_wzz_int)*sqrt(constA*constB)
                                D2 = Pint*sqrt(const2)/d
                            else:
                                D2 = -999

                    #=========================== Calculating gg discriminants ===========================

                    elif tag == "none":

                        if (abs(ZZflav)==11*11*11*11 or abs(ZZflav)==2*11*11*11*11 or abs(ZZflav)==2*11*11*2*11*11): const1 = DbkgkinSpline4e.Eval(M)
                        elif (abs(ZZflav)==11*11*13*13 or abs(ZZflav)==2*11*11*13*13 or abs(ZZflav)==2*11*11*2*13*13): const1 = DbkgkinSpline2e2mu.Eval(M)
                        elif (abs(ZZflav)==13*13*13*13 or abs(ZZflav)==2*13*13*13*13 or abs(ZZflav)==2*13*13*2*13*13): const1 = DbkgkinSpline4mu.Eval(M)
                        else: 
                            branchdict["Dbkg"].append(-999) # = np.append(branchdict["Dbkg"], -999)
                            branchdict["Dbsi"].append(-999) # = np.append(branchdict["Dbsi"], -999)
                            break

                        var = [t.p_GG_SIG_ghg2_1_ghz1_1_JHUGen, t.p_QQB_BKG_MCFM]
                        if any(v < 0 for v in var): 
                            D1 = -999

                        else:
                            d = (var[0]+const1*var[1])
                            if d != 0:
                                D1 = var[0]/d
                            else:
                                D1 = -999

                        if (abs(ZZflav)==11*11*11*11 or abs(ZZflav)==2*11*11*11*11 or abs(ZZflav)==2*11*11*2*11*11): const2 = DggbkgkinSpline4e.Eval(M)
                        elif (abs(ZZflav)==11*11*13*13 or abs(ZZflav)==2*11*11*13*13 or abs(ZZflav)==2*11*11*2*13*13): const2 = DggbkgkinSpline2e2mu.Eval(M)
                        elif (abs(ZZflav)==13*13*13*13 or abs(ZZflav)==2*13*13*13*13 or abs(ZZflav)==2*13*13*2*13*13): const2 = DggbkgkinSpline4mu.Eval(M)
                        else: 
                                branchdict["Dbkg"].append(-999) # = np.append(branchdict["Dbkg"], -999)
                                branchdict["Dbsi"].append(-999) # = np.append(branchdict["Dbsi"], -999)
                                break


                        var = [t.p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM, t.p_GG_BKG_MCFM, t.p_GG_BSI_kappaTopBot_1_ghz1_1_MCFM, t.pConst_GG_SIG_kappaTopBot_1_ghz1_1_MCFM, t.pConst_GG_BKG_MCFM]
                        if any(v == 0 for v in var[3:]) or any(v < 0 for v in var[:3]): 
                            D2 = -999

                        else:
                            d = (var[0]+const2*var[1])
                            if (d != 0) and (var[3]*var[4]*const2 >= 0):
                                D2 = (var[2]*(1./var[3]+1./var[4])-var[0]/var[3]-var[1]/var[4])*sqrt(var[3]*var[4]*const2)/d
                            else:
                                D2 = -999
                        
                    #=========================== Saving calculated discriminants ===========================
                    
                    branchdict["Dbkg"].append(D1) # = np.append(branchdict["Dbkg"], D1)
                    branchdict["Dbsi"].append(D2) # = np.append(branchdict["Dbsi"], D2)
                        
                    break
                    
                #=========================== Loop over branches to copy values ===========================

                while t.GetEntry(ent):
                    for branchcand in branchlist:
                        if branchcand in treebranches:
                            exec("global datatypecode; datatypecode = datatype[t.GetLeaf('"+branchcand+"').GetTypeName()]")
                            exec("global branchval; branchval = t." + branchcand)
                            branchdict[branchcand].append(branchval) # = np.append(branchdict[branchcand], branchval)
                        else:
                            branchdict[branchcand].append(-999) # = np.append(branchdict[branchcand], -999)
                    break
        
        Path("/".join(tagtreefilename.split("/")[:-1])).mkdir( 0o755, True, True )
        ftagtree = ROOT.TFile(tagtreefilename, "CREATE")
        tagtree = ROOT.TTree("eventTree", "eventTree")
        
        print("\n================ Create new TTree and fill with branches ================\n")
        
        for key in branchdict.keys():
            print(len(branchdict[key]), key)
            exec("value"+key.replace('-', 'm')+" = array('"+datatypecode.replace("F","f")+"', [0])")
            exec("tagtree.Branch('"+key.replace('-', 'm')+"', value"+key.replace('-', 'm')+", '"+key.replace('-', 'm')+"/"+datatypecode+"')")
        
        print("\n================ Fill with branches in event loop ================\n")
        
        for elem in trange(len(branchdict[key])):
            for key in branchdict.keys():
                exec("value"+key.replace('-', 'm') +"[0] = branchdict[key][elem]")
            tagtree.Fill()

        print("\n================ Check new TTree branches ================\n")
            
        tagtree.Print()
        tagtree.Write("", ROOT.TObject.kOverwrite)
        ftagtree.Close()

        print("\n", Counter(branchdict["EventTag"]))

if __name__ == "__main__":
    main(sys.argv[1:])
