import os
import ROOT

# This initializes discriminat splines for categorizations #
def init_spline():
  D2jetZHSpline = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_DjjZH13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  D2jetWHSpline = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_DjjWH13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  D2jetVBFSpline = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_DjjVBF13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  D1jetVBFSpline = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_DjVBF13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  DbkgkinSpline4e = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_Dbkgkin_4e13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  DbkgkinSpline4mu = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_Dbkgkin_4mu13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  DbkgkinSpline2e2mu = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_Dbkgkin_2e2mu13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  DbkgjjEWQCDSpline4lHadVH = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_HadVHTagged_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  DbkgjjEWQCDSpline2l2lHadVH = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_HadVHTagged_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  DbkgjjEWQCDSpline4lJJVBF = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_DbkgjjEWQCD_4l_JJVBFTagged_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  DbkgjjEWQCDSpline2l2lJJVBF = ROOT.TFile("Discriminants/SmoothKDConstant_m4l_DbkgjjEWQCD_2l2l_JJVBFTagged_13TeV.root").Get("sp_gr_varReco_Constant_Smooth")
  array_splines=[D2jetZHSpline,D2jetWHSpline,D2jetVBFSpline,DbkgkinSpline4e,DbkgkinSpline4mu,DbkgkinSpline2e2mu,DbkgjjEWQCDSpline4lHadVH,DbkgjjEWQCDSpline2l2lHadVH,DbkgjjEWQCDSpline4lJJVBF,DbkgjjEWQCDSpline2l2lJJVBF]
  return array_splines


def getDVBF2jetsConstant(DVBF2jetsSpline,ZZMass):
  return DVBF2jetsSpline.Eval(ZZMass)
def getDVBF1jetConstant(DVBF1jetSpline,ZZMass):
  return DVBF1jetSpline.Eval(ZZMass)
def getDWHhConstant(DWHhSpline, ZZMass):
  return DWHhSpline.Eval(ZZMass)
def getDZHhConstant(DZHhSpline,ZZMass):
  return DZHhSpline.Eval(ZZMass)
def getDVBF2jetsWP( ZZMass,  useQGTagging):
  if (useQGTagging):
    assert 0
    return 0.363
  else:
    return 0.46386

def getDVBF1jetWP( ZZMass,  useQGTagging):
  if (useQGTagging):
    assert 0 
    return 0.716
  else:
    #return 0.37605;
    return 0.58442
def getDWHhWP(ZZMass,useQGTagging):
  if (useQGTagging):
    assert(0)
    return 0.965
  else:
    return 0.88384

def getDZHhWP(ZZMass,useQGTagging):
  if (useQGTagging):
    assert 0
    return 0.9952
  else:
    return 0.91315

def getDVBF2jetsConstant_shiftWP(DVBF2jetsSpline, ZZMass,  useQGTagging,  newWP) :
  oldc = getDVBF2jetsConstant(DVBF2jetsSpline,ZZMass)
  oldWP = getDVBF2jetsWP(ZZMass, useQGTagging)
  return oldc * (oldWP/newWP) * ((1-newWP)/(1-oldWP))

def getDVBF1jetConstant_shiftWP(DVBF1jetSpline, ZZMass,  useQGTagging,  newWP) :
  oldc = getDVBF1jetConstant(DVBF1jetSpline,ZZMass)
  oldWP = getDVBF1jetWP(ZZMass, useQGTagging)
  return oldc * (oldWP/newWP) * ((1-newWP)/(1-oldWP))

def getDWHhConstant_shiftWP(DWHhSpline, ZZMass,  useQGTagging,  newWP):
  oldc = getDWHhConstant(DWHhSpline,ZZMass)
  oldWP = getDWHhWP(ZZMass, useQGTagging)
  return oldc * (oldWP/newWP) * ((1-newWP)/(1-oldWP))

def getDZHhConstant_shiftWP( ZZMass,  useQGTagging,  newWP):
  oldc = getDZHhConstant(ZZMass)
  oldWP = getDZHhWP(ZZMass, useQGTagging)
  return oldc * (oldWP/newWP) * ((1-newWP)/(1-oldWP))

def getDbkgVBFdecConstant( ZZflav,  ZZMass): # ZZflav==id1*id2*id3*id4
  if (abs(ZZflav)==11*11*11*11 or abs(ZZflav)==2*11*11*11*11 or abs(ZZflav)==2*11*11*2*11*11): 
    return DbkgVBFdecSpline4l.Eval(ZZMass)  
  if (abs(ZZflav)==11*11*13*13 or abs(ZZflav)==2*11*11*13*13 or abs(ZZflav)==2*11*11*2*13*13): 
    return DbkgVBFdecSpline2l2l.Eval(ZZMass)
  if (abs(ZZflav)==13*13*13*13 or abs(ZZflav)==2*13*13*13*13 or abs(ZZflav)==2*13*13*2*13*13): 
    return DbkgVBFdecSpline4l.Eval(ZZMass)
  print "Invalid ZZflav " + str(ZZflav)
  assert 0
  return 0

def getDbkgVHdecConstant( ZZflav,  ZZMass): # ZZflav==id1*id2*id3*id4
  if (abs(ZZflav)==11*11*11*11 or abs(ZZflav)==2*11*11*11*11 or abs(ZZflav)==2*11*11*2*11*11): 
    return DbkgVHdecSpline4l.Eval(ZZMass);
  if (abs(ZZflav)==11*11*13*13 or abs(ZZflav)==2*11*11*13*13 or abs(ZZflav)==2*11*11*2*13*13): 
    return DbkgVHdecSpline2l2l.Eval(ZZMass)
  if (abs(ZZflav)==13*13*13*13 or abs(ZZflav)==2*13*13*13*13 or abs(ZZflav)==2*13*13*2*13*13): 
    return DbkgVHdecSpline4l.Eval(ZZMass)
  print "Invalid ZZflav " + str(ZZflav)
  assert 0
  return 0

def getDbkgkinConstant( ZZflav,  ZZMass): # ZZflav==id1*id2*id3*id4
  if (abs(ZZflav)==11*11*11*11 or abs(ZZflav)==2*11*11*11*11 or abs(ZZflav)==2*11*11*2*11*11):
     return DbkgkinSpline4e.Eval(ZZMass)
  if (abs(ZZflav)==11*11*13*13 or abs(ZZflav)==2*11*11*13*13 or abs(ZZflav)==2*11*11*2*13*13):
     return DbkgkinSpline2e2mu.Eval(ZZMass)
  if (abs(ZZflav)==13*13*13*13 or abs(ZZflav)==2*13*13*13*13 or abs(ZZflav)==2*13*13*2*13*13):
     return DbkgkinSpline4mu.Eval(ZZMass)
  print "Invalid ZZflav " + str(ZZflav) 
  assert 0
  return 0

def getDbkgConstant( ZZflav,  ZZMass):
  return getDbkgkinConstant(ZZflav, ZZMass)

