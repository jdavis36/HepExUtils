import ROOT

def init_LEPSF():
   # 2016 Electrons
   fipEleNotCracks_2016 = "LepSF/ElectronSF_Legacy_2016_NoGap.root"
   root_file = ROOT.TFile.Open(fipEleNotCracks_2016,"READ")
   h_Ele_notCracks_2016 = root_file.Get("EGamma_SF2D").Clone()

   fipEleCracks_2016 = "LepSF/ElectronSF_Legacy_2016_Gap.root"
   root_file = ROOT.TFile.Open(fipEleCracks_2016,"READ")
   h_Ele_Cracks_2016 = root_file.Get("EGamma_SF2D").Clone()

   fipEleReco_highPt_2016 = "LepSF/Ele_Reco_2016.root"
   root_file = ROOT.TFile.Open(fipEleReco_highPt_2016,"READ")
   h_Ele_Reco_highPT_2016 = root_file.Get("EGamma_SF2D").Clone()

   fipEleReco_lowPt_2016 = "LepSF/Ele_Reco_LowEt_2016.root"
   root_file = ROOT.TFile.Open(fipEleReco_lowPt_2016,"READ")
   h_Ele_Reco_lowPT_2016 = root_file.Get("EGamma_SF2D").Clone()

   # 2017 Electrons
   fipEleNotCracks_2017 = "LepSF/ElectronSF_Legacy_2017_NoGap.root"
   root_file = ROOT.TFile.Open(fipEleNotCracks_2017,"READ")
   h_Ele_notCracks_2017 = root_file.Get("EGamma_SF2D").Clone()

   fipEleCracks_2017 = "LepSF/ElectronSF_Legacy_2017_Gap.root"
   root_file = ROOT.TFile.Open(fipEleCracks_2017,"READ")
   h_Ele_Cracks_2017 = root_file.Get("EGamma_SF2D").Clone()

   fipEleReco_highPt_2017 = "LepSF/Ele_Reco_2017.root"
   root_file = ROOT.TFile.Open(fipEleReco_highPt_2017,"READ")
   h_Ele_Reco_highPT_2017 = root_file.Get("EGamma_SF2D").Clone()

   fipEleReco_lowPt_2017 = "LepSF/Ele_Reco_LowEt_2017.root"
   root_file = ROOT.TFile.Open(fipEleReco_lowPt_2017,"READ")
   h_Ele_Reco_lowPT_2017 = root_file.Get("EGamma_SF2D").Clone()


   # 2018 Electrons
   fipEleNotCracks_2018 = "LepSF/ElectronSF_Legacy_2018_NoGap.root"
   root_file = ROOT.TFile.Open(fipEleNotCracks_2018,"READ")
   h_Ele_notCracks_2018 = root_file.Get("EGamma_SF2D").Clone()

   fipEleCracks_2018 = "LepSF/ElectronSF_Legacy_2018_Gap.root"
   root_file = ROOT.TFile.Open(fipEleCracks_2018,"READ")
   h_Ele_Cracks_2018 = root_file.Get("EGamma_SF2D").Clone()

   fipEleReco_highPt_2018 = "LepSF/Ele_Reco_2018.root"
   root_file = ROOT.TFile.Open(fipEleReco_highPt_2018,"READ")
   h_Ele_Reco_highPT_2018 = root_file.Get("EGamma_SF2D").Clone()

   fipEleReco_lowPt_2018 = "LepSF/Ele_Reco_LowEt_2018.root"
   root_file = ROOT.TFile.Open(fipEleReco_lowPt_2018,"READ")
   h_Ele_Reco_lowPT_2018 = root_file.Get("EGamma_SF2D").Clone()

   # 2016 Muons
   fipMu_2016 = "LepSF/final_HZZ_muon_SF_2016RunB2H_legacy_newLoose_newIso_paper.root"
   root_file = ROOT.TFile.Open(fipMu_2016,"READ")
   h_Mu_SF_2016  = root_file.Get("FINAL").Clone()
   h_Mu_Unc_2016 = root_file.Get("ERROR").Clone()

   # 2017 Muons
   fipMu_2017 = "LepSF/final_HZZ_muon_SF_2017_newLooseIso_mupogSysts_paper.root"
   root_file = ROOT.TFile.Open(fipMu_2017,"READ")
   h_Mu_SF_2017  = root_file.Get("FINAL").Clone()
   h_Mu_Unc_2017 = root_file.Get("ERROR").Clone()

   # 2018 Muons
   fipMu_2018 = "LepSF/final_HZZ_muon_SF_2018RunA2D_ER_newLoose_newIso_paper.root"
   root_file = ROOT.TFile.Open(fipMu_2018,"READ")
   h_Mu_SF_2018  = root_file.Get("FINAL").Clone()
   h_Mu_Unc_2018 = root_file.Get("ERROR").Clone()

   print "[LeptonSFHelper] SF maps opened from root files."

def getSF(year, flav, pt, eta, SCeta, isCrack):

   RecoSF = 1.0
   SelSF = 1.0
   SF = 1.0

   print  "year = " + str(year) + " flav = " + str(flav) + " pt = " + str(pt) + " eta = " + str(eta) + " SCeta = " + str(SCeta) + " isCrack = " + str(isCrack) + "\n"

   # Electron reconstruction SFs
   if(abs(flav) == 11):
      if(year == 2016):
         if(pt < 20.):
            RecoSF = h_Ele_Reco_lowPT_2016.GetBinContent(h_Ele_Reco_lowPT_2016.GetXaxis().FindBin(SCeta),h_Ele_Reco_lowPT_2016.GetYaxis().FindBin(15.)) # FIXME: the histogram contains 1 pt bin only
         else:
            RecoSF = h_Ele_Reco_highPT_2016.GetBinContent(h_Ele_Reco_highPT_2016.GetXaxis().FindBin(SCeta),h_Ele_Reco_highPT_2016.GetYaxis().FindBin(min(pt,499.)))
      elif(year == 2017):
         if(pt < 20.):
            RecoSF = h_Ele_Reco_lowPT_2017.GetBinContent(h_Ele_Reco_lowPT_2017.GetXaxis().FindBin(SCeta),h_Ele_Reco_lowPT_2017.GetYaxis().FindBin(15.)) # FIXME: the histogram contains 1 pt bin only
         else:
            RecoSF = h_Ele_Reco_highPT_2017.GetBinContent(h_Ele_Reco_highPT_2017.GetXaxis().FindBin(SCeta),h_Ele_Reco_highPT_2017.GetYaxis().FindBin(min(pt,499.)))
      elif(year == 2018):
         if(pt < 20.):
            RecoSF = h_Ele_Reco_lowPT_2018.GetBinContent(h_Ele_Reco_lowPT_2018.GetXaxis().FindBin(SCeta),h_Ele_Reco_lowPT_2018.GetYaxis().FindBin(15.)) # FIXME: the histogram contains 1 pt bin only
         else:
            RecoSF = h_Ele_Reco_highPT_2018.GetBinContent(h_Ele_Reco_highPT_2018.GetXaxis().FindBin(SCeta),h_Ele_Reco_highPT_2018.GetYaxis().FindBin(min(pt,499.)))
      else:
         print "Ele SFs for " + str(year) + " is not supported!";
         return;
      
      # Electron HZZ selection SF
      if(year == 2016):
         if(isCrack):
            SelSF = h_Ele_Cracks_2016.GetBinContent(h_Ele_Cracks_2016.FindFixBin(SCeta, min(pt,499)))
         else:
            SelSF = h_Ele_notCracks_2016.GetBinContent(h_Ele_notCracks_2016.FindFixBin(SCeta, min(pt,499)))
      elif(year == 2017):
         if(isCrack):
            SelSF = h_Ele_Cracks_2017.GetBinContent(h_Ele_Cracks_2017.FindFixBin(SCeta, min(pt,499)))
         else:
            SelSF = h_Ele_notCracks_2017.GetBinContent(h_Ele_notCracks_2017.FindFixBin(SCeta, min(pt,499)))
      elif(year == 2018):
         if(isCrack):
            SelSF = h_Ele_Cracks_2018.GetBinContent(h_Ele_Cracks_2018.FindFixBin(SCeta, min(pt,499)))
         else:
            SelSF = h_Ele_notCracks_2018.GetBinContent(h_Ele_notCracks_2018.FindFixBin(SCeta, min(pt,499)))
      else:
        print "Ele SFs for " + str(year) + " is not supported!"
        return;
      SF = RecoSF*SelSF
   
   if(abs(flav) == 13 ):
      if(year == 2016):
         SelSF = h_Mu_SF_2016.GetBinContent(h_Mu_SF_2016.GetXaxis().FindBin(eta),h_Mu_SF_2016.GetYaxis().FindBin(min(pt,199)))  #last bin contains the overflow
      elif(year == 2017):
         SelSF = h_Mu_SF_2017.GetBinContent(h_Mu_SF_2017.GetXaxis().FindBin(eta),h_Mu_SF_2017.GetYaxis().FindBin(min(pt,199)))  #last bin contains the overflow
      elif(year == 2018):
         SelSF = h_Mu_SF_2018.GetBinContent(h_Mu_SF_2018.GetXaxis().FindBin(eta),h_Mu_SF_2018.GetYaxis().FindBin(min(pt,199)))  #last bin contains the overflow
      else:
         print "Ele SFs for " + str(year) + " is not supported!"
         return;
      SF = SelSF

   return SF




def fixleptonscalefactor(year, LepLepId, LepPt, LepEta, dataMCWeight):
  init_LEPSF()
  updatedSF = (
    getSF(year, LepLepId[0], LepPt[0], LepEta[0], LepEta[0], False) *
    getSF(year, LepLepId[1], LepPt[1], LepEta[1], LepEta[1], False) *
    getSF(year, LepLepId[2], LepPt[2], LepEta[2], LepEta[2], False) *
    getSF(year, LepLepId[3], LepPt[3], LepEta[3], LepEta[3], False)
  )
  return updatedSF / dataMCWeight
