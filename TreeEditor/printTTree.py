#!/usr/bin/python3
import ROOT

f = ROOT.TFile("/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst/RunIILegacy/200205_CutBased/MC_2017/OffshellAC/gg/ggTo2e2mu_0PMH125_MCFM701/ZZ4lAnalysis.root")
t = f.Get("ZZTree/candTree")

print(t.Print())
