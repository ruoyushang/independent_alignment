
import sys,ROOT
import array
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False) # without this, the histograms returned from a function will be non-type
ROOT.gStyle.SetPaintTextFormat("0.3f")


def MakeEigenvectorPlot(Hist,name):
    
    c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
    pad3 = ROOT.TPad("pad3","pad3",0,0.8,1,1)
    pad3.SetBottomMargin(0.0)
    pad3.SetTopMargin(0.03)
    pad3.SetBorderMode(1)
    pad1 = ROOT.TPad("pad1","pad1",0,0,1,0.8)
    pad1.SetBottomMargin(0.15)
    pad1.SetTopMargin(0.0)
    pad1.SetBorderMode(0)
    pad1.Draw()
    pad3.Draw()

    pad1.cd()

    max_y = 0.
    max_idx = 0.
    for i in range(0,6):
        if max_y<Hist[i].GetMaximum():
            max_y = Hist[i].GetMaximum()
            max_idx = i
    Hist[max_idx].Draw("E")
    for i in range(0,6):
        Hist[i].SetLineColor(i+1)
        Hist[i].Draw("E same")
    c_both.SaveAs('%s.pdf'%(name))


myfile = open('Q_eigenvector.txt')
sensors = 0
edges = 0
modes = 0
misalignment = 0
eigenvector = []
eigenvector += [ROOT.TH1D("eigenvector_1","",16,0,16)]
eigenvector += [ROOT.TH1D("eigenvector_2","",16,0,16)]
eigenvector += [ROOT.TH1D("eigenvector_3","",16,0,16)]
eigenvector += [ROOT.TH1D("eigenvector_4","",16,0,16)]
eigenvector += [ROOT.TH1D("eigenvector_5","",16,0,16)]
eigenvector += [ROOT.TH1D("eigenvector_6","",16,0,16)]
for line in myfile:
    misalignment = float(line)
    eigenvector[sensors].SetBinContent(edges+1,misalignment)
    sensors += 1
    if sensors==6:
        sensors = 0
        edges += 1
    if edges==15:
        MakeEigenvectorPlot(eigenvector,'mode%s'%(modes))
        eigenvector[0].Reset()
        eigenvector[1].Reset()
        eigenvector[2].Reset()
        eigenvector[3].Reset()
        eigenvector[4].Reset()
        eigenvector[5].Reset()
        edges = 0
        modes += 1
