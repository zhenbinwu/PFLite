from ROOT import *
from sys import argv
import os

NORM=False
LOG=False

gStyle.SetOptStat(False)
gROOT.SetBatch(True)

try: fileNames=argv[1:]
except:
    print "No files specified"
    exit()

if not os.path.exists("plots"): os.mkdir("plots")

files=[]
for fileName in fileNames:
    files.append(TFile(fileName))

keys=files[0].GetListOfKeys()

c=TCanvas()
for key in keys:
    if type(files[0].Get(key.GetName()))!=type(TH1F()) and \
       type(files[0].Get(key.GetName()))!=type(TProfile()) : continue

    l=TLegend(.49,.69,.89,.89)
    l.SetFillStyle(0)
    hists=[]
    max=0

    for lp in range(len(files)):
        file=files[lp]
        fileName=fileNames[lp]

        file.cd()
        h=file.Get(key.GetName())
        #h.Rebin(4)
        hists.append(h)

        if NORM:h.Scale(1./h.Integral())

        if h.GetMaximum()>max: max=h.GetMaximum()

        h.SetLineColor(1+lp)
        h.SetLineWidth(3)
        legname = fileNames[lp].split('/')[1]
        legname = legname.split('.')[0]
        l.AddEntry(h,legname.replace('_', ' '),"l")

    for lp in range(len(hists)):
        h=hists[lp]

        if lp==0:
            if LOG or key.GetName() == "jet_eta":
              c.SetLogy(True)
              h.SetMaximum(100*max)
            else:
              c.SetLogy(False)
              h.SetMaximum(2*max)
            h.Draw("hist")
            #if type(h) == type(TProfile):
                #h.Draw("p")
            #else:
                #h.Draw()
        else:
            h.Draw("histSAME")
            #if type(h) == type(TProfile):
                #h.Draw("SAMEp")
            #else:
                #h.Draw("SAME")
    l.Draw("SAME")
    c.SaveAs("plots/"+key.GetName()+".pdf")
    c.SaveAs("plots/"+key.GetName()+".png")

