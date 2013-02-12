import ROOT as r
from utils import dependence
r.gROOT.SetBatch(1)

c = r.TCanvas()
c.Divide(2,1)
c.GetPad(1).SetRightMargin(0.2)
c.GetPad(2).SetRightMargin(0.2)
leptons = ['el','mu']

for el in leptons :
    tfile = r.TFile.Open('data/stats_melded_%s_ph_c_20.root'%el)
    d = 'fitTopQueuedBin7TridiscriminantWTopQCD'
    samples = [key.GetName() for key in tfile.Get(d).GetListOfKeys()]
    fileName = "%s_dep.pdf"%el
    c.Print(fileName+'[')
    for s in samples :
        h = tfile.Get("%s/%s"%(d,s))
        h_ = h.RebinY(20)
        #h_dep = dependence(h_,"",3,True)
        h_dep = dependence(h_)
        h_dep.SetTitle(s)
        h_.SetTitle(s)
        h_.SetMinimum(0)
        c.cd(1)
        h_.Draw("colz")
        c.cd(2)
        h_dep.Draw("colz")
        c.Print(fileName)
    c.Print(fileName+']')
        
