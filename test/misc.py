import ROOT as r
import random

hist = r.TH1D('hist','hist',100,0,1)
for i in range(int(1e6)) : hist.Fill(random.gauss(0.5,0.1))

x = r.RooRealVar('x','x',0,1)

roohist = r.RooDataHist('roohist','roohist',r.RooArgList(x),hist)
pdf = r.RooHistPdf('pdf','pdf',r.RooArgSet(x), roohist, 0)

y = r.RooRealVar('y','y',0,20)
pois = r.RooPoisson('pois','pois',y, r.RooFit.RooConst(5.5))
pois2 = r.RooPoisson('pois2','pois2',y, r.RooFit.RooConst(6))
pois3 = r.RooPoisson('pois3','pois3',y, r.RooFit.RooConst(5))

pl = x.frame()
pdf.plotOn(pl)
pl.Draw()

raw_input()

plot = y.frame()

pois.plotOn(plot)
pois2.plotOn(plot,r.RooFit.LineColor(r.kRed))
pois3.plotOn(plot,r.RooFit.LineColor(r.kGreen))
plot.Draw()



raw_input()
