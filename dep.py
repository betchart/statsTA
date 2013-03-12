import ROOT as r
from utils import dependence,symmAnti
r.gROOT.SetBatch(1)

c = r.TCanvas()
c.Divide(2,1)
c.GetPad(1).SetRightMargin(0.2)
c.GetPad(2).SetRightMargin(0.2)
c_symm = r.TCanvas()
c_anti = r.TCanvas()
c_ratio = r.TCanvas()
leptons = ['el','mu']

for el in leptons :
    tfile = r.TFile.Open('data/stats_melded_%s_ph_pn_sn_jn_20.root'%el)
    #d = 'fitTopQueuedBin7TridiscriminantWTopQCD'
    #d = 'fitTopTanhRapiditySum_triD'
    #d = 'fitTopPtOverSumPt_triD'
    d = 'fitTopTanhAvgAbsSumRapidities_triD'
    samples = [key.GetName() for key in tfile.Get(d).GetListOfKeys()]
    fileName = "%s_dep.pdf"%el
    symmName = "%s_symm.pdf"%el
    antiName = "%s_anti.pdf"%el
    ratioName = "%s_ratio.pdf"%el
    c.Print(fileName+'[')
    c_symm.Print(symmName+'[')
    c_anti.Print(antiName+'[')
    c_ratio.Print(ratioName+'[')
    x2ndf = []
    for s in samples :
        h = tfile.Get("%s/%s"%(d,s))
        h_ = h.RebinY(20)
        if "QueuedBin" not in d : h.RebinX()
        #h_dep = dependence(h_,"",3,True)
        h_dep = dependence(h_,"",0.3)
        h_dep.SetTitle(s)
        h_.SetTitle(s)
        h_.SetMinimum(0)
        c.cd(1)
        h_.Draw("colz")
        c.cd(2)
        h_dep.Draw("colz")
        c.Print(fileName)

        oneD = h.ProjectionX()
        symm,anti = symmAnti(oneD)
        symm.SetMinimum(0)

        c_symm.cd()
        symm.Draw()
        c_symm.Print(symmName)

        c_anti.cd()
        anti.Draw()
        c_anti.Print(antiName)

        x2ndf.append("%s: %.2f  %d"%( s, sum([anti.GetBinContent(i)**2/anti.GetBinError(i)**2 for i in range(1,anti.GetNbinsX()+1)]), anti.GetNbinsX()))

        anti.Divide(symm)
        c_ratio.cd()
        anti.Draw()
        c_ratio.Print(ratioName)

    c.Print(fileName+']')
    c_symm.Print(symmName+']')
    c_anti.Print(antiName+']')
    c_ratio.Print(ratioName+']')
    print '\n'.join(x2ndf)
