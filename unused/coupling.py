import sys
import lib
import ROOT as r
r.gStyle.SetOptStat(0)
r.gStyle.SetOptFit(1)

tfile = r.TFile.Open(sys.argv[1])
c = r.TCanvas()
c.Divide(2,2)
c.cd(1)

zero = r.TF1('zero','0',-1,1)
one = r.TF1('one','1',-1,1)
one.SetLineWidth(1)
one.SetLineStyle(r.kDashed)

for item in ['ttqq','ttgg','ttqg','ttag']:
    h = tfile.Get('R03_pileUpRatios2_signalhists_/genTopQueuedBin5_fitTopQueuedBin5_TridiscriminantWTopQCD/%s'%item)

    for iC in range(1,1+h.GetZaxis().GetNbins()):
        h.GetZaxis().SetRange(iC,iC)
        hxy = h.Project3D(('exy%d'%iC))
        hxy.SetTitle('%s %d;gen;reco'%(item,iC))
        hxy_coup = lib.coupling(hxy)

        symm,anti = lib.coupling_symmAnti(hxy_coup)
        anti_devs = r.TH1D('antidevs%s%d'%(item,iC),'', 50, -5, 5)
        anti_rel = anti.Clone('_rel')

        for iX in range(2+anti.GetNbinsX()):
            for iY in range(2+anti.GetNbinsY()):
                con,err = anti.GetBinContent(iX,iY), anti.GetBinError(iX,iY)
                if not err: continue
                anti_rel.SetBinContent(iX,iY, con/err)
                anti_devs.Fill(con/err)
 
        anti.Divide(symm)
        anti_devs.Fit('gaus','Q')
        c.cd(1)
        symm.Draw('colz')
        ##
        c.cd(2)
        anti_rel.Draw('colz')
        ##
        c.cd(3)
        yorig = hxy_coup.ProjectionY()
        y = symm.ProjectionY()
        y.SetMinimum(0)
        y.SetMaximum(1.1*max(y.GetMaximum(),yorig.GetMaximum()))
        y.SetLineColor(r.kRed)
        y.SetLineWidth(2)
        y.Draw('hist')
        one.Draw('same')
        yorig.SetMarkerSize(1)
        yorig.SetMarkerStyle(20)
        yorig.Draw('same')
        ##
        c.cd(4)
        ydiff = yorig.Clone('diff')
        ydiff.Add(y,-1)
        ydiff.Draw()
        zero.Draw('same')
        c.Update()
        print item,iC, [anti_devs.GetFunction('gaus').GetParameter(i) for i in [1,2]] 
        raw_input()
