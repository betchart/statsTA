import sys
import utils
import ROOT as r
r.gStyle.SetOptStat(0)
r.gStyle.SetOptFit(1)

c = r.TCanvas()
tfile = r.TFile.Open(sys.argv[1])




for item in ['ttqq','ttgg','ttqg','ttag']:
    print item
    h = tfile.Get('R03_pileUpRatios2_signalhists_/genTopQueuedBin5_fitTopQueuedBin5_TridiscriminantWTopQCD/%s'%item)
    hxy = h.Project3D('xy')
    hxy_coup = utils.coupling(hxy)
    symm,anti = utils.coupling_symmAnti(hxy_coup)
    anti_devs = r.TH1D('antidevs%s'%item,'', 50, -5, 5)

    for iX in range(2+anti.GetNbinsX()):
        for iY in range(2+anti.GetNbinsY()):
            con,err = anti.GetBinContent(iX,iY), anti.GetBinError(iX,iY)
            if not err: continue
            anti.SetBinContent(iX,iY, con/err)
            anti_devs.Fill(con/err)
            
    anti_devs.Fit('gaus')
    hxy_coup.Draw('colz')
    c.Update()
    raw_input()
    symm.Draw('colz')
    c.Update()
    raw_input()
    anti.Draw('colz')
    c.Update()
    raw_input()
    anti_devs.Draw()
    c.Update()
    raw_input()
    x = symm.ProjectionX()
    x.SetMinimum(0)
    x.Draw()
    c.Update()
    raw_input()
