import ROOT as r
import utils

def flows(h):
    nbins = h.GetNbinsX()
    utils.combineBinContentAndError(h,1,0)
    utils.combineBinContentAndError(h,nbins,nbins+1)

r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L tdrstyle.C")
r.setTDRStyle()
r.gStyle.SetOptFit(0)
r.tdrStyle.SetPadRightMargin(0.06)

tfiles = [r.TFile.Open('data/extra_control_top_mu_ph_sn_jn_%d.root'%d) for d in [20,30]]
tt = tuple([tfile.Get('qRecoilPt/tt') for tfile in tfiles]  )
gg = tuple([tfile.Get('qRecoilPt/ttgg') for tfile in tfiles])
qq = tuple([tfile.Get('qRecoilPt/ttqq') for tfile in tfiles])
qg = tuple([tfile.Get('qRecoilPt/ttqg') for tfile in tfiles][::-1])
for hs in [tt,gg,qq,qg]: 
    for h in hs:
        h.UseCurrentStyle()
        h.SetLineWidth(3)
        h.GetXaxis().SetTitle('recoiling quark p_{T} (GeV)')
        h.GetYaxis().SetTitle('events / bin / pb^{-1}')
        flows(h)

for color,(h1,h2) in zip([r.kGray,r.kRed,r.kBlue,r.kGreen],[tt,qq,gg,qg]) :
    h1.SetLineColor(color)
    h2.SetLineColor(color)
    h2.SetLineStyle(r.kDashed)

c=r.TCanvas()
leg = r.TLegend(0.5,0.7,0.9,0.9)
leg.SetBorderSize(0)
leg.SetFillColor(r.kWhite)
leg.AddEntry(tt[0],'total t#bar{t}q (POWHEG)','l')
tt[0].Draw('hist')
for (h,k),label in zip([qg,gg,qq],
                       ['qg#rightarrow t#bar{t}q','gg#rightarrow t#bar{t}q', 'q#bar{q}#rightarrow t#bar{t}q']):
    h.Draw('hist same')
    k.Draw('hist same')
    leg.AddEntry(h,label,'l')
qg[1].Draw('hist same')
leg.Draw()

c.Print('graphics/qRecoilPt.pdf(')
c.SetLogy()
c.Print('graphics/qRecoilPt.pdf)')

