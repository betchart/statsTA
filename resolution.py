import utils
import ROOT as r
r.gROOT.SetBatch(True)
from shortest import shortest

r.gROOT.ProcessLine(".L tdrstyle.C")
r.setTDRStyle()
r.tdrStyle.SetErrorX(0.5);

filenames = {'qq':'tmp/%s/ttj_ph.wQQ.pu.sf_plots.root',
             'gg':'tmp/%s/ttj_ph.wGG.pu.sf_plots.root',
             'qg':'tmp/%s/ttj_ph.wQG.pu.sf_plots.root',
             'ag':'tmp/%s/ttj_ph.wAG.pu.sf_plots.root',
         }

tfiles_el = dict([(k,r.TFile.Open(v%'el')) for k,v in filenames.items() ])
tfiles_mu = dict([(k,r.TFile.Open(v%'el')) for k,v in filenames.items() ])

c = r.TCanvas()
#c.SetGrid()
outname = 'graphics/res.pdf'

def combined(n,(path_el,path_mu),var):
    el = tfiles_el[n].Get(path_el+var)
    mu = tfiles_el[n].Get(path_mu+var)
    el.Scale(el.GetEffectiveEntries()/el.Integral())
    mu.Scale(mu.GetEffectiveEntries()/mu.Integral())
    comb = el.Clone('combined')
    comb.Add(mu)
    return comb


def do2(paths,var):
    colors = [r.kRed,r.kBlack,r.kBlue,r.kGreen]
    rmsFOMs = {}
    for i,n in enumerate(tfiles_el):
        two = '/label/kinFitLook/resolutions'
        fit = combined(n,paths,var+'_reco')
        fit2 = combined(n,paths,two+var+'_reco')
        un = combined(n,paths,var+'_unfit')
        un2 = combined(n,paths,two+var+'_unfit')
        rmsFOMs[n] = dict([(item, (eval(item).GetRMS(),shortest(eval(item)).FOM)) for item in ['fit','un','fit2','un2']])
        for h in [fit,un,fit2,un2]: 
            h.Scale(1./h.Integral(0,2+h.GetNbinsX()))
            h.Rebin()
            h.SetMaximum(0.20 if 'XL' in var else 0.12)
            h.UseCurrentStyle()
            h.SetMarkerSize(0)
            h.SetMarkerStyle(0)
            h.GetYaxis().SetTitle('probability')
            h.GetXaxis().SetTitle(h.GetXaxis().GetTitle().replace('XT','X_{T}').replace('XL','X_{L}'))

        diff_plus = fit.Clone('diff_plus')
        diff_minus = fit.Clone('diff_minus')
        diff_plus2 = fit.Clone('diff_plus2')
        diff_minus2 = fit.Clone('diff_minus2')
        diff_plus.Reset()
        diff_minus.Reset()
        diff_plus2.Reset()
        diff_minus2.Reset()
        for i in range(0,2+fit.GetNbinsX()):
            a = (fit.GetBinContent(i) + un.GetBinContent(i))/2
            d = (fit.GetBinContent(i) - un.GetBinContent(i))/2
            a2 = (fit2.GetBinContent(i) + un2.GetBinContent(i))/2
            d2 = (fit2.GetBinContent(i) - un2.GetBinContent(i))/2
            if d>0:
                diff_plus.SetBinContent(i,a)
                diff_plus.SetBinError(i,d)
                diff_plus2.SetBinContent(i,a2)
                diff_plus2.SetBinError(i,d2)
            else: 
                diff_minus.SetBinContent(i,fit.GetBinContent(i))
                diff_minus2.SetBinContent(i,fit2.GetBinContent(i))
        s = "%.3f / %.3f"
        print n, s%(shortest(fit).FOM,shortest(un).FOM)
        fit.SetTitle(n)
        diff_plus.SetFillColor(r.kBlack)
        diff_minus.SetFillColor(r.kBlack)
        diff_plus2.SetFillColor(r.kRed)
        diff_minus2.SetFillColor(r.kRed)
        for h in [fit2,un2,diff_plus2,diff_minus2]: 
            h.SetLineColor(r.kRed)
            h.SetMarkerColor(r.kRed)
        diff_plus.Draw('')
        diff_plus.Draw('E2 same')
        diff_minus2.Draw('hist same')
        diff_minus.Draw('hist same')
        diff_plus2.Draw('E2 same')
        fit2.Draw('hist same')
        un2.Draw('hist same')
        fit.Draw('hist same')
        un.Draw('hist same')
        initialstate = {'qq':"q#bar{q}",'gg':'gg','qg':'qg','ag':'#bar{q}g'}[n]
        label = r.TLatex(-1.7,0.18,"%s#rightarrow{}^{}t#bar{t}"%initialstate)
        label.Draw()
        leg = r.TLegend(0.65,0.65,0.95,0.9)
        print rmsFOMs[n]
        leg.AddEntry("","RMS | HSI^{68}","")
        [leg.AddEntry(eval(item)," %.2f | %.2f"%tuple([round(f,2) for f in rmsFOMs[n][item.replace('diff_plus','fit')]]), 'F') for i,item in enumerate(['diff_plus','un','diff_plus2','un2'])]
        leg.SetBorderSize(0)
        leg.SetFillColor(r.kWhite)
        leg.Draw()
        c.Print(outname)
    return rmsFOMs

def do(paths):
    return {'XL':do2(paths,'/XL'),
            'XT':do2(paths,'/XT')}

path_el = tfiles_el['qq'].FindObjectAny("XL_reco").GetDirectory().GetPath().split(":")[1]
path_mu = tfiles_mu['qq'].FindObjectAny("XL_reco").GetDirectory().GetPath().split(":")[1]
c.Print(outname+'[')
foms = do((path_el,path_mu))
c.Print(outname+']')
