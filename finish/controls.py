import math
import lib
import ROOT as r
r.gROOT.SetBatch(True)
r.gROOT.ProcessLine(".L tdrstyle.C")
r.setTDRStyle()
r.tdrStyle.SetErrorX(r.TStyle().GetErrorX())
r.tdrStyle.SetPadTopMargin(0.065)
r.TGaxis.SetMaxDigits(3)
#r.tdrStyle.SetPadRightMargin(0.06)
from inputs import channel_data

d_xs_wj = 0.758705 #0.75986
d_xs_tt = 0.165759 #0.166398
factor_qcd = {'el':3.1996, #3.20075,
              'mu':1.25478 #1.25533
}

channels = dict([(('_'.join([lep,par])), 
                  channel_data(lep, par, signal='fitTopQueuedBin5_TridiscriminantWTopQCD')) 
                 for lep in ['el', 'mu'] for par in ['top','QCD']])

for n,c in channels.items():
    lep,par = n.split('_')
    c.samples.update(channel_data(lep, par, signal='fitTopQueuedBin5_TridiscriminantWTopQCD', getTT=True).samples)
    for item in ['gg','qq','qg','ag'] : del c.samples['tt'+item]
    c.samples['tt'].eff /= 4
    c.samples['tt'].xs *= (1+d_xs_tt)
    c.samples['wj'].xs *= (1+d_xs_wj)
    print c


tfile = {'mu_QCD':r.TFile.Open('data/control_QCD_mu_ph_sn_jn_20.root'),
         'mu_top':r.TFile.Open('data/control_top_mu_ph_sn_jn_20.root'),
         'el_QCD':r.TFile.Open('data/control_QCD_el_ph_sn_jn_20.root'),
         'el_top':r.TFile.Open('data/control_top_el_ph_sn_jn_20.root')}

names = set(k.GetName().replace('el','%(lep)s').replace('mu','%(lep)s') for f in tfile.values() for k in f.GetListOfKeys() if 'Moment' not in k.GetName())
rebins = dict([(n,1) for n in names])
rebins['LiCSV']=4
rebins['fitTopTanhRapiditySum']=2

colors = {'tt':r.kBlue, 'wj':r.kGreen, 'mj':r.kRed, 'st':r.kGray, 'dy':r.kGray}

zero = r.TF1('zero','0',-10000,10000)

labels = {'chia':'#chi_{a}', 
          'lepMetMt':'M_{T} (GeV)',
          'ProbabilityHTopMasses':'P_{MSD}', 
          'MSD':'MSD',
          'TopRatherThanWProbability':'P_{CSV}', 
          'LiCSV':'L_{i}^{CSV} / max(L^{CSV})',
          'fitTopTanhRapiditySum':'tanh|y_{t#bar{t}}|',
          'fitTopSumP4.mass':'m_{t#bar{t}} (GeV)',
          'rawHadTopMass':'m_{bpq} (GeV)',
          'rawHadWMass':'m_{pq} (GeV)', 
          'jetPti0':'jet0 p_{T} (GeV)',
          'jetPti1':'jet1 p_{T} (GeV)',
          'jetPti2':'jet2 p_{T} (GeV)',
          'jetPti3':'jet3 p_{T} (GeV)',
          'jetAdjustedP4.absEtai0': '#eta jet0', 
          'jetAdjustedP4.absEtai1': '#eta jet1', 
          'jetAdjustedP4.absEtai2': '#eta jet2', 
          'jetAdjustedP4.absEtai3': '#eta jet3',
          'metAdjustedP4.pt':'#slash{E}_{T} (GeV)', 
          'lepP4.pti0':'%(lep)s p_{T} (GeV)',
          'lepP4.absEtai0':'%(lep)s #eta',
          'TridiscriminantWTopQCD':'#Delta',
          'residualCDF_lepLfitTopRecoIndex':'cdf(residual): MET1',
          'residualCDF_lepSfitTopRecoIndex':'cdf(residual): MET2',
          'residualCDF_lepBfitTopRecoIndex':'cdf(residual): leptonic b jet',
          'residualCDF_lepTfitTopRecoIndex':'cdf(residual): leptonic top mass',
          'residualCDF_hadTfitTopRecoIndex':'cdf(residual): hadronic top mass',
          'residualCDF_hadPfitTopRecoIndex':'cdf(residual): light jet 1',
          'residualCDF_hadQfitTopRecoIndex':'cdf(residual): light jet 2',
          'residualCDF_hadBfitTopRecoIndex':'cdf(residual): hadronic b jet',
          'residualCDF_hadWfitTopRecoIndex':'cdf(residual): hadronic W mass',
          'residualCDF_lepWfitTopRecoIndex':'cdf(residual): leptonic W mass',
}

def flows(h):
    bins = h.GetNbinsX()
    lib.combineBinContentAndError(h,binToContainCombo=1,binToBeKilled=0)
    lib.combineBinContentAndError(h,binToContainCombo=bins,binToBeKilled=bins+1)

def hists(key,name,rbin=1):
    hs = {}
    ch = channels[key]
    f = tfile[key]
    d = name
    for s,S in ch.samples.items():
        if s=='data': continue
        h = f.Get(d+'/'+s)
        n = ch.lumi * S.xs * S.eff
        h.Scale(n/h.Integral(0,h.GetNbinsX()+1))
        h.Rebin(rbin)
        flows(h)
        hs[s]=h
    data = f.Get(d+'/'+'data')
    data.SetMinimum(0)
    data.Rebin(rbin)
    flows(data)
    return data,hs

def getRatio(data,hists):
    den = data.Clone()
    den.Reset()
    for h in filter(None,hists.values()) : den.Add(h)
    ratio = lib.ratioHistogram(data,den,0.08)
    ratio.SetTitle(';'+data.GetXaxis().GetTitle()+';'+'Ratio')
    ratio.SetMarkerSize(0)
    ratio.SetMinimum(0.7)
    ratio.SetMaximum(1.3)
    ratio.SetLabelSize(0.21)
    ratio.SetLabelOffset(0.04)
    ratio.GetXaxis().SetTitleSize(0.25)
    ratio.GetXaxis().SetTitleOffset(0.92)
    ratio.SetTickLength(0.1)
    ratio.SetTickLength(0.01,'Y')
    ratio.GetYaxis().SetLabelSize(0.1)
    ratio.GetYaxis().SetTitleSize(0.1)
    ratio.GetYaxis().SetTitleOffset(0.4)
    return ratio

def logH(h):
    bins = [h.GetBinContent(i) for i in range(h.GetNbinsX()+2)]
    errs = [h.GetBinError(i) for i in range(h.GetNbinsX()+2)]
    logh = h.Clone('log')
    for i,(b,e) in enumerate(zip(bins,errs)):
        if not b: continue
        logh.SetBinContent(i,math.log(b))
        logh.SetBinError(i, 0.5 * (math.log(b+e)-math.log(max(0.5,b-e))))
    w = 0.5* (math.log(h.GetMaximum())-math.log(h.GetMinimum()))
    logh.SetMinimum(-w)
    logh.SetMaximum(w)
    logh.GetYaxis().SetTitle('log '+h.GetYaxis().GetTitle()+'  ')
    return logh

def cprep(c,ratio=False):
    c.Clear()
    if not ratio: c.Divide(1,1)
    else:
        split = 0.23
        c.Divide(1,2)
        bo = 0.01
        to = 0.975
        le = 0.01
        ri = 1-0.01
        c.cd(1).SetBottomMargin(0)
        c.cd(2).SetTopMargin(0)
        c.cd(2).SetBottomMargin(0.56)
        c.cd(1).SetPad(le,bo+split,ri,to)
        c.cd(2).SetPad(le,bo,ri,bo+split)
    

def make(gname):
    r.tdrStyle.SetPadRightMargin(0.11 if 'mass' in gname else 0.06)
    r.tdrStyle.SetPadLeftMargin(0.14 if 'mass' in gname else 0.15)
    ggname = (gname%{'lep':'lep'}).replace('[','').replace(']','')
    if ggname not in labels: print 'No label for ', ggname; return
    fname = 'graphics/control/'+ ggname.replace('.','-') + '.pdf'
    c = r.TCanvas()

    c.Print(fname+'[')
    for lep in ['el','mu']:
        name = gname%{'lep':lep}
        
        sig = hists('%s_top'%lep, name, rbin=rebins[gname])
        bkg = hists('%s_QCD'%lep, name, rbin=rebins[gname])
        bkg[1]['dy'] = None
        mj = bkg[0].Clone('mj')
        for n,color in colors.items():
            if n=='mj':
                mj.SetFillColor(color)
                mj.SetLineColor(color)
                continue
            s = sig[1][n]
            b = bkg[1][n]
            s.SetFillColor(color)
            s.SetLineColor(color)
            if b:
                b.SetFillColor(color)
                b.SetLineColor(color)
                mj.Add(b,-1)
        mj.Scale(factor_qcd[lep])
        sig[1]['mj'] = mj
        bkg[1]['mj'] = None

        def draw(data,hists):
            doRatio = all(hists.values())
            cprep(c,doRatio)
            data.UseCurrentStyle()
            data.GetXaxis().SetTitle(labels[ggname]%{'lep':{'el':'electron','mu':'muon'}[lep]})
            if doRatio:
                c.cd(2)
                ratio = getRatio(data,hists)
                logr = logH(ratio)
                logr.Draw()
                zero.Draw('same')
                logr.Draw('same')
                data.SetLabelSize(0)
                data.GetXaxis().SetTitleSize(0)
                data.GetYaxis().SetTitleSize(0.079)
                data.GetYaxis().SetTitleOffset(0.95)
                data.GetYaxis().SetLabelSize(0.065)
            stack = r.THStack('stack','')
            for item in ['dy','st','mj','wj','tt']:
                if item in hists and hists[item]: stack.Add(hists[item])
            mx = max(i.GetMaximum() for i in [data,stack])
            for item in [data,stack]: item.SetMaximum(1.05*mx)
            c.cd(1)
            data.Draw()
            stack.Draw('hist same')
            data.Draw('same')
            c.Print(fname)
            return

        draw(*bkg)
        draw(*sig)
    c.Print(fname+']')

for name in names: 
    make(name)
