import ROOT as r
r.gROOT.SetBatch(True)
r.gROOT.ProcessLine(".L tdrstyle.C")
r.setTDRStyle()

from inputs import channel_data

d_xs_wj = 0.7822232858137581
d_xs_tt = 0.15925221429457714
factor_elqcd = 3.1447404463336044
factor_muqcd = 1.2367604160238135


channels = dict([(('_'.join([lep,par])), 
                  channel_data(lep, par, signal='fitTopQueuedBin7TridiscriminantWTopQCD')) 
                 for lep in ['el', 'mu'] for par in ['top','QCD']])

for n,c in channels.items():
    lep,par = n.split('_')
    c.samples.update(channel_data(lep, par, signal='fitTopQueuedBin7TridiscriminantWTopQCD', getTT=True).samples)
    for item in ['gg','qq','qg','ag'] : del c.samples['tt'+item]
    c.samples['tt'].eff /= 4
    c.samples['tt'].xs *= (1+d_xs_tt)
    c.samples['wj'].xs *= (1+d_xs_wj)
    print c



#N = lumi*xs*eff*(1+delta)

#channel
#xs eff

tfile = {'mu_QCD':r.TFile.Open('data/control_QCD_mu_ph_sn_jn_20.root'),
         'mu_top':r.TFile.Open('data/control_top_mu_ph_sn_jn_20.root'),
         'el_QCD':r.TFile.Open('data/control_QCD_el_ph_sn_jn_20.root'),
         'el_top':r.TFile.Open('data/control_top_el_ph_sn_jn_20.root')}

names = set(k.GetName().replace('el','%(lep)s').replace('mu','%(lep)s') for f in tfile.values() for k in f.GetListOfKeys() if 'Moment' not in k.GetName())

def hists(key,name):
    hs = {}
    ch = channels[key]
    f = tfile[key]
    d = name
    for s,S in ch.samples.items():
        if s=='data': continue
        h = f.Get(d+'/'+s)
        n = ch.lumi * S.xs * S.eff
        h.Scale(n/h.Integral(0,h.GetNbinsX()+1))
        hs[s]=h
    data = f.Get(d+'/'+'data')
    data.SetMinimum(0)
    return data,hs

def make(gname):
    fname = 'graphics/'+ (gname%{'lep':'lep'}).replace('[','').replace(']','') + '.pdf'
    c = r.TCanvas()
    c.Print(fname+'[')
    for lep in ['el','mu']:
        name = gname%{'lep':lep}
        
        sig = hists('%s_top'%lep, name)
        bkg = hists('%s_QCD'%lep, name)
        bkg[1]['dy'] = None
        mj = bkg[0].Clone('mj')
        sigStack = r.THStack('sigstack','')
        bkgStack = r.THStack('bkgstack','')
        for i,n in enumerate(sig[1]):
            s = sig[1][n]
            b = bkg[1][n]
            s.SetFillColor(i+2)
            s.SetLineColor(i+2)
            sigStack.Add(s)
            if b:
                b.SetFillColor(i+2)
                b.SetLineColor(i+2)
                bkgStack.Add(b)
                mj.Add(b,-1)
        mj.Scale(eval('factor_%sqcd'%lep))
        mj.SetFillColor(r.kGreen)
        mj.SetLineColor(r.kGreen)
        sigStack.Add(mj)
        bkg[0].Draw()
        bkgStack.Draw('hist same')
        bkg[0].Draw('same')
        c.Print(fname)
        sig[0].Draw()
        sigStack.Draw('hist same')
        sig[0].Draw('same')
        c.Print(fname)
    c.Print(fname+']')

for name in names: make(name)
