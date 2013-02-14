import random, ROOT as r, sys

samples = ['wj','dy','st','ttgg','ttqg','ttqq','ttag','qcd','data']

class sample_data(object) :
    def __init__(self, signalDistribution, xs = None, lumi = None, selectionEfficiency = 1.0, preselectionFraction = 1.0 ) :
        self.data = signalDistribution
        self.data.Scale(1./self.data.Integral())

        self.eff = selectionEfficiency
        self.lumi = lumi
        self.xs = xs
        self.frac = preselectionFraction

        assert (xs==None)^(lumi==None)
        
    def __str__(self) : return ';  '.join([('  xs: %8.2f'%self.xs if self.xs else 'lumi: %8.2f'%self.lumi),
                                           ('eff: %.4f'%self.eff).rjust(7), 
                                           ('f: %.4f'%self.frac).rjust(8)])
    @property
    def key(self) : return (self.lumi, self.xs, self.frac)

class channel_data(object) :
    def __init__(self,lepton, filePattern="data/stats_melded_%s_ph_c_20.root",
                 signal="fitTopQueuedBin7TridiscriminantWTopQCD/", 
                 preselection="allweighted/" ) :

        tfile = r.TFile.Open(filePattern%lepton)

        self.samples = {}
        for s in samples:
            pre = tfile.Get(preselection+s)
            data = tfile.Get(signal+s)
            xs = tfile.Get('xsHisto/'+s).GetBinContent(1) if s!='data' else None
            lumi = tfile.Get('lumiHisto/'+s).GetBinContent(1) if s=='data' else None

            self.samples[s] = sample_data( data, xs, lumi, 
                                           selectionEfficiency = (data.Integral()/pre.Integral() if pre else 0), 
                                           preselectionFraction = 1.0 if s[:2]!='tt' else pre.Integral()/tfile.Get(preselection+'tt').Integral()
                                           )
    def __str__(self) :
        return '\n'.join( '%s:: %s'%(s.rjust(5),str(data)) for s,data in sorted(self.samples.items(), key = lambda (s,d):d.key, reverse=True))
        
        
for lepton in ['el','mu'] : 
    print
    print lepton
    print channel_data(lepton)

sys.exit(0)

luminosity = {'el': (5100,None),
              'mu': (5096,0.04)} # 1/pb

xs = {'tt' : ( 149.600, 1.0), # (pb, %)
      'wj' : (1911.800, 2.0),
      'mj' : (   16.000, None),
      'st' : (  71.968, 0.04),
      'dy' : (2475.000, 0.04)}

def ttcomps() :
    tfile = r.TFile.Open('data/stats_top_electron_ph.root')
    num = tfile.Get(counts)
    comps = ['qq','ag','gg','qg']
    cnts = [num.Get('tt'+c).Integral() for c in comps]
    tfile.Close()
    return zip(comps, [c/sum(cnts) for c in cnts])

components = dict( [(item,[('',1.0)]) for item in ['wj','mj','st','dy']] +
                   [('tt', ttcomps())] )

def getEfficiencies(chan) :
    tfile = r.TFile.Open('data/stats_top_%s_ph.root'%chan)
    num = tfile.Get(signals)
    den = tfile.Get(counts)
    effs = dict([(samp,[(comp,
                         num.Get(samp+comp).Integral()/
                         den.Get(samp+comp).Integral())
                        for comp,_ in comps]) 
                 for samp,comps in components.items() 
                 if samp!="mj"])
    tfile.Close()
    effs['mj'] = [('', None )]
    return effs

efficiency = {'mu' : getEfficiencies("muon"),
              'el' : getEfficiencies("electron")}

def histogram(dist, chan, samp, comp = '', d=1) : 
    fileName = 'data/stats_%s_%s_ph.root'%('melded' if samp=='mj' else 'top', 
                                           {'mu':'muon','el':'electron'}[chan])
    f = r.TFile.Open(fileName)
    h = f.Get(signals).Get((samp+comp) if samp!='mj' else 'QCD.multijet' )
    gram = h.Clone() if d==2 else h.ProjectionX('_px_'+chan+samp) if dist=='d3' else h.ProjectionY('_py_'+chan+samp)
    gram.SetDirectory(0)
    f.Close()
    rebin = {'d3':1,'ptpt':4}
    #if d==1 : gram.Rebin(rebin[dist]) 
    if d==2:   gram.Rebin2D(rebin['d3'],rebin['ptpt'])
    return gram

