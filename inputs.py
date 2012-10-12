import random, ROOT as r

counts = "beamHaloCSCLooseHaloId"
signals = "triD_v_sqtsumptopt"

luminosity = {'el': (5100,None),
              'mu': (5096,0.04)} # 1/pb

xs = {'tt' : ( 149.600, 1.0), # (pb, %)
      'wj' : (1911.800, 2.0),
      'mj' : (   1.000, None),
      'st' : (  71.968, 0.04),
      'dy' : (2475.000, 0.04)}

components = dict( [(item,[('',1.0)]) for item in ['wj','mj','st','dy']] +
                   [('tt',[('qq',1.2359e-01),
                           ('ag',3.7104e-02),
                           ('gg',6.1509e-01),
                           ('qg',2.2421e-01)]
                     )] )

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
    gram = h.Clone() if d==2 else h.ProjectionX('_px_'+chan) if dist=='d3' else h.ProjectionY('_py_'+chan)
    gram.SetDirectory(0)
    f.Close()
    if d==2 : gram.Rebin2D(1,4)
    return gram

