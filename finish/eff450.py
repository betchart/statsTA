import sys
import math
import ROOT as r

if not len(sys.argv)>1:
    print 'Usage: eff450.py <data/statsfile.root>'
    exit()

tfile = r.TFile.Open(sys.argv[1])

d1 = tfile.Get('R04_genPdfWeights53_symmAnti_fitTopQueuedBin7inantisymmpartsoftt_TridiscriminantWTopQCD/fitTopQueuedBin7TridiscriminantWTopQCD')
d2 = tfile.Get('R16_genPdfWeights53_symmAnti_fitTopQueuedBin7inantisymmpartsoftt_TridiscriminantWTopQCD_mass_450.00<=fitTopSumP4.mass/fitTopQueuedBin7TridiscriminantWTopQCD')

samples = ['data', 'tt', 'wj', 'st', 'dy']

def eff(name):
    one = d1.Get(name)
    two = d2.Get(name)
    n = one.GetEffectiveEntries()
    N,P = [h.Integral() for h in [one,two]]
    p = P/N
    q = 1-p
    return 100*p, 100*math.sqrt(p*q/n)

space = "   &   "
print r'$m_{\ttbar}>450\GeV$', space,  d2.Get('data').Integral()/1000, space, space.join('%.4f(%.4f)'%eff(s) for s in samples[1:-1 if 'QCD' in sys.argv[1] else None])

