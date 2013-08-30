import math
from numpy import cumsum, ediff1d

class shortest(object):
    '''Functor for shortest interval containing some probability.

    Half the shortest interval containing 68% of the distribution is
    suggested as a figure of merit (@property FOM) based on the
    coincidence with RMS for a normal distribution.
    '''

    def __init__(self, H, ndigits=6):
        self.ndigits = ndigits
        self.edges = [H.GetBinLowEdge(i) for i in range(1, 2+H.GetNbinsX())]
        self.cum = cumsum([H.GetBinContent(i) for i in range(2+H.GetNbinsX())])
        self.cum /= self.cum[-1]
        self.slopes = ediff1d(self.cum)[:-1] / ediff1d(self.edges)

    def key(self, (iLo, iHi)):
        return (round(self.edges[iHi]-self.edges[iLo], self.ndigits),
                1 / (min(self.slopes[i] for i in [iLo, iHi-1]) or 1e-6))

    def __call__(self, frac):
        iHis = [next(i+j+1 for i, ci in enumerate(self.cum[j+1:]) if ci-cj > frac)
                for j, cj in enumerate(self.cum) if cj+frac < self.cum[-2]]

        iLo, iHi = min(enumerate(iHis), key=self.key)

        trim = ((self.cum[iHi] - self.cum[iLo] - frac) /
                min(self.slopes[i] for i in [iLo, iHi-1]))

        return self.edges[iHi] - self.edges[iLo] - trim

    @property
    def FOM(self): return 0.5*self(self.xFOM)
    xFOM = math.erf(1/math.sqrt(2))


if __name__ == '__main__':
    '''Compare shortest interval FOM for several distributions of unit RMS.'''

    import math
    import ROOT as r
    r.gROOT.SetBatch(1)
    funcs = {'uniform': lambda x: 1. if abs(x) < 0.5*math.sqrt(12) else 0,
             'triang': lambda x: max(0, math.sqrt(6) - abs(x)),
             'normal': lambda x: math.exp(-0.5*x**2),
             'laplace': lambda x: math.exp(-math.sqrt(2)*abs(x))}

    hists = dict([(h, r.TH1D(h, '', 300, -6, 6)) for h in funcs])
    [hists[n].SetBinContent(i, f(hists[n].GetBinCenter(i)))
     for n, f in funcs.items()
     for i in range(0, 2+hists[n].GetNbinsX())]

    shorts = dict([(n, shortest(h)) for n, h in hists.items()])
    keys = sorted(funcs, key=lambda f: -shorts[f].FOM)

    print '#', '\t'.join([' ']+keys)
    print '#', '\t'.join(['RMS'] + ["%.3f" % hists[k].GetRMS() for k in keys])
    print '#', '\t'.join(['FOM'] + ["%.3f" % shorts[k].FOM for k in keys])
    print

    c = r.TCanvas('','',800,800)
    for i, k in enumerate(keys[::-1]):
        h = hists[k]
        h.Scale(1./h.Integral(0, 2+h.GetNbinsX()))
        h.SetLineColor(i+1)
        h.Draw('hist' + ('same' if i else ''))
    c.Print('shortest.pdf')
