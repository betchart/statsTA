import ROOT as r
import utils


class sample_data(object):
    def __init__(self, signalDistributions, xs=None, sigma=None,
                 selectionEfficiency=1.0, preselectionFraction=1.0):
        self.xs = xs
        self.xs_sigma = sigma
        self.eff = selectionEfficiency
        self.frac = preselectionFraction
        self.datas = (signalDistributions if all(signalDistributions) else
                      ((signalDistributions[0],) +
                       tuple(utils.symmAnti(signalDistributions[0]))))

        self.alphaMax = utils.alphaMax(*self.datas[1:])

        for d in filter(None, self.datas):
            if d.GetNbinsY() == 100: d.RebinY(20)
            if d.GetNbinsX() > 80: d.RebinX()

        self.datasX = tuple(d.ProjectionX() if d else None for d in self.datas)
        self.datasY = tuple(d.ProjectionY() if d else None for d in self.datas)
        for d in self.datasX + self.datasY + self.datas: d.SetDirectory(0)

    def subtract(self,other):
        assert self.xs == other.xs
        assert self.xs_sigma == other.xs_sigma
        assert self.frac == other.frac
        assert self.eff >= other.eff

        self.eff -= other.eff
        for group in ['datas','datasX','datasY']:
            for d,od in zip(getattr(self,group),getattr(other,group)):
                d.Add(od,-1)

    def __str__(self):
        return ('data' if not self.xs else
                ';  '.join([('xs: % 8.2f' % self.xs),
                            ('eff: % .5f' % self.eff).rjust(7),
                            ('f: %.4f' % self.frac).rjust(8),
                            ('d: %.4f' % self.xs_sigma if self.xs_sigma else
                             'd: None   ').rjust(8)]))

    @property
    def key(self): return (self.xs, self.frac)


class channel_data(object):
    __samples__ = ['data', 'wj', 'dy', 'st', 'ttgg', 'ttqg', 'ttqq', 'ttag', 'tt']
    __xs_uncertainty__ = {'tt': 1.0, 'wj': 2.0, 'st': 0.04, 'dy': 0.04}

    def __init__(self, lepton, partition, tag = 'ph_sn_jn_20',
                 signal="", sigPrefix="", dirPrefix="R02", getTT=False,
                 prePre = False, hackZeroBins=False):
        filePattern="data/stats_%s_%s_%s.root"
        tfile = r.TFile.Open(filePattern % (partition, lepton, tag))

        self.lepton = lepton
        self.lumi = tfile.Get('lumiHisto/data').GetBinContent(1)
        self.lumi_sigma = 0.04
        self.hackZeroBins = hackZeroBins

        def full(pf) :
            return next((ky.GetName() + '/' for ky in tfile.GetListOfKeys()
                         if pf == ky.GetName().split('_')[0]),
                        '')
        fullDirName = full(dirPrefix)

        paths = (fullDirName + sigPrefix + signal,
                 fullDirName + signal)

        prepaths = (full('R01') + (sigPrefix if prePre else '') + 'genTopDeltaBetazRel; genTopPhiBoost',
                    'meweighted/')

        self.samples = {}
        for s in self.__samples__[4 if getTT else 0:None if getTT else -1]:
            self.add(s, tfile, paths, prepaths)
        tfile.Close()

    def add(self, s, tfile, paths, prepaths):
        def get(s,ps):
            return next(iter(filter(None, [utils.get(tfile,p+s) for p in ps])), None)

        pre = get( s, prepaths)
        if not pre and not s == 'data': return
        doSymmAnti = s[:2] == 'tt' and 'QueuedBin' in paths[0]

        datas = (get('/' + s,paths).Clone(self.lepton + '_' + s),
                 get('_symm/' + s,paths).Clone(self.lepton + '_symm_' + s)
                 if doSymmAnti else None)
        if doSymmAnti:
            datas += (datas[0].Clone(self.lepton + '_anti_' + s),)
            datas[2].Add(datas[1], -1)

        for d in filter(None, datas): d.SetDirectory(0)
        xs = tfile.Get('xsHisto/' + s).GetBinContent(1) if s != 'data' else None
        delta = (self.__xs_uncertainty__[s[:2]]
                 if s[:2] in self.__xs_uncertainty__ else None)

        named = \
            {'selectionEfficiency': (datas[0].Integral() / pre.Integral() if pre else 0),
             'preselectionFraction': (1.0 if s[:2] != 'tt' else
                                      pre.Integral() / get( 'tt', prepaths ).Integral())
             }
        self.samples[s] = sample_data(datas, xs, delta, **named)
        # RooFit breaks with zero bins in data but not the pdf:
        #  this hack fixes one breaking case.
        if self.hackZeroBins and s=='data':
            for iX in range(1,1+self.samples[s].datas[0].GetNbinsX()):
                for iY in range(1,1+self.samples[s].datas[0].GetNbinsY()):
                    if not self.samples[s].datas[0].GetBinError(iX,iY):
                        self.samples[s].datas[0].SetBinError(iX,iY,1)
                        self.samples[s].datas[0].SetBinContent(iX,iY,1)

    def subtract(self, other):
        for s,samp in self.samples.items():
            samp.subtract(other.samples[s])

    def __str__(self):
        return ('%s : %.2f/pb  (%.2f)\n' % (self.lepton, self.lumi, self.lumi_sigma) +
                '\n'.join('%s:: %s' % (s.rjust(5), str(data))
                          for s, data in sorted(self.samples.items(),
                                                key=lambda (s, d): d.key,
                                                reverse=True)
                          if data.xs))


if __name__ == '__main__':
    channels = dict([(lep, channel_data(lep, 'top')) for lep in ['el', 'mu']])
    print
    print channels['el']
    print
    print channels['mu']
