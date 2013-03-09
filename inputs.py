import ROOT as r,utils

# can symmetrize:
#  DY
#  QCD

class sample_data(object) :
    def __init__(self, signalDistributions, xs = None, sigma = None, selectionEfficiency = 1.0, preselectionFraction = 1.0 ) :
        self.xs = xs
        self.xs_sigma = sigma
        self.eff = selectionEfficiency
        self.frac = preselectionFraction
        self.datas = signalDistributions if all(signalDistributions) else ((signalDistributions[0],) + tuple(utils.symmAnti(signalDistributions[0])))
        norm = self.datas[0].Integral()
        for d in filter(None,self.datas) :
            #d.Scale(1./norm)
            d.RebinY(20)
        
        self.datasX = tuple(d.ProjectionX() if d else None for d in self.datas)
        self.datasY = tuple(d.ProjectionY() if d else None for d in self.datas)
        for d in self.datasX+self.datasY+self.datas : d.SetDirectory(0)

    def __str__(self) : return 'data' if not self.xs else ';  '.join([('xs: % 8.2f'%self.xs),
                                                                      ('eff: % .5f'%self.eff).rjust(7), 
                                                                      ('f: %.4f'%self.frac).rjust(8),
                                                                      ('d: %.4f'%self.xs_sigma if self.xs_sigma else 'd: None   ').rjust(8)])
    @property
    def key(self) : return ( self.xs, self.frac )

class channel_data(object) :
    __samples__ = ['wj','dy','st','ttgg','ttqg','ttqq','ttag','qcd','data']
    __xs_uncertainty__ = {'tt':1.0,'wj':2.0,'st':0.04,'dy':0.04}

    def __init__(self,lepton, filePattern="data/stats_melded_%s_ph_pn_sn_jn_20.root",
                 signal="fitTopQueuedBin7TridiscriminantWTopQCD", 
                 preselection="allweighted/" ) :
        tfile = r.TFile.Open(filePattern%lepton)

        self.lepton = lepton
        self.lumi = tfile.Get('lumiHisto/data').GetBinContent(1)
        self.lumi_sigma = 0.04
        self.samples = {}

        for s in self.__samples__:
            pre = tfile.Get(preselection+s)
            datas = ( tfile.Get(signal+     '/'+s).Clone(lepton+'_'+s),
                      tfile.Get(signal+'_symm/'+s).Clone(lepton+'_symm_'+s) if s[:2]=='tt' else None,
                      tfile.Get(signal+'_anti/'+s).Clone(lepton+'_anti_'+s) if s[:2]=='tt' else None)
            for d in filter(None,datas) : d.SetDirectory(0)
            xs = tfile.Get('xsHisto/'+s).GetBinContent(1) if s!='data' else None
            delta = self.__xs_uncertainty__[s[:2]] if s[:2] in self.__xs_uncertainty__ else None

            self.samples[s] = sample_data( datas, xs, delta,
                                           selectionEfficiency = (datas[0].Integral()/pre.Integral() if pre else 0),
                                           preselectionFraction = 1.0 if s[:2]!='tt' else pre.Integral()/tfile.Get(preselection+'tt').Integral()
                                           )
        tfile.Close()

    def __str__(self) :
        return '%s : %.2f/pb  (%.2f)\n'%(self.lepton,self.lumi,self.lumi_sigma)+'\n'.join( '%s:: %s'%(s.rjust(5),str(data)) for s,data in sorted(self.samples.items(), key = lambda (s,d):d.key, reverse=True) if data.xs)


if __name__=='__main__' :
    channels = dict([(lep,channel_data(lep)) for lep in ['el','mu']])
    print
    print channels['el']
    print
    print channels['mu']
