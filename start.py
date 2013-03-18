import sys,math,model,roo,ROOT as r
from falphaLSlice import falphaLSlice
r.gROOT.SetBatch(1)


class topAsymmFit(object) :
    @roo.quiet
    def __init__(self, dist, provar, tag ) :
        self.model = model.topModel( dist = dist )
        self.import_data(self.model.w)
        for item in ['falphaL','d_lumi','d_xs_dy','d_xs_st'] : self.model.w.arg(item).setConstant()

        if 'QueuedBin' not in dist :
            for item in ['alphaL','alphaT'] :
                self.model.w.arg(item).setConstant()
                self.model.w.arg(item).setVal(1)

        print '\n'.join(str(i) for i in ['',self.model.channels['el'],'',self.model.channels['mu'],''])
        #self.plot_fracs(self.model.w)
        self.defaults(self.model.w)
        self.print_fracs(self.model.w)
        self.print_n(self.model.w)

        self.fitArgs = [r.RooFit.Extended(True),
                        r.RooFit.ExternalConstraints(self.model.w.argSet('constraints')),
                        r.RooFit.Constrain(self.model.w.argSet('d_lumi,d_xs_st,d_xs_dy,d_xs_tt,d_xs_wj')),
                        r.RooFit.NumCPU(4),
                        r.RooFit.PrintLevel(-1),
                        ]
        central = self.model.w.pdf('model').fitTo( self.model.w.data('data'), *(self.fitArgs+[r.RooFit.Save(True)]))
        central.Print()

        self.print_fracs(self.model.w)
        self.print_n(self.model.w)
        print
    
        with open('data/falphaL_scan_%s_%s.txt'%(dist,provar),'w') as output :
            slices,lo,hi = 40,-0.2,0.4
            print >> output, '#'+'\t'.join(falphaLSlice.columns())
            for i in range(slices+1) :
                print 'slice %d'%i
                falphaL = hi - i*(hi-lo)/slices
                self.model.w.arg('falphaL').setVal(falphaL)
                sl = falphaLSlice(  self.model.w, self.fitArgs , 'falphaT' )
                print >> output, str(sl)
                output.flush()

    def print_fracs(self,w) :
        for item in ['lumi_mu','lumi_el','f_gg','f_qg','f_qq','f_ag']+['xs_'+i for i in self.model.channels['el'].samples if i!='data'] : print "%s: %.04f"%(item, w.arg(item).getVal())
        print

    def print_n(self,w) :
        length = 24
        tots = {'el':0,'mu':0}
        print (' ').join(i.rjust(8) for i in ['']+tots.keys())
        for xs in ['tt','wj','mj','st','dy'] :
            if xs=='data' : continue
            print xs.rjust(length/3),
            for chan in tots :
                val = w.arg('expect_%s_%s'%(chan,xs)).getVal()
                tots[chan]+=val
                print ("%d"%val).rjust(length/3),
            print
        print '-'*(3+length)
        print ' '.join(["tot".ljust(length/3)] +
                       [("%d"%t).rjust(length/3) for chan,t in tots.items()])
        print


    def defaults(self,w ) :
        w.var('d_qq').setVal(0)
        w.var('R_ag').setVal(w.arg('f_ag_hat').getVal()/w.arg('f_qq_hat').getVal())

    @roo.quiet
    def import_data(self,w) :
        obs_ = w.argSet(','.join(self.model.observables))
        obs = r.RooArgList(obs_)

        datas = [(lepton, r.RooDataHist('data_'+lepton,'N_{obs}^{%s}'%lepton,obs,chan.samples['data'].datas[0])) for lepton,chan in self.model.channels.items()]

        [roo.wimport(w, dat) for _,dat in datas]
        args = [r.RooFit.Import(*dat) for dat in datas]
        roo.wimport(w, r.RooDataHist('data', 'N_{obs}', 
                                     obs,
                                     r.RooFit.Index(w.arg('channel')),
                                     *args
                                     ),
                )


    def plot_fracs(self,w, logy = False) :
        r.gErrorIgnoreLevel = r.kWarning
        c = r.TCanvas()
        c.Print('fractions.pdf[')
        c.SetLogy(logy)
        plt = w.var('d_qq').frame()
        plt.SetMaximum(1.0)
        for v in range(0,94) :
            val = v/100.+0.07
            plt.SetTitle("R_ag = %.3f"%(val))
            w.var('R_ag').setVal(val)
            w.arg('f_gg').plotOn(plt,r.RooFit.LineWidth(1),r.RooFit.LineColor(r.kBlue))
            w.arg('f_qg').plotOn(plt,r.RooFit.LineWidth(1),r.RooFit.LineColor(r.kRed))
            w.arg('f_qq').plotOn(plt,r.RooFit.LineWidth(1),r.RooFit.LineColor(r.kGreen))
            w.arg('f_ag').plotOn(plt,r.RooFit.LineWidth(1),r.RooFit.LineColor(r.kViolet))
            plt.Draw()
            c.Print('fractions%s.pdf'%('_logy' if logy else ''), 'pdf')
        r.gErrorIgnoreLevel = r.kInfo
        c.Print('fractions.pdf]')

if __name__=='__main__' :
    distributions = ['fitTopQueuedBin7TridiscriminantWTopQCD',
                     'fitTopPtOverSumPt_triD',
                     'fitTopTanhRapiditySum_triD',
                     'fitTopTanhAvgRapidity_triD',
                     ]
    dist = distributions[0]
    variables = ['alphaL','R_ag']
    provar = variables[0]

    topAsymmFit( dist, provar, tag=None)
