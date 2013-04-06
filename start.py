import sys,math,model,roo,numpy as np,ROOT as r
from falphaLSlice import falphaLSlice
from paraboloid import paraboloid
r.gROOT.SetBatch(1)


class topAsymmFit(object) :
    @roo.quiet
    def __init__(self, dist, provar, tag ) :
        self.model = model.topModel( dist = dist , asymmetry = 'QueuedBin' in dist)
        self.import_data(self.model.w)

        print '\n'.join(str(i) for i in ['',self.model.channels['el'],'',self.model.channels['mu'],''])
        #self.plot_fracs(self.model.w)
        self.print_fracs(self.model.w)
        self.print_n(self.model.w)

        self.fitArgs = [r.RooFit.Extended(True),
                        r.RooFit.NumCPU(4),
                        r.RooFit.PrintLevel(-1),
                        ]
        for i in range(5):
            central = self.model.w.pdf('model').fitTo( self.model.w.data('data'), *(self.fitArgs+[r.RooFit.Save(True)]))
        central.Print()

        self.print_fracs(self.model.w)
        self.print_n(self.model.w)
        print

        faL,faLe = self.model.w.arg('falphaL').getVal(), self.model.w.arg('falphaL').getError()
        faT,faTe = self.model.w.arg('falphaT').getVal(), self.model.w.arg('falphaT').getError()

        nll = self.model.w.pdf('model').createNLL( self.model.w.data('data'), *self.fitArgs[:-1] )
        pll = nll.createProfile(self.model.w.argSet('falphaL,falphaT'))

        def point(faL_,faT_) :
            self.model.w.arg('falphaL').setVal(faL+dfaL)
            self.model.w.arg('falphaT').setVal(faT+dfaT)
            return faL_,faT_,pll.getVal()

        points = [point(faL+dfaL,faT+dfaT) for dfaL,dfaT in [(0,0),(-faLe,0),(faLe,0),(0,-faTe),(0,faTe),(faLe/2,faTe/2)]]
        parb = paraboloid(points)
        oneSigmaNLL = 1.14
        twoSigmaNLL = 3.0
        oneSigmas = parb.dxy(oneSigmaNLL)
        print oneSigmas
        print math.sqrt(oneSigmas[0,0]), math.sqrt(oneSigmas[1,1])
        param1sigma = parb.parametricEllipse(oneSigmaNLL)
        param2sigma = parb.parametricEllipse(twoSigmaNLL)
        scales = [self.model.w.arg(a).getVal() for a in ['Ac_y_ttqq','Ac_y_ttqg']]
        with open('points.txt','w') as wfile :
            for t in np.arange(0,2*math.pi+0.00001,math.pi/50) :
                point1 = param1sigma.dot([math.cos(t),math.sin(t),1])
                point2 = param2sigma.dot([math.cos(t),math.sin(t),1])
                point1/=point1[2]
                point2/=point2[2]
                seq = list(parb.xymin) + list(point1[:2]) + list(point2[:2]) + scales
                print>>wfile, '\t'.join(str(f) for f in seq)
        print "Wrote points.txt"

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
                     'TridiscriminantQQggQg_triD'
                     ]
    if len(sys.argv)>1 :
        iDist = sys.argv[1]
    else:
        print '\n'.join('%d %s'%i for i in enumerate(distributions))
        iDist = raw_input("Which?")
    dist = distributions[int(iDist)] if int(iDist) in range(len(distributions)) else ''
    if not dist :
        print 'not a number'
        exit(0)
    else: print "fitting",dist

    variables = ['alphaL','R_ag']
    provar = variables[0]

    topAsymmFit( dist, provar, tag=None)
