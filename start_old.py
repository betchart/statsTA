import sys,math,model,roo,ROOT as r
from dqqSlice import dqqSlice
r.gROOT.SetBatch(1)


class topAsymmFit(object) :
    @roo.quiet
    def __init__(self, dist, provar, tag ) :
        self.model = model.topModel( dist = dist )
        self.import_data(self.model.w)
        for item in ['d_qq','d_lumi','d_xs_dy','d_xs_st'] : self.model.w.arg(item).setConstant()

        if 'QueuedBin' not in dist :
            for item in ['alphaL','alphaT'] :
                self.model.w.arg(item).setConstant()
                self.model.w.arg(item).setVal(1)

        print '\n'.join(str(i) for i in ['',self.model.channels['el'],'',self.model.channels['mu'],''])
        #self.plot_fracs(self.model.w)
        self.defaults(self.model.w)
        self.print_fracs(self.model.w)
        self.print_n(self.model.w)
        #self.draw(w,'before')

        self.fitArgs = [r.RooFit.Extended(True),
                        r.RooFit.ExternalConstraints(self.model.w.argSet('constraints')),
                        r.RooFit.Constrain(self.model.w.argSet('d_lumi,d_xs_st,d_xs_dy,d_xs_tt,d_xs_wj')),
                        #r.RooFit.Optimize(False),
                        r.RooFit.NumCPU(4),
                        r.RooFit.PrintLevel(-1),
                        ]
        self.model.w.arg('d_qq').setVal(0)
        #self.model.w.arg('R_ag').setConstant()
        #self.model.w.arg('alphaL').setConstant()
        #self.model.w.arg('alphaT').setConstant()
        central = self.model.w.pdf('model').fitTo( self.model.w.data('data'), *(self.fitArgs+[r.RooFit.Save(True)]))
        central.Print()

        self.print_fracs(self.model.w)
        self.print_n(self.model.w)
        print
        print
        self.model.w.Print()
        return
    
        with open('data/dqq_scan_%s_%s.txt'%(dist,provar),'w') as output :
            slices,lo,hi = 40,-0.99,0.4
            print >> output, '#'+'\t'.join(dqqSlice.columns())
            for i in range(slices+1) :
                print 'slice %d'%i
                dqq = hi - i*(hi-lo)/slices
                self.model.w.arg('d_qq').setVal(dqq)
                #sl = dqqSlice(  self.model.w, self.fitArgs , 'alphaL' )
                sl = dqqSlice(  self.model.w, self.fitArgs , provar )
                print >> output, str(sl)
                output.flush()
                if False :
                #if dqq==0 :
                    self.print_fracs(self.model.w)
                    self.print_n(self.model.w)
                    c = r.TCanvas()
                    c.Divide(2,2)
                    fileName = 'graphics/ensembleTest%03d.pdf'%(1000*dqq)
                    c.Print(fileName+'[')
                    pullMean = {}
                    pullSigma = {}
                    for alphaL in [-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3] :
                        fit = self.model.w.pdf('model').fitTo( self.model.w.data('data'), *(self.fitArgs+[r.RooFit.Save(True)]))
                        self.model.w.arg('alphaL').setVal(alphaL)
                        mcstudy = r.RooMCStudy( self.model.w.pdf('model'), self.model.w.argSet(','.join(self.model.observables+['channel'])),
                                                r.RooFit.Binned(True), r.RooFit.Extended(True), r.RooFit.FitOptions(*self.fitArgs) )

                        mcstudy.generateAndFit(1000)
                        bins = r.RooFit.Bins(40)
                        frame4 = mcstudy.plotNLL(bins)
                        c.cd(4)
                        frame4.Draw()
                        arg = self.model.w.arg('alphaL')
                        for i,name in enumerate(['Param','Error','Pull']) :
                            frame = getattr(mcstudy, "plot"+name)(*([arg,bins]+([r.RooFit.FitGauss(True)] if name=='Pull' else [])))
                            c.cd(i+1)
                            frame.Draw()
                            #if name=='Pull':
                            # where is the f#$%ing programatic access to the gaussian fit to the pull distribution?
                        c.Print(fileName)
                    #
                    c.Print(fileName+']')

        #self.defaults(self.model.w)
        #self.print_fracs(self.model.w)
        #self.print_n(self.model.w)

        #mcstudy = r.RooMCStudy( self.model.w.pdf('model'), self.model.w.argSet(','.join(self.model.observables+['channel'])),
        #                        r.RooFit.Binned(True), r.RooFit.Extended(True), r.RooFit.FitOptions(*self.fitArgs) )
        #
        #
        #self.print_fracs(self.model.w)
        #self.print_n(self.model.w)

    def plotMCStudy(self,mcstudy) :
        c = r.TCanvas()
        c.Divide(2,2)
        c.Print('graphics/plots.pdf[')
        bins = r.RooFit.Bins(40)
        frame4 = mcstudy.plotNLL(bins)
        c.cd(4)
        frame4.Draw()
        for var in ['alphaL','alphaT','R_ag','d_lumi']+['d_xs_%s'%s for s in ['tt','wj','dy','st']]+['eff_%s_qcd'%l for l in ['el','mu']] :
            arg = self.model.w.arg(var)
            for i,name in enumerate(['Param','Error','Pull']) :
                frame = getattr(mcstudy, "plot"+name)(*([arg,bins]+([r.RooFit.FitGauss(True)] if name=='Pull' else [])))
                c.cd(i+1)
                frame.Draw()
            c.Print('graphics/plots.pdf')
        c.Print('graphics/plots.pdf]')
        del c


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

    def test(self, w) :
        if not w.arg('fsum') :
            factory(w, "sum::fsum(f_qq,f_ag,f_qg,f_gg)")
            factory(w, "expr::fprod('@0/@1 * @2/@3 - @4/@5 * @6/@7',{f_gg,f_gg_hat,f_qq,f_qq_hat,f_ag,f_ag_hat,f_qg,f_qg_hat})")
        print w.arg('fsum').getVal(), w.arg('fprod').getVal()


    @roo.quiet
    def draw(self, w, fileName) :
        if not hasattr(self,'maxi') : self.maxi = {}
        leg = r.TLegend(0.7,0.5,0.99,0.99)
        c = r.TCanvas()
        cats = self.model.channels.keys()
        data = w.data('data')
        fakedata = w.pdf('model').generate(r.RooArgSet(w.argSet(','.join(self.observables+['channel']))),
                                           r.RooFit.AllBinned() )
        c.Divide(2,len(cats))
        for i,cat in enumerate(cats) :
            for j,var in enumerate(self.observables) :
                pad = c.cd(i+1+2*j)
                plt = w.var(var).frame()
                cut = r.RooFit.Cut("channel==channel::%s"%cat)
                comps = [('dy','_symm'),('st','_both'),('qcd','_symm'),('wj','_both'),('ttgg','_both'),('ttag',''),('ttqg',''),('ttqq','')]
                colors = [r.kGray, r.kGray, r.kGreen, r.kRed, r.kBlue-6, r.kBlue+2, r.kBlue+1, r.kViolet]

                data.plotOn(plt, cut, r.RooFit.MarkerSize(0.5), r.RooFit.Name('data'+cat+var))
                for iComps in range(len(comps)) :
                    print '.',
                    col = colors[iComps]
                    compstr = ','.join(['%s_%s%s'%((cat,)+sub) for sub in comps[:iComps+1]])
                    w.pdf('model').plotOn( plt,
                                           r.RooFit.ProjWData(w.argSet(''),fakedata),
                                           r.RooFit.Slice(w.cat('channel'),cat),
                                           r.RooFit.Components(compstr),
                                           r.RooFit.LineWidth(0),
                                           r.RooFit.LineColor(col),
                                           r.RooFit.FillColor(col),
                                           r.RooFit.DrawOption('F'),
                                           r.RooFit.Name(comps[iComps][0]+cat+var),
                                           r.RooFit.MoveToBack()
                                           )
                args = w.argSet(','.join(['d_qq','R_ag','alphaL'
                                          ]))
                w.pdf('model').plotOn( plt,
                                       r.RooFit.ProjWData(w.argSet(''),fakedata),
                                       r.RooFit.Slice(w.cat('channel'),cat),
                                       r.RooFit.LineWidth(1),
                                       r.RooFit.LineColor(r.kOrange),
                                       r.RooFit.Name('model'+cat+var)
                                       )
                if j : data.statOn(plt, r.RooFit.What("N"), cut)
                if i+1==len(cats) and not j : w.pdf('model').paramOn(plt,
                                                                     r.RooFit.Layout(0.6,0.99,0.99),
                                                                     r.RooFit.Parameters(args),
                                                                     )
                if (cat,var) not in self.maxi : self.maxi[(cat,var)] = plt.GetMaximum()
                plt.SetTitle("" if not j else ("Channel: "+cat))
                plt.SetMaximum(self.maxi[(cat,var)]*1.1)
                plt.Draw()
                #if not (i or j) : leg.AddEntry('data'+cat+var,'Data', "LP")
                #for iComps in reversed(range(len(comps))) :
                #    if not (i or j) : leg.AddEntry(comps[iComps]+cat+var,comps[iComps],"f")
        print '!'
        c.cd(1)
        #leg.Draw()
        c.Print('%s.pdf'%fileName)

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

    @roo.quiet
    def contourProfileNLL(self, w, result, levels = [1.15, 3.0], nSpokes = 40, scale_dqq = 0.05, scale_Rag = 0.1, tolerance = 0.015) :
        nll = w.pdf('model').createNLL(*self.fitArgs)
        prof = r.RooProfileLL('prof','',nll,w.argSet('d_qq,R_ag'))

        def final() :
            fit = result.floatParsFinal()
            return fit.at(fit.index('d_qq')).getVal(), fit.at(fit.index('R_ag')).getVal()
        d_qq,R_ag = final()
        w.arg('d_qq').Print()
        w.arg('R_ag').Print()
        prof.Print()
        
        cache = {(None,None):None}
        previous = dict([(lev,1) for lev in levels])
        def getPoint( phi, level, last = (0,0), guess = 1) :
            if phi > math.pi : return (None,None)
            coord = ( R_ag + scale_Rag*guess*math.cos(phi),
                      d_qq + scale_dqq*guess*math.sin(phi) )
            [w.var(item).setVal(coord[i]) for i,item in enumerate(['R_ag','d_qq'])]
            if coord not in cache : cache[coord] = prof.getVal()
            line_guess = guess + (level - cache[coord]) * (guess - last[0]) / (cache[coord] - last[1])
            simple_parab_guess = guess*math.sqrt(level/cache[coord])

            if last[0] :
                b = (last[1]/last[0]**2 - cache[coord]/guess**2) / (last[0] - guess)
                a = last[1]/last[0]**2 - b * last[0]
                def pg(x) : 
                    cub = a*x**2 + b*x**3
                    return x if 2*abs(cub-level)<tolerance else pg(x*math.sqrt(level/cub))
                parab_guess = pg(guess)
            else :
                parab_guess = 0

            print "%.5f  %.2f"%(guess, cache[coord])
            previous[level] = guess
            return ( coord if abs(cache[coord]-level) < tolerance else
                     getPoint(phi, level, last = (guess,cache[coord]), guess = parab_guess if parab_guess else simple_parab_guess   ) )

        files = [open('data/nll_%.3f.txt'%lev,'w') for lev in levels]
        for iPhi in range(nSpokes) :
            print iPhi
            points = [getPoint( iPhi*2*math.pi/nSpokes, lev, guess = previous[lev]) for lev in levels]
            for p,f in zip(points,files) :
                print >>f, p[0],p[1],cache[p]
                f.flush()
        for f in files : f.close()

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