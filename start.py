import sys,math,inputs, ROOT as r
r.gROOT.SetBatch(1)

def roo_quiet(function) :
    def wrapped(*args,**kwargs) :
        msg = r.RooMsgService.instance() 
        cache = msg.globalKillBelow()            # get current message level
        #msg.setGlobalKillBelow(r.RooFit.WARNING) # suppress messages
        msg.setGlobalKillBelow(r.RooFit.ERROR) # suppress messages
        val = function(*args,**kwargs)                 # call function
        msg.setGlobalKillBelow(cache)            # resume prior message level
        return val
    return wrapped
    
@roo_quiet
def wimport(w, *args, **kwargs) : getattr(w, "import")(*args,**kwargs)

@roo_quiet
def wimport_const(w, name, value) : getattr(w, "import")(r.RooConstVar(*(2*[name]+[value])))

def factory(w, command) :  w.factory(command)

class topAsymmFit(object) :
    @roo_quiet
    def __init__(self) :
        self.observables = ['queuedbins','tridiscr']
        self.channels = dict((lepton,inputs.channel_data(lepton)) for lepton in ['el','mu'])
        print
        print self.channels['el']
        print
        print self.channels['mu']
        print
        self.ttcomps = ('qq','ag','gg','qg')

        w = self.setUp(useData = True)
        w.data('data').Print()
        w.data('data_el').Print()
        w.data('data_mu').Print()

        self.plot_fracs(w)
        self.defaults(w)
        self.print_fracs(w)
        self.print_n(w)
        #self.draw(w,'before')

        dqq = w.arg('d_qq')
        dqq.setConstant()
        output = open('dqq_scan.txt','w')
        for i in range(120) :
            dqq.setVal(0.2 - i*0.01)
            self.fitArgs = [w.data('data'),
                            r.RooFit.Extended(True),
                            r.RooFit.ExternalConstraints(r.RooArgSet(w.pdf('constraints'))),
                            #r.RooFit.Optimize(False),
                            r.RooFit.NumCPU(4)]
            result = w.pdf('model').fitTo(*(self.fitArgs+
                                            [r.RooFit.PrintLevel(-1),
                                             r.RooFit.Save(True),
                                             r.RooFit.PrintEvalErrors(-1),
                                             ]))
            print>>output, dqq.getVal(), w.arg('alphaL').getVal(), w.arg('alphaL').getError(), w.arg('f_gg').getVal(), w.arg('f_qg').getVal(), w.arg('f_qq').getVal(), w.arg('f_qg').getVal()
            #self.print_fracs(w)
            #self.print_n(w)
            #result.Print()
        output.close()
        #self.draw(w,'fit')
        #self.contourProfileNLL(w, result)

    def print_fracs(self,w) :
        for item in ['lumi_mu','lumi_el','f_gg','f_qg','f_qq','f_ag']+['xs_'+i for i in self.channels['el'].samples if i!='data'] : print "%s: %.04f"%(item, w.arg(item).getVal())
        print

    def print_n(self,w) :
        length = 24
        tots = {'el':0,'mu':0}
        print (' ').join(i.rjust(8) for i in ['']+tots.keys())
        for xs in self.channels['el'].samples :
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

    def setUp(self,useData) :
        w = r.RooWorkspace('Workspace')
        init_sequence = ['fractions','constraints','efficiencies','shapes','model','data']
        for item in init_sequence :
            print item,
            sys.stdout.flush()
            getattr(self, 'import_'+item)(*([w]+([useData] if item=='data' else [])))
            print '.'
        return w

    def import_fractions(self,w) :
        [wimport(w, r.RooConstVar("f_%s_hat"%comp,"#hat{f}_{%s}"%comp, self.channels['el'].samples['tt'+comp].frac)) for comp in self.ttcomps]
        factory(w, "R_ag[%f,0.07,1]"%(self.channels['el'].samples['ttag'].frac/self.channels['el'].samples['ttqq'].frac) )
        #factory(w, "R_ag[%f,%f,%f]"%(3*(self.channels['el'].samples['ttag'].frac/self.channels['el'].samples['ttqq'].frac,) ))
        factory(w, "d_qq[-0.9999999,1]")
        factory(w, "expr::f_qq('(1+@0)*@1',{d_qq,f_qq_hat})")
        factory(w, "prod::f_ag(R_ag,f_qq)")
        factory(w, "expr::f_qg('(1-@0-@1)/(1+@2*@3*@4/(@5*@6))',{f_qq,f_ag,R_ag,f_gg_hat,f_qq_hat,f_ag_hat,f_qg_hat})")
        factory(w, "expr::f_gg('1-@0-@1-@2',{f_qq,f_ag,f_qg})")

    def import_constraints(self,w) :
        factory(w, "Gaussian::lumi_constraint( d_lumi[0,-0.2,0.2], 0, %f)"%self.channels['el'].lumi_sigma)
        for channel in self.channels.values() :
            wimport_const( w, 'lumi_%s_hat'%channel.lepton, channel.lumi )
            factory(w, "expr::lumi_%s('(1+@0)*@1', {d_lumi, lumi_%s_hat})"%(channel.lepton,channel.lepton))

        xs_constraints = dict([(samp[:2],(data.xs,data.xs_sigma)) for samp,data in self.channels['el'].samples.items() if data.xs>0])

        for sample,(xs,delta) in xs_constraints.items() :
            wimport_const(w, 'xs_%s_hat'%sample, xs)
            factory(w, "Gaussian::xs_%s_constraint( d_xs_%s[0,-1,1.5], 0, %f)"%(sample,sample,delta))
            factory(w, "expr::xs_%s('(1+@0)*@1',{d_xs_%s, xs_%s_hat})"%(sample,sample,sample))
        
        factory(w, "PROD::constraints(%s)"%', '.join([s+'_constraint' for s in ['lumi']+['xs_'+t for t in xs_constraints]]))

        [factory(w, "prod::xs_tt%s(f_%s,xs_tt)"%(comp,comp)) for comp in self.ttcomps]
        wimport_const(w, 'xs_qcd', 100)
        
        
    def import_efficiencies(self,w) :
        [wimport(w, 
                 r.RooConstVar(*( 2*['eff_%s_%s'%(channel.lepton,sample)] + [data.eff] ) ) if data.eff>0 else
                 r.RooRealVar (*( 2*['eff_%s_%s'%(channel.lepton,sample)] + [0.01,0,1]) ))
         for channel in self.channels.values()
         for sample,data in channel.samples.items()
         if sample!='data']

    def import_shapes(self,w) :
        factory(w, "channel[%s]"%','.join("%s=%d"%(s,i) for i,s in enumerate(self.channels)))
        [factory(w, "%s[-1,1]"%obs) for obs in self.observables]

        [self.import_shape(w,channel.lepton,sample,data)
         for channel in self.channels.values()
         for sample,data in channel.samples.items()
         if sample!='data' ]
        
    @roo_quiet
    def import_shape(self,w,lepton,sample,data) :
        name = '_'.join([lepton,sample])
        arglist = r.RooArgList(*[w.var(o) for o in self.observables])
        argset = r.RooArgSet(arglist)

        for i,label in enumerate(['both','symm']) :
            wimport(w, r.RooDataHist('%s_sim_%s'%(name,label), '', arglist, data.datas[i]))
            wimport(w, r.RooHistPdf('%s_%s'%(name,label),'', argset, w.data('%s_sim_%s'%(name,label))))

        factory(w, "prod::expect_%s(%s)"%(name, ','.join(['lumi_'+lepton,
                                                          'xs_'+sample,
                                                          'eff_'+name])))

    def import_model(self,w) :
        wimport_const(w, 'alphaT', 1.0)
        factory(w, "alphaL[1, -15, 20]")

        [factory(w, "SUM::%(n)s( alphaT * %(n)s_both, %(n)s_symm )"%{'n':lepton+'_ttag'}) for lepton in self.channels]
        [factory(w, "SUM::%(n)s( alphaT * %(n)s_both, %(n)s_symm )"%{'n':lepton+'_ttqg'}) for lepton in self.channels]
        [factory(w, "SUM::%(n)s( alphaL * %(n)s_both, %(n)s_symm )"%{'n':lepton+'_ttqq'}) for lepton in self.channels]

        [factory(w, "SUM::model_%s( %s )"%(lepton, ','.join([ 'expect_%s_%s * %s_%s%s'%(lepton, key, lepton, key, value)
                                                              for key,value in {'qcd' :'_symm',
                                                                                'dy'  :'_symm',
                                                                                'wj'  :'_both',
                                                                                'st'  :'_both',
                                                                                'ttgg':'_both',
                                                                                'ttag':'',
                                                                                'ttqg':'',
                                                                                'ttqq':''}.items()
                                                              ]) ) ) for lepton in self.channels]

        factory(w, "SIMUL::model(channel, %s)"%', '.join("%(chan)s=model_%(chan)s"%{'chan':lepton} for lepton in self.channels))
        factory(w, "PROD::constrained_model( model, constraints )")


    @roo_quiet
    def import_data(self,w, useData) :
        obs_ = w.argSet(','.join(self.observables))
        obs = r.RooArgList(obs_)

        [w.var(v).setBins( getattr( self.channels['el'].samples['data'].datas[0], 'GetNbins'+X)() ) for v,X in zip(self.observables,'XY')]

        datas = [(lepton, r.RooDataHist('data_'+lepton,'N_{obs}^{%s}'%lepton,obs,chan.samples['data'].datas[0])) for lepton,chan in self.channels.items()]

        gens = [(lepton,w.pdf('model_'+lepton).generateBinned( obs_,
                                                               349030 if chan=='mu' else 338286, #r.RooFit.Extended(),
                                                               r.RooFit.Name('data_'+lepton))) for lepton in self.channels]

        [wimport(w, dat[1]) for dat in (datas if useData else gens)]
        args = [r.RooFit.Import(*dat) for dat in (datas if useData else gens)]
        wimport(w, r.RooDataHist('data', 'N_{obs}', 
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


    @roo_quiet
    def draw(self, w, fileName) :
        if not hasattr(self,'maxi') : self.maxi = {}
        leg = r.TLegend(0.7,0.5,0.99,0.99)
        c = r.TCanvas()
        cats = self.channels.keys()
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

    @roo_quiet
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

        files = [open('nll_%.3f.txt'%lev,'w') for lev in levels]
        for iPhi in range(nSpokes) :
            print iPhi
            points = [getPoint( iPhi*2*math.pi/nSpokes, lev, guess = previous[lev]) for lev in levels]
            for p,f in zip(points,files) :
                print >>f, p[0],p[1],cache[p]
                f.flush()
        for f in files : f.close()

if __name__=='__main__' : topAsymmFit()
