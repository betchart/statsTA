import sys
sys.argv.append('-b')
import inputs, ROOT as r

def roo_quiet(function) :
    def wrapped(*args,**kwargs) :
        msg = r.RooMsgService.instance() 
        cache = msg.globalKillBelow()            # get current message level
        msg.setGlobalKillBelow(r.RooFit.WARNING) # suppress messages
        function(*args,**kwargs)                 # call function
        msg.setGlobalKillBelow(cache)            # resume prior message level
    return wrapped
    
@roo_quiet
def wimport(w, *args, **kwargs) : getattr(w, "import")(*args,**kwargs)

def factory(w, command) : w.factory(command)

class topAsymmFit(object) :
    def __init__(self) :
        w = self.setUp(useData = True)
        w.data('data').Print()
        w.data('data_el').Print()
        w.data('data_mu').Print()
        self.print_fracs(w)
        self.print_n(w)
        self.draw(w,'before')
        result = w.pdf('model').fitTo(w.data('data'),
                                      r.RooFit.Extended(True),
                                      r.RooFit.ExternalConstraints(r.RooArgSet(w.pdf('constraints'))),
                                      r.RooFit.PrintLevel(-1),
                                      r.RooFit.Save(True),
                                      r.RooFit.PrintEvalErrors(-1),
                                      r.RooFit.NumCPU(4)
                                      )
        self.print_fracs(w)
        self.print_n(w)
        result.Print() 
        self.draw(w,'fit')
        #self.test()
        return

        mc = self.setUpModel(w)
        plc = r.RooStats.ProfileLikelihoodCalculator(w.data('data'), mc)
        plc.SetConfidenceLevel(.90)
        
        interval = plc.GetInterval()
        limits = dict([(a,(interval.LowerLimit(w.arg(a)),interval.UpperLimit(w.arg(a)))) for a in ['d_qq','d_ag']])
        print 'cl',interval.ConfidenceLevel()
        print limits

    def print_fracs(self,w) :
        for item in ['lumi_mu','lumi_el','f_gg','f_qg','f_qq','f_ag']+['xs_'+i for i in inputs.xs] : print "%s: %.04f"%(item, w.arg(item).getVal())
        print

    def print_n(self,w) :
        length = 24
        tots = {'el':0,'mu':0}
        print (' ').join(i.rjust(8) for i in ['']+tots.keys())
        for xs in inputs.xs :
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


    def setUp(self,useData) :
        w = r.RooWorkspace('Workspace')
        init_sequence = ['fractions','constraints','efficiencies','shapes','model','data']
        for item in init_sequence : 
            print item,
            sys.stdout.flush()
            getattr(self, 'import_'+item)(*([w]+([useData] if item=='data' else [])))
            print '.'
        return w

    def setUpModel(self, w) :
        mc = r.RooStats.ModelConfig(w)
        mc.SetObservables(r.RooArgSet(w.var('d3'),w.var('ptpt'),w.arg('channel')))
        poi = r.RooArgSet(w.arg('d_qq'),w.arg('d_ag'))
        mc.SetParametersOfInterest(poi)
        nuis = r.RooArgSet(*[w.arg(i) for i in ['d_lumi_mu','eff_el_mj','eff_mu_mj','global_R_el','global_R_mu']+['d_xs_'+j for j in inputs.xs if j!='mj']])
        mc.SetNuisanceParameters(nuis)
        ##mc.SetConstraintParameters(r.RooArgSet(*[w.arg(i) for i in 'd_lumi,d_xs_st,d_xs_dy'.split(',')]))
        mc.SetPdf(w.pdf('constrained_model'))
        return mc

    def import_fractions(self,w) :
        comps,fracs = zip(*inputs.components['tt'])
        assert comps == ('qq','ag','gg','qg')
        [wimport(w, r.RooConstVar("f_%s_hat"%comp,"#hat{f}_{%s}"%comp, frac)) for comp,frac in zip(comps,fracs)]
        [factory(w, "d_%s[0,-1,3]"%comp) for comp in comps[:2]]
        args = sum(zip(*[['f_%s_hat'%comp,'d_'+comp] for comp in comps]),())
        factory(w, "expr::d_gg('((1+@5)*(@3-@4*@0-@5*@1) - (1+@4)*@3) / ((1+@4)*@3 + (1+@5)*@2)',{%s})"%(', '.join(args[:-2])))
        factory(w, "expr::d_qg('-(@0*@4 + @1*@5 + @2*@6)/@3', {%s})"%(', '.join(args[:-1])))
        [factory(w, "expr::f_%s('(1+@0)*@1',{d_%s,f_%s_hat})"%tuple([comp]*3)) for comp in comps]

    def import_constraints(self,w) :
        def gdelta( w, var, (mean,rel_unc), limits = (-1.5,1.5), delta = None) :
            if not delta : delta = 'd_%s[0%s]'%(var, (",%f,%f"%limits) if rel_unc else '')
            factory(w, "expr::%s('(1+@0)*@1',{ %s, %s_hat[%f]})"%(var,delta,var,mean))
            if rel_unc : factory(w, "Gaussian::%s_constraint( d_%s, 0, %f)"%(var,var,rel_unc))
            return
        gdelta( w, "lumi_mu", inputs.luminosity['mu'], limits = (-0.2,0.2))
        gdelta( w, "lumi_el", inputs.luminosity['el'], limits = (-0.2,0.2), delta='d_lumi_mu')
        [gdelta( w, "xs_"+samp, xs) for samp,xs in inputs.xs.items() if samp!='mj']
        gdelta( w, 'xs_mj', inputs.xs['mj'], limits = (-1,1) )
        factory(w, "PROD::constraints(%s)"%', '.join([s+'_constraint' for s in ['lumi_mu']+['xs_'+t for t in inputs.components if t!='mj']]))
        
    def import_shapes(self,w) :
        factory(w, "channel[%s]"%','.join("%s=%d"%(s,i) for i,s in enumerate(inputs.efficiency)))
        factory(w, "d3[-1,1]")
        factory(w, "ptpt[0,1]")
        [self.import_shape(w,chan,samp,comp)
         for chan,samps in inputs.efficiency.items()
         for samp,comps in samps.items()
         for comp,frac in comps]
        [factory(w, "expr::expect_%s_tt('%s',{%s})"%(chan,
                                                    '+'.join('@%d'%i for i in range(len(inputs.components['tt']))),
                                                    ','.join('expect_%s_tt%s'%(chan,c) for c,_ in inputs.components['tt'])))
         for chan in inputs.efficiency]

    @staticmethod
    @roo_quiet
    def import_shape(w,chan,samp,comp) :
        lab = ['d3','ptpt']
        name = '_'.join([chan,samp+comp])
        [wimport(w, r.RooDataHist('%s_%s_sim'%(i,name), '', r.RooArgList(w.var(i)), inputs.histogram(i,chan,samp,comp))) for i in lab]
        [factory(w, "HistPdf::%s(%s,%s)"%('%s_%s'%(i,name),i, '%s_%s_sim'%(i,name))) for i in lab]
        factory(w, "PROD::%s(%s)"%( '%s_%s'%('-'.join(lab),name), ','.join('%s_%s'%(i,name) for i in lab)))
        factory(w, "expr::%s('%s',{%s})"%('expect_'+name,
                                         '*'.join('@%d'%i for i in range(5 if comp else 4)),
                                         ','.join(['lumi_'+chan,'global_R_'+chan,'xs_'+samp,'eff_'+name,'f_'+comp][:None if comp else -1])))

    def import_efficiencies(self,w) :
        [wimport(w,
                 r.RooConstVar(*(2*['_'.join(filter(None,['eff',chan,samp+comp]))]+[eff])) if eff else
                 r.RooRealVar(*(2*['_'.join(filter(None,['eff',chan,samp+comp]))] + [0.01,0,1])))
         for chan,samps in inputs.efficiency.items()
         for samp,comps in samps.items() 
         for comp,eff in comps]
        [wimport(w,r.RooRealVar('global_R_%s'%chan, 'R_{%s}'%chan, 1, 0.5,2)) for chan in inputs.efficiency]

    def import_model(self,w) :
        [factory(w, "SUM::model_%s( %s )"%(chan,','.join('expect_%s_%s * d3-ptpt_%s_%s'%(2*(chan,samp+comp)) 
                                                        for samp,comps in samp.items() 
                                                        for comp,_ in comps))) for chan,samp in inputs.efficiency.items()]
        factory(w, "SIMUL::model(channel, %s)"%', '.join("%s=%s"%(chan,'model_'+chan) for chan in inputs.efficiency))
        factory(w, "PROD::constrained_model( model, constraints )")

    @roo_quiet
    def import_data(self,w, useData) :
        obs_ = w.argSet('d3,ptpt')
        obs = r.RooArgList(obs_)
        hists = [(chan, inputs.histogram(None,chan,chan if chan=='mu' else 'EleHad.2011',d=2)) for chan in inputs.efficiency]
        [w.var(v).setBins(hists[0][1].GetNbinsX() if v=='d3' else hists[0][1].GetNbinsY()) for v in ['d3','ptpt']]
        datas = [(chan, r.RooDataHist('data_'+chan,'N_{obs}^{%s}'%chan,obs,hist)) for chan,hist in hists]
        gens = [(chan,w.pdf('model_'+chan).generateBinned(obs_,
                                                          39055 if chan=='mu' else 19528, #r.RooFit.Extended(),
                                                          r.RooFit.Name('data_'+chan))) for chan in inputs.efficiency]
        [wimport(w, dat[1]) for dat in (datas if useData else gens)]
        args = [r.RooFit.Import(*dat) for dat in (datas if useData else gens)]
        wimport(w, r.RooDataHist('data', 'N_{obs}', 
                                 obs,
                                 r.RooFit.Index(w.arg('channel')),
                                 *args
                                 ),
                )

    def test(self, w) :
        wimport(w, r.RooFormulaVar('sum','','@0*@1+@2*@3+@4*@5+@6*@7',r.RooArgList(*sum([[w.arg('d_'+i),w.arg('f_'+i+'_hat')] for i in ['gg','qg','qq','ag']],[]))))
        wimport(w, r.RooFormulaVar('prod','','(1+@0)*(1+@2)-(1+@1)*(1+@3)',r.RooArgList(*[w.arg('d_'+i) for i in ['gg','qg','qq','ag']])))


    @roo_quiet
    def draw(self, w, fileName) :
        if not hasattr(self,'maxi') : self.maxi = {}
        c = r.TCanvas()
        cats = inputs.efficiency.keys()
        data = w.data('data')
        c.Divide(2,len(cats))
        for i,cat in enumerate(cats) :
            for j,var in enumerate(['d3','ptpt']) :
                c.cd(i+1+2*j)
                plt = w.var(var).frame()
                cut = r.RooFit.Cut("channel==channel::%s"%cat)
                data.plotOn(plt, cut, r.RooFit.MarkerSize(0.5))
                if (cat,var) not in self.maxi : self.maxi[(cat,var)] = plt.GetMaximum()
                if j : data.statOn(plt, r.RooFit.What("N"), cut)
                comps = ['st','dy','mj','wj','ttgg','ttqg','ttag','ttqq']
                colrs = [r.kGray, r.kGray, r.kGreen, r.kRed, r.kBlue+3, r.kBlue+2, r.kBlue+1, r.kViolet]
                w.pdf('model').plotOn( plt,
                                       r.RooFit.ProjWData(w.argSet(''),data),
                                       r.RooFit.Slice(w.cat('channel'),cat),
                                       r.RooFit.LineWidth(1),
                                       r.RooFit.LineColor(r.kOrange)
                                       )
                for iComps in range(len(comps)) :
                    col = colrs[iComps]
                    compstr = ','.join(['d3-ptpt_%s_%s'%(cat,sub) for sub in comps[:iComps+1]])
                    w.pdf('model').plotOn( plt,
                                           r.RooFit.ProjWData(w.argSet(''),data),
                                           r.RooFit.Slice(w.cat('channel'),cat),
                                           r.RooFit.Components(compstr),
                                           r.RooFit.LineWidth(0),
                                           r.RooFit.LineColor(col),
                                           r.RooFit.FillColor(col),
                                           r.RooFit.DrawOption('F'),
                                           r.RooFit.MoveToBack()
                                           )
                args = w.argSet(','.join(['d_'+t for t in ['gg','qg','qq','ag']] + ['global_R_mu','global_R_el']))
                if i+1==len(cats) and not j : w.pdf('model').paramOn(plt,
                                                                     r.RooFit.Layout(0.6,0.99,0.99),
                                                                     r.RooFit.Parameters(args),
                                                                     )
                plt.SetTitle("" if not j else ("Channel: "+cat))
                plt.SetMaximum(self.maxi[(cat,var)] * (1.1 if var=='d3' else 1.3))
                plt.Draw()
        c.Print('%s.pdf'%fileName)

if __name__=='__main__' : topAsymmFit()
