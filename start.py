import inputs,sys
import ROOT as r

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

class topAsymmFit(object) :
    def __init__(self) :
        self.useData = True
        self.keep = []
        w = r.RooWorkspace('Workspace')
        init_sequence = ['fractions','constraints','efficiencies','shapes','model','data']
        for item in init_sequence : 
            print item,
            sys.stdout.flush()
            getattr(self, 'import_'+item)(w)
            print '.'


        def print_fracs() :
            for item in ['lumi_mu','lumi_el','f_gg','f_qg','f_qq','f_ag']+['xs_'+i for i in inputs.xs] : print "%s: %.04f"%(item, w.arg(item).getVal())
        def print_n() :
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

        print_fracs() ; print
        print_n() ; print
        result = w.pdf('model').fitTo(w.data('data'),
                                      r.RooFit.Extended(True),
                                      r.RooFit.ExternalConstraints(r.RooArgSet(w.pdf('constraints'))),
                                      r.RooFit.PrintLevel(-1),
                                      r.RooFit.Save(True),
                                      r.RooFit.PrintEvalErrors(-1))
        print
        print_fracs();print
        print_n();print
        w.data('data').Print()
        w.data('data_el').Print()
        w.data('data_mu').Print()
        result.Print()
        #wimport(w, r.RooFormulaVar('sum','','@0*@1+@2*@3+@4*@5+@6*@7',r.RooArgList(*sum([[w.arg('d_'+i),w.arg('f_'+i+'_hat')] for i in ['gg','qg','qq','ag']],[]))))
        #wimport(w, r.RooFormulaVar('prod','','(1+@0)*(1+@2)-(1+@1)*(1+@3)',r.RooArgList(*[w.arg('d_'+i) for i in ['gg','qg','qq','ag']])))

        mc = r.RooStats.ModelConfig(w)
        mc.SetObservables(r.RooArgSet(w.var('d3'),w.var('ptpt'),w.arg('channel')))
        poi = r.RooArgSet(w.arg('d_qq'),w.arg('d_ag'))
        mc.SetParametersOfInterest(poi)
        nuis = r.RooArgSet(*[w.arg(i) for i in ['d_lumi_mu','eff_el_mj','eff_mu_mj','global_R_el','global_R_mu']+['d_xs_'+j for j in inputs.xs if j!='mj']])
        mc.SetNuisanceParameters(nuis)
        ##mc.SetConstraintParameters(r.RooArgSet(*[w.arg(i) for i in 'd_lumi,d_xs_st,d_xs_dy'.split(',')]))
        mc.SetPdf(w.pdf('constrained_model'))

        plc = r.RooStats.ProfileLikelihoodCalculator(w.data('data'), mc)
        plc.SetTestSize(.05);
        interval = plc.GetInterval()
        dqq_LL = interval.LowerLimit(w.arg('d_qq'))
        dqq_UL = interval.UpperLimit(w.arg('d_qq'))
        dag_LL = interval.LowerLimit(w.arg('d_ag'))
        dag_UL = interval.UpperLimit(w.arg('d_ag'))
        
        print 'cl',interval.ConfidenceLevel()
        print 'll-dqq',dqq_LL
        print 'ul-dqq',dqq_UL
        print 'll-dag',dag_LL
        print 'ul-dag',dag_UL
        #w.Print()

    def import_fractions(self,w) :
        compfracs = inputs.components['tt']
        assert zip(*compfracs)[0] == ('qq','ag','gg','qg')
        fhats = [r.RooConstVar('f_%s_hat'%comp,'#hat{f}_{%s}'%comp, frac) for comp,frac in compfracs]
        deltas = [r.RooRealVar('d_%s'%comp,'#delta_{%s}'%comp, 0, -1, 3 ) for comp,frac in compfracs[:2] ]
        deltas.append( r.RooFormulaVar('d_gg','#delta_{gg}', '((1+@5)*(@3-@4*@0-@5*@1) - (1+@4)*@3) / ((1+@4)*@3 + (1+@5)*@2)', r.RooArgList(*(fhats+deltas)) ) )
        deltas.append( r.RooFormulaVar('d_qg','#delta_{qg}', '-(@0*@4 + @1*@5 + @2*@6)/@3', r.RooArgList(*(fhats+deltas)) ) )
        fs = [r.RooFormulaVar('f_%s'%comp, 'f_{%s}'%comp, '(1+@0)*@1', r.RooArgList(delta,fhat)) for (comp,frac),fhat,delta in zip(compfracs,fhats,deltas)]
        [wimport(w,item) for item in fs[:3]]
        wimport(w,fs[3], r.RooFit.RecycleConflictNodes())

    def import_constraints(self,w) :
        self.gaussian_delta( w, "lumi_mu", inputs.luminosity['mu'], symbol = "Lmu", units = "(1/pb)" , limits = (-0.2,0.2))
        self.gaussian_delta( w, "lumi_el", inputs.luminosity['el'], symbol = "Lel", units = "(1/pb)" , limits = (-0.2,0.2), delta=w.var('d_lumi_mu'))
        [self.gaussian_delta( w, "xs_"+samp, xs, symbol = "#sigma_{%s}"%samp, units = "(pb)") for samp,xs in inputs.xs.items() if samp!='mj']
        self.gaussian_delta( w, 'xs_mj', inputs.xs['mj'], symbol = "#sigma_{%s}"%samp, units = "(pb)", limits = (-1,1) )
        
        #F = r.RooConstVar('F','F',10000)
        #wimport(w, r.RooGenericPdf("ordering_constraint","P(f_{gg}-f_{qg})", "(exp(@2*(@0-@1))/(1+exp(@2*(@0-@1))))*(exp(@2*(1-@1-@0))/(1+exp(@2*(1-@1-@0))))", r.RooArgList(w.arg('f_gg'),w.arg('f_qg'),F)))
        #F = r.RooFormulaVar('F','F',"@0/@1",r.RooArgList(w.arg('f_gg'),w.arg('f_qg')))
        #Fmc = r.RooConstVar('Fmc','F_{mc}',2.75)
        #sFmc = r.RooConstVar('sFmc','#sigma_{Fmc}',0.2)
        #wimport(w, r.RooGaussian("ordering_constraint","P(f_{gg}-f_{qg})", F,Fmc,sFmc))
        constr = r.RooProdPdf('constraints', 'l_{constraints}',
                              r.RooArgList(*[w.pdf('%s_constraint'%item)
                                             for item in (['lumi_mu','ordering'][:1]+
                                                          ['xs_%s'%s for s in inputs.components if s!='mj'])]))
        wimport(w, constr)

    @staticmethod
    def gaussian_delta( workspace, var, (mean,rel_unc), limits = (-1.5,1.5), symbol = "", units = "", delta = None) :
        sym = symbol if symbol else var
        hat = r.RooConstVar('%s_hat'%var,'#hat{%s}%s'%(sym,units), mean)
        if not delta :
            delta = ( r.RooRealVar('d_%s'%var,'#delta_{%s}'%sym, 0, *limits ) if rel_unc else
                      r.RooConstVar('d_%s'%var,'#delta_{%s}'%sym, 0))
        func = r.RooFormulaVar(var, sym+units, '(1+@0)*@1', r.RooArgList(delta, hat) )
        wimport(workspace,func)
        if not rel_unc : return
        err = r.RooConstVar('%s_rel_unc'%var,'#omega_{%s}/#hat{%s}'%(sym,sym), rel_unc )
        wimport(workspace, r.RooGaussian('%s_constraint'%var,'P(#delta_{%s})'%sym, delta, r.RooFit.RooConst(0), err) )


    def import_shapes(self,w) :
        ch = r.RooCategory('channel','chan')
        for chan in inputs.efficiency : ch.defineType(chan)
        [wimport(w,item) for item in [ch,
                                      r.RooRealVar('d3','D_{3}',-1,1),
                                      r.RooRealVar('ptpt','t#bar{t}.pt/(t.pt+#bar{t}.pt)',0,1) ]]
        [self.import_shape(w,chan,samp,comp)
         for chan,samps in inputs.efficiency.items()
         for samp,comps in samps.items()
         for comp,frac in comps]
        for chan in inputs.efficiency : wimport(w, r.RooFormulaVar('expect_%s_tt'%chan, '#lambda_{tt}^{%s}'%chan,
                                                                   '+'.join('@%d'%i for i in range(len(inputs.components['tt']))),
                                                                   r.RooArgList(*[w.arg('expect_%s_tt%s'%(chan,c)) for c,_ in inputs.components['tt']])) )

    @staticmethod
    @roo_quiet
    def import_shape(w,chan,samp,comp) :
        lab = ['d3','ptpt']
        var  = [w.var(i) for i in lab]
        hist = [inputs.histogram(i,chan,samp,comp) for i in lab]
        rdh = [r.RooDataHist('%s_%s_%s_sim'%(i,chan,samp+comp), '', r.RooArgList(a), h) for i,a,h in zip(lab,var,hist)]
        pdf = [r.RooHistPdf('%s_%s_%s'%(i,chan,samp+comp),'P(%s)_%s'%(lab,samp+comp),r.RooArgSet(a), dh) for i,a,dh in zip(lab,var,rdh)]
        wimport(w, r.RooProdPdf('%s_%s_%s'%('-'.join(lab),chan,samp+comp),'P(%s,%s)_%s^{%s}'%tuple(lab+[samp+comp,chan]),*pdf))
        wimport(w, r.RooFormulaVar( 'expect_%s_%s'%(chan,samp+comp), '#lambda_{%s}^{%s}'%(samp+comp,chan), '*'.join(['@%d'%i for i in range(5 if comp else 4)]), 
                                    r.RooArgList(*[w.arg(i) for i in ['lumi_'+chan,'global_R_'+chan,'xs_%s'%samp,'_'.join(['eff',chan,samp+comp])] + ([] if not comp else ['f_%s'%comp])]) ) )
        return hist

    def import_efficiencies(self,w) :
        [wimport(w,
                 r.RooConstVar(*(2*['_'.join(filter(None,['eff',chan,samp+comp]))]+[eff])) if eff else
                 r.RooRealVar(*(2*['_'.join(filter(None,['eff',chan,samp+comp]))] + [0.01,0,1])))
         for chan,samps in inputs.efficiency.items()
         for samp,comps in samps.items() 
         for comp,eff in comps]
        [wimport(w,r.RooRealVar('global_R_%s'%chan, 'R_{%s}'%chan, 1, 0.5,2)) for chan in inputs.efficiency]

    def import_model(self,w) :
        chmodels = [ r.RooAddPdf('model_%s'%chan,'model_{%s}'%chan, 
                                 r.RooArgList(*[w.pdf('d3-ptpt_%s_%s'%(chan,samp+comp)) for samp,comps in samps.items() for comp,_ in comps]), 
                                 r.RooArgList(*[w.arg('expect_%s_%s'%(chan,samp+comp)) for samp,comps in samps.items() for comp,_ in comps]) )
                     for chan,samps in inputs.efficiency.items() ]
        model = r.RooSimultaneous('model','d3-ptpt-ch', r.RooArgList(*chmodels), w.arg('channel')) 
        wimport(w, r.RooProdPdf('constrained_model','d3-ptpt-ch constrained', r.RooArgList(model,w.pdf('constraints')) ))

    @roo_quiet
    def import_data(self,w) :
        obs_ = r.RooArgSet(w.var('d3'), w.var('ptpt'))
        obs = r.RooArgList(w.var('d3'), w.var('ptpt'))

        gens = [(chan,w.pdf('model_'+chan).generateBinned(obs_, 39055 if chan=='mu' else 19528,#r.RooFit.Extended(),
                                                            r.RooFit.Name('genBinned_'+chan))) for chan in inputs.efficiency]
        hists = [(chan, inputs.histogram(None,chan,chan if chan=='mu' else 'EleHad.2011',d=2)) for chan in inputs.efficiency]
        datas = [(chan, r.RooDataHist('data_'+chan,'N_{obs}^{%s}'%chan,obs,hist)) for chan,hist in hists]
        self.keep.append(hists)
        [wimport(w, dat[1]) for dat in (datas if self.useData else gens)]
        args = [r.RooFit.Import(*dat) for dat in (datas if self.useData else gens)]
        wimport(w, r.RooDataHist('data', 'N_{obs}', 
                                 obs,
                                 r.RooFit.Index(w.arg('channel')),
                                 *args
                                 ),
                )

if __name__=='__main__' : topAsymmFit()
