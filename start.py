import inputs
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
        w = r.RooWorkspace('Workspace')
        init_sequence = ['constraints','fractions','efficiencies','shapes','model','data']
        for item in init_sequence : getattr(self, 'import_'+item)(w)

        wimport(w, r.RooFormulaVar('sum','','@0*@1+@2*@3+@4*@5+@6*@7',r.RooArgList(*sum([[w.arg('d_'+i),w.arg('f_'+i+'_hat')] for i in ['gg','qg','qq','ag']],[]))))
        wimport(w, r.RooFormulaVar('prod','','(1+@0)*(1+@2)-(1+@1)*(1+@3)',r.RooArgList(*[w.arg('d_'+i) for i in ['gg','qg','qq','ag']])))
        w.Print()

    def import_constraints(self,w) :
        self.gaussian_delta( w, "lumi", inputs.luminosity, symbol = "L", units = "(1/pb)" )
        [self.gaussian_delta( w, "xs_"+samp, xs, symbol = "#sigma_{%s}"%samp, units = "(pb)")
         for samp,xs in inputs.xs.items() if samp!='mj']
        self.gaussian_delta( w, 'xs_mj', inputs.xs['mj'], symbol = "#sigma_{%s}"%samp, units = "(pb)",
                             limits = (-0.99,100) )
        constr = r.RooProdPdf('constraints', 'l_{constraints}',
                              r.RooArgList(*[w.pdf('%s_constraint'%item)
                                             for item in (['lumi']+
                                                          ['xs_%s'%s for s in inputs.samples if s!='mj'])]))
        wimport(w, constr)

    @staticmethod
    def gaussian_delta( workspace, var, (mean,rel_unc), limits = (-0.5,0.5), symbol = "", units = "") :
        sym = symbol if symbol else var
        hat = r.RooConstVar('%s_hat'%var,'#hat{%s}%s'%(sym,units), mean)
        delta = r.RooRealVar('d_%s'%var,'#delta_{%s}'%sym, 0, *limits )
        func = r.RooFormulaVar(var, sym+units, '(1+@0)*@1', r.RooArgList(delta, hat) )
        wimport(workspace,func)
        if not rel_unc : return
        err = r.RooConstVar('%s_rel_unc'%var,'#omega_{%s}/#hat{%s}'%(sym,sym), rel_unc )
        wimport(workspace, r.RooGaussian('%s_constraint'%var,'P(#delta_{%s})'%sym, delta, r.RooFit.RooConst(0), err) )


    def import_fractions(self,w) :
        compfracs = inputs.components['tt']
        assert zip(*compfracs)[0] == ('gg','qg','qq','ag')
        fhats = [r.RooConstVar('f_%s_hat'%comp,'#hat{f}_{%s}'%comp, frac) for comp,frac in compfracs]
        deltas = [r.RooRealVar('d_%s'%comp,'#delta_{%s}'%comp, 0, -0.5, 0.5 ) for comp,frac in compfracs[:2] ]
        deltas.append( r.RooFormulaVar('d_qq','#delta_{qq}', '((1+@5)*(@3-@4*@0-@5*@1) - (1+@4)*@3) / ((1+@4)*@3 + (1+@5)*@2)', r.RooArgList(*(fhats+deltas)) ) )
        deltas.append( r.RooFormulaVar('d_ag','#delta_{ag}', '-(@0*@4 + @1*@5 + @2*@6)/@3', r.RooArgList(*(fhats+deltas)) ) )
        fs = [r.RooFormulaVar('f_%s'%comp, 'f_{%s}'%comp, '(1+@0)*@1', r.RooArgList(delta,fhat)) for (comp,frac),fhat,delta in zip(compfracs,fhats,deltas)]
        [wimport(w,item) for item in fs[:3]]
        wimport(w,fs[3], r.RooFit.RecycleConflictNodes())


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

    @staticmethod
    @roo_quiet
    def import_shape(w,chan,samp,comp) :
        lab = ['d3','ptpt']
        var  = [w.var(i) for i in lab]
        hist = [inputs.histogram(i,chan,samp,comp) for i in lab]
        rdh = [r.RooDataHist('%s_%s_%s_sim'%(i,chan,samp+comp), '', r.RooArgList(a), h) for i,a,h in zip(lab,var,hist)]
        pdf = [r.RooHistPdf('%s_%s_%s'%(i,chan,samp+comp),'P(%s)_%s'%(lab,samp+comp),r.RooArgSet(a), dh) for i,a,dh in zip(lab,var,rdh)]
        wimport(w, r.RooProdPdf('%s_%s_%s'%('-'.join(lab),chan,samp+comp),'P(%s,%s)_%s^{%s}'%tuple(lab+[samp+comp,chan]),*pdf))
        wimport(w, r.RooFormulaVar( 'expect_%s_%s'%(chan,samp+comp), '#lambda_{%s}^{%s}'%(samp+comp,chan), '*'.join(['@%d'%i for i in range(4 if comp else 3)]), 
                                    r.RooArgList(*[w.arg(i) for i in ['lumi','xs_%s'%samp,'_'.join(['eff',chan,samp+comp])] + ([] if not comp else ['f_%s'%comp])]) ) )

    def import_efficiencies(self,w) :
        [wimport(w, r.RooConstVar(*(2*['_'.join(filter(None,['eff',chan,samp+comp]))]+[frac])))
         for chan,samps in inputs.efficiency.items()
         for samp,comps in samps.items() 
         for comp,frac in comps]

    def import_model(self,w) :
        chmodels = [ r.RooAddPdf('model_%s'%chan,'model_{%s}'%chan, 
                                 r.RooArgList(*[w.pdf('d3-ptpt_%s_%s'%(chan,samp+comp)) for samp,comps in samps.items() for comp,_ in comps]), 
                                 r.RooArgList(*[w.arg('expect_%s_%s'%(chan,samp+comp)) for samp,comps in samps.items() for comp,_ in comps]) )
                     for chan,samps in inputs.efficiency.items() ]
        wimport(w, r.RooSimultaneous('model','d3-ptp-ch', r.RooArgList(*chmodels), w.arg('channel')) )

    @roo_quiet
    def import_data(self,w) :
        datas = [(chan,inputs.histogram(None,None,None,None,2)) for chan in inputs.efficiency]
        args = [r.RooFit.Import(*dat) for dat in datas]
        wimport(w, r.RooDataHist('data', 'N_{obs}', 
                                 r.RooArgList(w.var('d3'), w.var('ptpt')),
                                 r.RooFit.Index(w.arg('channel')),
                                 *args
                                 ),
                )

if __name__=='__main__' : topAsymmFit()
