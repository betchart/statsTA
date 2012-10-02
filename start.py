import inputs
import ROOT as r

def roo_quiet(function) :
    def wrapped(*args,**kwargs) :
        r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.WARNING) #suppress info messages
        function(*args,**kwargs)
        r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.DEBUG) #re-enable all messages
    return wrapped
    
@roo_quiet
def wimport(w, *args, **kwargs) : getattr(w, "import")(*args,**kwargs)

def gaussian_delta( workspace, var, (mean,rel_unc), 
                    limits = (-0.5,0.5), symbol = "", units = "") :
    sym = symbol if symbol else var

    hat = r.RooConstVar('%s_hat'%var,'#hat{%s}%s'%(sym,units), mean)
    delta = r.RooRealVar('delta_%s'%var,'#delta_{%s}'%sym, 0, *limits )
    func = r.RooFormulaVar(var, sym+units, '(1+@0)*@1', r.RooArgList(delta, hat) )
    wimport(workspace,func)

    if not rel_unc : return
    err = r.RooConstVar('%s_rel_unc'%var,'#omega_{%s}/#hat{%s}'%(sym,sym), rel_unc )
    wimport(workspace, r.RooGaussian('%s_constraint'%var,'P(#delta_{%s})'%sym, 
                                     delta, r.RooFit.RooConst(0), err) )

def import_constraints(w) :
    gaussian_delta( w, "lumi", inputs.luminosity, symbol = "L", units = "(1/pb)" )
    [gaussian_delta( w, "xs_"+samp, xs, symbol = "#sigma_{%s}"%samp, units = "(pb)")
     for samp,xs in inputs.xs.items() if samp!='mj']
    gaussian_delta( w, 'xs_mj', inputs.xs['mj'], symbol = "#sigma_{%s}"%samp, units = "(pb)",
                    limits = (-0.99,100) )
    constr = r.RooProdPdf('constraints', 'l_{constraints}',
                          r.RooArgList(*[w.pdf('%s_constraint'%item)
                                         for item in (['lumi']+
                                                      ['xs_%s'%s for s in inputs.samples if s!='mj'])]))
    wimport(w, constr)

def import_efficiencies(w) :
    [wimport(w, r.RooConstVar(*(2*['_'.join(filter(None,['eff',chan,samp,comp]))]+[frac])))
     for chan,samps in inputs.efficiency.items()
     for samp,comps in samps.items()
     for comp,frac in comps]

def import_fractions(w) :
    compfracs = inputs.components['tt']
    assert zip(*compfracs)[0] == ('gg','qg','qq','ag')
    fhats = [r.RooConstVar('f_%s_hat'%comp,'#hat{f}_{%s}'%comp, frac) for comp,frac in compfracs]
    deltas = [r.RooRealVar('delta_%s'%comp,'#delta_{%s}'%comp, 0, -0.5, 0.5 ) for comp,frac in compfracs[:2] ]
    deltas.append( r.RooFormulaVar('delta_qq','#delta_{qq}', '(@3*(@5-@4) + (1+@5)*(@0*@4 + @1*@5))/(@3*(1+@4) + @2*(1+@5))', r.RooArgList(*(fhats+deltas)) ) )
    deltas.append( r.RooFormulaVar('delta_ag','#delta_{ag}', '(@0*@4 + @1*@5 + @2*@6)/@3', r.RooArgList(*(fhats+deltas)) ) )
    fs = [r.RooFormulaVar('f_%s'%comp, 'f_{%s}'%comp, '(1+@0)*@1', r.RooArgList(delta,fhat)) for (comp,frac),fhat,delta in zip(compfracs,fhats,deltas)]
    [wimport(w,item) for item in fs[:3]]
    wimport(w,fs[3], r.RooFit.RecycleConflictNodes())


w = r.RooWorkspace('Workspace')
import_constraints(w)
import_efficiencies(w)
import_fractions(w)
w.Print()
#plot = w.var('delta_lumi').frame()
#w.pdf('lumi_constraint').plotOn(plot)
#plot.Draw()
#raw_input()
