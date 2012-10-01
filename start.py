import ROOT as r

def roo_quiet(function) :
    def wrapped(*args,**kwargs) :
        r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.WARNING) #suppress info messages
        function(*args,**kwargs)
        r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.DEBUG) #re-enable all messages
    return wrapped
    
@roo_quiet
def wimport(w, item) : getattr(w, "import")(item)

def gaussian_delta( workspace, var,  mean, rel_unc = None, limits = (None,None) , symbol = "", units = "" ) :
    sym = symbol if symbol else var

    hat = r.RooConstVar('%s_hat'%var,'#hat{%s}%s'%(sym,units), mean)
    delta = r.RooRealVar('delta_%s'%var,'#delta_{%s}'%sym, 0, *limits )
    func = r.RooFormulaVar(var, sym+units, '(1+@0)*@1', r.RooArgList(delta, hat) )
    wimport(workspace,func)

    if not rel_unc : return
    err = r.RooConstVar('%s_rel_unc'%var,'#omega_{%s}/#hat{%s}'%(sym,sym), rel_unc )
    constraint = r.RooGaussian('%s_constraint'%var,'P(#delta_{%s})'%sym, delta, r.RooFit.RooConst(0), err)
    wimport(workspace,constraint)

w = r.RooWorkspace('Workspace')
gaussian_delta( w, "lumi", 5008.000, 0.04, (-0.5,0.5), symbol = "L", units = "(1/pb)" )
gaussian_delta( w, "xs_tt", 149.600, 0.30, (-0.5,0.5), symbol = "#sigma_{tt}", units = "(pb)")
gaussian_delta( w, "xs_wj",1911.800, 0.30, (-0.5,0.5), symbol = "#sigma_{wj}", units = "(pb)")
gaussian_delta( w, "xs_mj",   1.000, None,(-0.99,1e6), symbol = "#sigma_{mj}", units = "(pb)")
gaussian_delta( w, "xs_st",  71.968, 0.01, (-0.5,0.5), symbol = "#sigma_{st}", units = "(pb)")
gaussian_delta( w, "xs_dy",2475.000, 0.01, (-0.5,0.5), symbol = "#sigma_{dy}", units = "(pb)")



w.Print()
#plot = w.var('delta_lumi').frame()
#w.pdf('lumi_constraint').plotOn(plot)
#plot.Draw()
#raw_input()
