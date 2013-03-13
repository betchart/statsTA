import roo,ROOT as r
from systematics import systematics
from parabola import parabola

class dqqSlice(object): 
    
    @staticmethod
    def columns() : return ['dqq','minNLL','alphaL','err',
                            'profile_up','profile_down',
                            'MINUIT_cov_quality','alphaT','R_ag',
                            'f_qq','f_qg','f_gg','f_ag','d_xs_tt','d_xs_wj','eff_el_qcd','eff_mu_qcd',
                            'A','B','C'
                            ]

    def __str__(self) :
        return '\t'.join('%+.4f'%self.results[item] if item in self.results else ' 0.0000' for item in self.columns())

    @staticmethod
    def getFitValue(parName, fit, error=False) :
        args = fit.floatParsFinal()
        i = args.index(parName)
        return args.at(i).getVal() if not error else (args.at(i).getVal(), args.at(i).getError())


    def __init__(self, workspace, fitArgs, varToProfile = 'alphaL' ) :
        res = {'dqq': workspace.arg('d_qq').getVal()}

        model = workspace.pdf('model')
        data = workspace.data('data')

        central = model.fitTo( data, *(fitArgs+[r.RooFit.Save(True)]))
        nll = model.createNLL(data, *fitArgs[:-1] )
        pll = nll.createProfile(workspace.argSet(varToProfile))

        res['MINUIT_cov_quality'] = central.covQual()

        aL,aLe = self.getFitValue(varToProfile,central,error=True)
        aLarg = workspace.arg(varToProfile)

        def rooeval(yvar,x) :
            aLarg.setVal(x)
            return yvar.getVal()

        prbl = parabola([(x,rooeval(pll,x)) for x in [aL+aLe,aL-aLe,aL]])
        pll_1sigma = 0.5
        aLarg.setVal(prbl.xmin)

        res['minNLL'] = nll.getVal()
        res[varToProfile] = prbl.xmin
        res['err'] = prbl.dx(pll_1sigma)
        res['A'],res['B'],res['C'] = prbl.ABC
        pll.getVal()
        for item in ['alphaT','f_qq','f_qg','f_ag','f_gg'] : res[item] = workspace.arg(item).getVal()
        for item in ['d_xs_tt','d_xs_wj','eff_el_qcd','eff_mu_qcd'] : res[item] = workspace.arg(item).getVal()

        for sign,lab in zip([-1,1],['down','up']) :
            xTrial = res[varToProfile]+sign*res['err']
            dY = rooeval(pll,xTrial) - pll_1sigma
            dX = -2*dY / prbl.slope(xTrial)
            bounds = sorted([xTrial, xTrial+dX])
            res['profile_%s'%lab] = xTrial if abs(dY)<2e-5 else pll.findRoot(aLarg, bounds[0], bounds[1], pll_1sigma)

        pll.IsA().Destructor(pll)
        nll.IsA().Destructor(nll)
        del central

        self.results = res
