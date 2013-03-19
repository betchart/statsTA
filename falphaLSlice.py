import roo,ROOT as r
from systematics import systematics
from parabola import parabola

class falphaLSlice(object): 
    
    @staticmethod
    def columns() : return ['falphaL','minNLL','falphaT','err',
                            'profile_up','profile_down',
                            'MINUIT_cov_quality','alphaT','alphaL','d_qq','R_ag',
                            'f_qq','f_qg','f_gg','f_ag','d_xs_tt','d_xs_wj','expect_el_mj','expect_mu_mj',
                            'Ac_y_ttqq','Ac_y_ttqg','Ac_phi_ttqq','Ac_phi_ttqg','Ac_phi_ttag',
                            'A','B','C'
                            ]

    def __str__(self) :
        return '\t'.join('%+.4f'%self.results[item] if item in self.results else ' 0.0000' for item in self.columns())

    @staticmethod
    def getFitValue(parName, fit, error=False) :
        args = fit.floatParsFinal()
        i = args.index(parName)
        return args.at(i).getVal() if not error else (args.at(i).getVal(), args.at(i).getError())


    def __init__(self, workspace, fitArgs, varToProfile = 'falphaT' ) :
        res = {'falphaL': workspace.arg('falphaL').getVal()}

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
        for item in ['alphaL','alphaT','f_qq','f_qg','f_ag','f_gg',
                     'Ac_y_ttqq','Ac_y_ttqg','Ac_phi_ttqq','Ac_phi_ttqg','Ac_phi_ttag',
                     'd_qq','R_ag','d_xs_tt','d_xs_wj','expect_el_mj','expect_mu_mj'
                     ] : res[item] = workspace.arg(item).getVal()

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
