import roo,ROOT as r
from systematics import systematics
from parabola import parabola

class dqqSlice(object): 
    
    @staticmethod
    def columns() : return ['dqq','minNLL','alphaL','alphaL_err','alphaL_2err',
                            'alphaL_profile_up','alphaL_profile_down',
                            'sys_d_lumi_up','sys_d_lumi_down',
                            'sys_d_xs_dy_up','sys_d_xs_dy_down',
                            'sys_d_xs_st_up','sys_d_xs_st_down',
                            'MINUIT_cov_quality','alphaT','R_ag'
                            ]

    def __str__(self) :
        return '\t'.join('%+.4f'%self.results[item] if item in self.results else ' 0.0000' for item in self.columns())


    @staticmethod
    def getFitValue(parName, fit, error=False) :
        args = fit.floatParsFinal()
        i = args.index(parName)
        return args.at(i).getVal() if not error else (args.at(i).getVal(), args.at(i).getError())


    def __init__(self, dqq, workspace, fitArgs ) :
        res = {'dqq':dqq}

        model = workspace.pdf('model')
        data = workspace.data('data')
        workspace.arg('d_qq').setVal(dqq)

        central = model.fitTo( data, *(fitArgs+[r.RooFit.Save(True)]))
        nll = model.createNLL(data, *fitArgs[:-1] )
        pll = nll.createProfile(workspace.argSet('alphaL'))

        res['MINUIT_cov_quality'] = central.covQual()

        aL,aLe = self.getFitValue('alphaL',central,error=True)
        aLarg = workspace.arg('alphaL')

        def rooeval(yvar,x) :
            aLarg.setVal(x)
            return yvar.getVal()

        prbl = parabola([(x,rooeval(pll,x)) for x in [aL+aLe,aL-aLe,aL]])
        pll_1sigma = 0.5
        pll_2sigma = 2.0
        aLarg.setVal(prbl.xmin)

        res['minNLL'] = nll.getVal()
        res['alphaL'] = prbl.xmin
        res['alphaL_err'] = prbl.dx(pll_1sigma)
        res['alphaL_2err'] = prbl.dx(pll_2sigma)
        pll.getVal()
        for item in ['alphaT','R_ag'] : res[item] = workspace.arg(item).getVal()

        for sign,lab in zip([-1,1],['down','up']) :
            xTrial = res['alphaL']+sign*res['alphaL_err']
            dY = rooeval(pll,xTrial) - pll_1sigma
            dX = -2*dY / prbl.slope(xTrial)
            bounds = sorted([xTrial, xTrial+dX])
            res['alphaL_profile_%s'%lab] = xTrial if abs(dY)<2e-5 else pll.findRoot(aLarg, bounds[0], bounds[1], pll_1sigma)

        pll.IsA().Destructor(pll)
        nll.IsA().Destructor(nll)

        for s,(nom,err) in systematics.items() :
            for sign,label in [(-1,'down'),(+1,'up')] :
                workspace.arg(s).setVal(nom+sign*err)
                nll = model.createNLL(data, *fitArgs[:-1] )
                pll = nll.createProfile(workspace.argSet('alphaL'))
                points = [(x,rooeval(pll,x+0.01)) for x in [aL+aLe,aL-aLe,aL]]
                pll.IsA().Destructor(pll)
                nll.IsA().Destructor(nll)
                res['sys_%s_%s'%(s,label)] = parabola(points).xmin

            workspace.arg(s).setVal(nom)

        del central

        self.results = res
