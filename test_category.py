import math,ROOT as r

def cats(constrain = True) :
    '''Test for binned category model.'''
    
    N = 100
    cats = ['f','g','h']
    xgauss = [(0.8,0.05),
              (0.2,0.05),
              (0.5,0.1)]

    w = r.RooWorkspace('Workspace')
    w.factory("chan[%s]"%', '.join('%s=%d'%(c,i) for i,c in enumerate(cats)))
    w.factory("Gaussian::constraint( N[%d,0,300], %d, 8)"%(N,N+20))
    w.factory('x[0,1]')
    [w.factory("ExtendPdf::%sE( Gaussian::%s( x, %f, %f), N)"%(cat, cat, mean, sigm)) for cat,(mean,sigm) in zip(cats,xgauss) ]
    w.factory("SIMUL::model( chan, %s )"%', '.join('%s=%sE'%(c,c) for c in cats))
    w.factory("PROD::constr_mod( model, constraint )")

    w.var('x').setBins(10)

    data = w.pdf('model').generate( r.RooArgSet(w.var('x'),w.arg('chan')),
                                    r.RooFit.AllBinned(),
                                    r.RooFit.Extended())
    fit = w.pdf('model').fitTo( *[data, 
                                  r.RooFit.Save(True),
                                  r.RooFit.PrintLevel(-1),
                                  r.RooFit.ExternalConstraints(w.argSet('constraint'))][:None if constrain else -1]
                                )
    fit.Print()

    fit2 = w.pdf('constr_mod').fitTo( *[data,
                                        r.RooFit.Save(True),
                                        r.RooFit.PrintLevel(-1),
                                        r.RooFit.Constrain(r.RooArgSet())][:-1 if constrain else None]
                                      )
    fit2.Print()

    c = r.TCanvas()
    c.Divide(len(cats))
    for i,cat in enumerate(cats) :
        c.cd(i+1)
        plt = w.var('x').frame()
        data.plotOn(plt, r.RooFit.Cut("chan==chan::%s"%cat))
        w.pdf('model').plotOn(plt,
                              r.RooFit.ProjWData(data),
                              r.RooFit.Slice(w.cat('chan'),cat)
                              )
        if i+1==len(cats) : w.pdf('model').paramOn(plt)
        plt.SetTitle("Channel: "+cat)
        plt.Draw()
    print fit.floatParsFinal()[0].getVal(), fit2.floatParsFinal()[0].getVal()
    print fit.floatParsFinal()[0].getError(), fit2.floatParsFinal()[0].getError()
    raw_input()
    

if __name__=='__main__' :
    #r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.WARNING)
    cats( constrain = False )
