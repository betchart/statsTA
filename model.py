import sys,roo,inputs,ROOT as r

class topAsymmModel(object) :
    @roo.quiet
    def __init__(self, w = None) :
        self.observables = ['queuedbins','tridiscr']
        self.channels = dict((lepton,inputs.channel_data(lepton)) for lepton in ['el','mu'])
        self.ttcomps = ('qq','ag','gg','qg')

        if not w : w = r.RooWorkspace('Workspace')
        init_sequence = ['fractions','constraints','efficiencies','shapes','model']
        for item in init_sequence :
            print item,
            sys.stdout.flush()
            getattr(self, 'import_'+item)(w)
            print '!'
        [w.var(v).setBins( getattr( self.channels['el'].samples['data'].datas[0], 'GetNbins'+X)() ) for v,X in zip(self.observables,'XY')]

        self.w = w

    def import_fractions(self,w) :
        [roo.wimport(w, r.RooConstVar("f_%s_hat"%comp,"#hat{f}_{%s}"%comp, self.channels['el'].samples['tt'+comp].frac)) for comp in self.ttcomps]
        roo.factory(w, "R_ag[%f,0.07,1]"%(self.channels['el'].samples['ttag'].frac/self.channels['el'].samples['ttqq'].frac) )
        roo.factory(w, "d_qq[-0.9999999,1]")
        roo.factory(w, "expr::f_qq('(1+@0)*@1',{d_qq,f_qq_hat})")
        roo.factory(w, "prod::f_ag(R_ag,f_qq)")
        roo.factory(w, "expr::f_qg('(1-@0-@1)/(1+@2*@3*@4/(@5*@6))',{f_qq,f_ag,R_ag,f_gg_hat,f_qq_hat,f_ag_hat,f_qg_hat})")
        roo.factory(w, "expr::f_gg('1-@0-@1-@2',{f_qq,f_ag,f_qg})")

    def import_constraints(self,w) :
        roo.factory(w, "Gaussian::lumi_constraint( d_lumi[0,-0.2,0.2], 0, %f)"%self.channels['el'].lumi_sigma)
        for channel in self.channels.values() :
            roo.wimport_const( w, 'lumi_%s_hat'%channel.lepton, channel.lumi )
            roo.factory(w, "expr::lumi_%s('(1+@0)*@1', {d_lumi, lumi_%s_hat})"%(channel.lepton,channel.lepton))

        xs_constraints = dict([(samp[:2],(data.xs,data.xs_sigma)) for samp,data in self.channels['el'].samples.items() if data.xs>0])

        for sample,(xs,delta) in xs_constraints.items() :
            roo.wimport_const(w, 'xs_%s_hat'%sample, xs)
            roo.factory(w, "Gaussian::xs_%s_constraint( d_xs_%s[0,-1,2], 0, %f)"%(sample,sample,delta))
            roo.factory(w, "expr::xs_%s('(1+@0)*@1',{d_xs_%s, xs_%s_hat})"%(sample,sample,sample))

        roo.factory(w, "Gaussian::alphaT_constraint( alphaT[1,0.5,1.5], 1, 0.1)")
        roo.factory(w, "PROD::constraints(%s)"%', '.join([s+'_constraint' for s in ['lumi','alphaT']+['xs_'+t for t in xs_constraints]]))

        [roo.factory(w, "prod::xs_tt%s(f_%s,xs_tt)"%(comp,comp)) for comp in self.ttcomps]
        roo.wimport_const(w, 'xs_qcd', 100)

    def import_efficiencies(self,w) :
        [roo.wimport(w,
                     r.RooConstVar(*( 2*['eff_%s_%s'%(channel.lepton,sample)] + [data.eff] ) ) if data.eff>0 else
                     r.RooRealVar (*( 2*['eff_%s_%s'%(channel.lepton,sample)] + [0.01,0,1]) ))
         for channel in self.channels.values()
         for sample,data in channel.samples.items()
         if sample!='data']

    def import_shapes(self,w) :
        roo.factory(w, "channel[%s]"%','.join("%s=%d"%(s,i) for i,s in enumerate(self.channels)))
        [roo.factory(w, "%s[-1,1]"%obs) for obs in self.observables]

        [self.import_shape(w,channel.lepton,sample,data)
         for channel in self.channels.values()
         for sample,data in channel.samples.items()
         if sample!='data' ]

    @roo.quiet
    def import_shape(self,w,lepton,sample,data) :
        name = '_'.join([lepton,sample])
        arglist = r.RooArgList(*[w.var(o) for o in self.observables])
        argset = r.RooArgSet(arglist)

        for i,label in enumerate(['both','symm']) :
            roo.wimport(w, r.RooDataHist('%s_sim_%s'%(name,label), '', arglist, data.datas[i]))
            roo.wimport(w, r.RooHistPdf('%s_%s'%(name,label),'', argset, w.data('%s_sim_%s'%(name,label))))

        roo.factory(w, "prod::expect_%s(%s)"%(name, ','.join(['lumi_'+lepton,
                                                              'xs_'+sample,
                                                              'eff_'+name])))

    def import_model(self,w) :
        roo.factory(w, "alphaL[1, -15, 20]")

        [roo.factory(w, "SUM::%(n)s( alphaT * %(n)s_both, %(n)s_symm )"%{'n':lepton+'_ttag'}) for lepton in self.channels]
        [roo.factory(w, "SUM::%(n)s( alphaT * %(n)s_both, %(n)s_symm )"%{'n':lepton+'_ttqg'}) for lepton in self.channels]
        [roo.factory(w, "SUM::%(n)s( alphaL * %(n)s_both, %(n)s_symm )"%{'n':lepton+'_ttqq'}) for lepton in self.channels]

        [roo.factory(w, "SUM::model_%s( %s )"%(lepton, ','.join([ 'expect_%s_%s * %s_%s%s'%(lepton, key, lepton, key, value)
                                                                  for key,value in {'qcd' :'_symm',
                                                                                    'dy'  :'_symm',
                                                                                    'wj'  :'_both',
                                                                                    'st'  :'_both',
                                                                                    'ttgg':'_both',
                                                                                    'ttag':'',
                                                                                    'ttqg':'',
                                                                                    'ttqq':''}.items()
                                                                  ]) ) ) for lepton in self.channels]

        roo.factory(w, "SIMUL::model(channel, %s)"%', '.join("%(chan)s=model_%(chan)s"%{'chan':lepton} for lepton in self.channels))
        roo.factory(w, "PROD::constrained_model( model, constraints )")
