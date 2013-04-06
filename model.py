import sys,roo,inputs,utils,ROOT as r

defaultDist = 'fitTopQueuedBin7TridiscriminantWTopQCD'

class topModel(object) :
    @roo.quiet
    def __init__(self, w = None, dist=defaultDist, asymmetry=True) :
        self.asymmetry=asymmetry
        self.observables = ['queuedbins','tridiscr']
        self.gen = inputs.channel_data('mu','top',signal='2_x_y',getTT=True, noRebin=True)
        self.channels = dict((lepton,inputs.channel_data(lepton,'top',signal=dist)) for lepton in ['el','mu'])
        self.channels_qcd = dict((lepton+'qcd',inputs.channel_data(lepton,'QCD',signal=dist)) for lepton in ['el','mu'])
        self.ttcomps = ('qq','ag','gg','qg')
        self.toSymmetrize = ['dy'] if dist==defaultDist else []

        if not w : w = r.RooWorkspace('Workspace')
        init_sequence = ['fractions','xs_lumi','efficiencies','shapes','qcd','asymmetry','model','expressions']
        for item in init_sequence :
            print item,
            sys.stdout.flush()
            getattr(self, 'import_'+item)(w)
            print '!'
        [w.var(v).setBins( getattr( self.channels['el'].samples['data'].datas[0], 'GetNbins'+X)() ) for v,X in zip(self.observables,'XY')]

        for item in ['d_lumi','d_xs_dy','d_xs_st'] : self.model.w.arg(item).setConstant()

        self.w = w

    def import_fractions(self,w) :
        [roo.wimport(w, r.RooConstVar("f_%s_hat"%comp,"#hat{f}_{%s}"%comp, self.channels['el'].samples['tt'+comp].frac)) for comp in self.ttcomps]
        roo.factory(w, "R_ag[%f,0.07,1]"%(self.channels['el'].samples['ttag'].frac/self.channels['el'].samples['ttqq'].frac) )
        roo.factory(w, "d_qq[-0.999999,1]")
        roo.factory(w, "expr::f_qq('(1+@0)*@1',{d_qq,f_qq_hat})")
        roo.factory(w, "prod::f_ag(R_ag,f_qq)")
        roo.factory(w, "expr::f_qg('(1-@0-@1)/(1+@2*@3*@4/(@5*@6))',{f_qq,f_ag,R_ag,f_gg_hat,f_qq_hat,f_ag_hat,f_qg_hat})")
        roo.factory(w, "expr::f_gg('1-@0-@1-@2',{f_qq,f_ag,f_qg})")

    def import_xs_lumi(self,w) :
        roo.factory(w, "d_lumi[0,-0.2,0.2]")
        for lepton,channel in self.channels.items() + self.channels_qcd.items():
            roo.wimport_const( w, 'lumi_%s_hat'%lepton, channel.lumi )
            roo.factory(w, "expr::lumi_%s('(1+@0)*@1', {d_lumi, lumi_%s_hat})"%(lepton,lepton))

        xs_constraints = dict([(samp[:2],(data.xs,data.xs_sigma)) for samp,data in self.channels['el'].samples.items() if data.xs>0])

        for sample,(xs,delta) in xs_constraints.items() :
            roo.wimport_const(w, 'xs_%s_hat'%sample, xs)
            roo.factory(w, "d_xs_%s[0,-1,2]"%sample)
            roo.factory(w, "expr::xs_%s('(1+@0)*@1',{d_xs_%s, xs_%s_hat})"%(sample,sample,sample))

        [roo.factory(w, "prod::xs_tt%s(f_%s,xs_tt)"%(comp,comp)) for comp in self.ttcomps]

    def import_efficiencies(self,w, channels = None) :
        if not channels : channels = self.channels
        [roo.wimport(w, r.RooConstVar(*( 2*['eff_%s_%s'%(lepton,sample)] + [data.eff] ) ))
         for lepton,channel in channels.items()
         for sample,data in channel.samples.items()
         if sample!='data']

    def import_shapes(self,w, channels = None) :
        if not channels :
            channels = self.channels
            roo.factory(w, "channel[%s]"%','.join("%s=%d"%(s,i) for i,s in enumerate(channels)))
            [roo.factory(w, "%s[-1,1]"%obs) for obs in self.observables]

        [self.import_shape(w,lepton,sample,data)
         for lepton,channel in channels.items()
         for sample,data in channel.samples.items()
         if sample!='data' ]

    @roo.quiet
    def import_shape(self,w,lepton,sample,data) :
        name = '_'.join([lepton,sample])
        arglist = r.RooArgList(*[w.var(o) for o in self.observables])
        argset = r.RooArgSet(arglist)

        for i,label in enumerate(['both','symm'][:None if sample in ['ttag','ttqg','ttqq','dy'] else -1]) :
            roo.wimport(w, r.RooDataHist('%s_sim_%s'%(name,label), '', arglist, data.datas[i]))
            roo.wimport(w, r.RooHistPdf('%s_%s'%(name,label),'', argset, w.data('%s_sim_%s'%(name,label))))

        roo.factory(w, "prod::expect_%s(%s)"%(name, ','.join(['lumi_'+lepton,
                                                              'xs_'+sample,
                                                              'eff_'+name,
                                                              '-1',
                                                              'factor_%s'%lepton][:None if 'qcd' in lepton else -2])))

    def import_qcd(self,w) :
        [roo.factory(w, "factor_%s[1,0,10]"%lepton) for lepton in self.channels_qcd]

        self.import_efficiencies(w,self.channels_qcd)
        self.import_shapes(w,self.channels_qcd)
        arglist = r.RooArgList(*[w.var(o) for o in self.observables])
        argset = r.RooArgSet(arglist)
        for lepton,channel in self.channels_qcd.items() :
            roo.wimport(w, r.RooDataHist('%s_data_sim_both'%lepton, '', arglist, channel.samples['data'].datas[0]))
            roo.wimport(w, r.RooHistPdf('%s_data_both'%lepton, '', argset, w.data('%s_data_sim_both'%lepton)))
            N = channel.samples['data'].datas[0].Integral()
            roo.factory(w, "expr::expect_%s_data('@0*%f*@1',{lumi_%s,factor_%s})"%(lepton,N/channel.lumi,lepton,lepton))

    def import_asymmetry(self,w) :
        if not self.asymmetry : return
        roo.factory(w, "falphaL[0.1, -15, 15]")
        roo.factory(w, "falphaT[0.2, -15, 15]")
        roo.factory(w, "expr::alphaL('@0/@1',{falphaL,f_qq})")
        roo.factory(w, "expr::alphaT('@0/@1',{falphaT,f_qg})")

        [roo.factory(w, "SUM::%(n)s( alphaT * %(n)s_both, %(n)s_symm )"%{'n':lepton+'_ttag'}) for lepton in self.channels.keys()+self.channels_qcd.keys()]
        [roo.factory(w, "SUM::%(n)s( alphaT * %(n)s_both, %(n)s_symm )"%{'n':lepton+'_ttqg'}) for lepton in self.channels.keys()+self.channels_qcd.keys()]
        [roo.factory(w, "SUM::%(n)s( alphaL * %(n)s_both, %(n)s_symm )"%{'n':lepton+'_ttqq'}) for lepton in self.channels.keys()+self.channels_qcd.keys()]

        assert self.gen.samples['tt'].datas[0].GetXaxis().GetTitle() == 'genTopPhiBoost'
        assert self.gen.samples['tt'].datas[0].GetYaxis().GetTitle() == 'genTopDeltaBetazRel'

        for name,data in self.gen.samples.items() :
            roo.wimport(w, r.RooConstVar(*(2*['Ac_y_'+name]+[utils.asymmetry(data.datasY[0])]))) #A_c^y(**)
            roo.wimport(w, r.RooConstVar(*(2*['Ac_phi_'+name]+[utils.asymmetry(data.datasX[0])]))) #A_c^\phi(**)
            w.arg('Ac_y_'+name).Print()
            w.arg('Ac_phi_'+name).Print()

    def import_model(self,w) :
        which = dict((i,'_both') for i in ['dy','wj','st','ttgg','ttqq','ttqg','ttag'])
        if self.asymmetry : which.update({'dy':'_symm','ttqq':'','ttqg':'','ttag':''})

        [roo.factory(w, "SUM::model_%s( expect_%sqcd_data * %sqcd_data_both,  %s  )"%(lepton, lepton, lepton,
                                                                                      ','.join([ 'expect_%s_%s * %s_%s%s'%(lepton+part, key, lepton+part, key, value)
                                                                                                 for key,value in which.items()
                                                                                                 for part in ['','qcd']
                                                                                                 if not (key=='dy' and 'qcd' in part)]) ) )
         for lepton in self.channels]

        roo.factory(w, "SIMUL::model(channel, %s)"%', '.join("%(chan)s=model_%(chan)s"%{'chan':lepton} for lepton in self.channels))


    def import_expressions(self,w) :
        [roo.factory(w, "sum::expect_%s_tt(%s)"%(lep,','.join(['expect_%s_tt%s'%(lep,f) for f in ['gg','qg','qq','ag']]))) for lep in self.channels]
        [roo.factory(w, "sum::expect_%s_notqcd(%s)"%(lep,','.join(['expect_%s_%s'%(lep,samp) for samp in ['wj','st','ttgg','ttag','ttqg','ttqq']]))) for lep in self.channels_qcd]
        [roo.factory(w, "sum::expect_%(n)s_mj(expect_%(n)sqcd_data,expect_%(n)sqcd_notqcd)"%{'n':lep}) for lep in self.channels]
