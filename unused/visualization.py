



    @roo.quiet
    def visualize1D(self, canvas=None):
        r.gStyle.SetTitleX(0.2)
        r.gStyle.SetTitleY(0.98)

        if not canvas:
            canvas = r.TCanvas()
            canvas.Divide(5,2,0,0)
        
        w = self.w
        
        for j,lep in enumerate(['el','mu']):
            w.arg('observable').SetTitle('Unrolled X_{T}(7) by X_{L}(7)' if self.asymmetry else
                                         'tanh|t#bar{t}.y|')
            for i in range(5):
                canvas.cd(5*j + i + 1)
                f = w.arg('observable').frame()
                f.SetTitle('%s + jets, Tridiscriminant Bin %d'%({'el':'e','mu':'#mu'}[lep],i+1))
                f.SetLineColor(r.kWhite)
                if not i: f.SetTitleOffset(-1,'Y')
                w.arg('tridiscr').setVal(-0.8 + i*0.4)
                tridiscrBinWidth = 0.4

                mod = w.pdf('model_%s'%lep)
                args = [r.RooFit.Project(r.RooArgSet()),
                        r.RooFit.Normalization(tridiscrBinWidth, r.RooAbsReal.Relative),
                        r.RooFit.LineWidth(1),
                        ]
                dargs = [f,
                         r.RooFit.Cut('%f<=tridiscr && tridiscr<%f'%(-1+i*0.4,-0.6+i*0.4)),
                         r.RooFit.MarkerSize(0.5),
                         r.RooFit.XErrorSize(0)]
                w.data('data_%s'%lep).plotOn(*dargs)

                pf = '' if self.asymmetry else '_both'
                stack = [('_ttqq'+pf,),('_ttqg'+pf,'_ttag'+pf),('_ttgg_both',),('_wj_both',),
                         ('qcd_*',), ('_dy*','_st_both')]
                colors = [r.kViolet,r.kBlue-7,r.kBlue+2,r.kGreen+1,r.kRed,r.kGray]
                for iStack in range(len(stack)):
                    comps = lep+(','+lep).join(sum(stack[iStack:],()))
                    mod.plotOn(f,
                               r.RooFit.Components(comps),
                               r.RooFit.LineColor(colors[iStack]),
                               r.RooFit.FillColor(colors[iStack]),
                               r.RooFit.DrawOption('F'),
                               *args)

                w.data('data_%s'%lep).plotOn(*dargs)
                mod.plotOn(f, r.RooFit.LineColor(r.kOrange), *args)

                f.SetMaximum(4500 if self.asymmetry  else 4000)
                f.Draw()

        return canvas

