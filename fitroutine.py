import math
import array
import numpy as np
import ROOT as r
import roo
import inputs
import model
from enclosing_ellipse import enclosing_ellipse
from asymmNames import genNameX,genNameY


class fit(object):
    def __init__(self, label, signal, profileVars, R0_,
                 d_lumi, d_xs_dy, d_xs_st, tag, genPre, sigPre, dirIncrement, genDirPre, 
                 quiet = False, hackZeroBins=False, defaults = {},
                 log=None, fixSM=False,altData=None, lumiFactor=1.0):

        for item in ['label','quiet','fixSM','profileVars'] : setattr(self,item,eval(item))
        self.log = log if log else sys.stdout
        if type(R0_) == tuple:
            diffR0_ = R0_[1]
            R0_ = R0_[0]
        else: diffR0_ = None
        prePre = dirIncrement in [0,4,5]
        channels = dict([((lep,part),
                          inputs.channel_data(lep, part, tag, signal, sigPre, 
                                              "R%02d" % (R0_ + dirIncrement),
                                              genDirPre, prePre = prePre,
                                              hackZeroBins=hackZeroBins and 'QCD'==part))
                         for lep in ['el', 'mu']
                         for part in ['top', 'QCD']
                         ])
        channels['gen'] = inputs.channel_data('mu', 'top', tag,
                                              '%s; %s'%(genNameX,genNameY),
                                              sigPrefix = sigPre if dirIncrement in [0,4,5] else '',
                                              dirPrefix=genDirPre, genDirPre=genDirPre, 
                                              getTT=True, prePre = prePre)

        #for n,c in channels.items() :
        #    if n == 'gen': continue
        #    print n
        #    print c.asymmStr()
        #    print

        if diffR0_ :
            for lepPart,chan in channels.items():
                if type(lepPart) != tuple: continue
                lep,part = lepPart
                chan.subtract(inputs.channel_data(lep,part,tag,signal,sigPre,
                                                  "R%02d" % (diffR0_ + dirIncrement),
                                                  genDirPre, prePre = prePre ))

        print "###", label
        print>>self.log, "###", label
        self.model = model.topModel(channels, asymmetry=True)
        self.model.w.arg('lumi_factor').setVal(lumiFactor)
        for k,v in defaults.items(): self.model.w.arg(k).setVal(v)
        for item in ['d_lumi', 'd_xs_dy', 'd_xs_st']: self.model.w.arg(item).setVal(eval(item))

        self.fitArgs = [r.RooFit.Extended(True), r.RooFit.NumCPU(1),
                        r.RooFit.PrintLevel(-1)]
        self.model.import_data(altData)

        if fixSM: self.doSM()
        else: self.doFit()

    @roo.quiet
    def doSM(self):
        fixVars = ['R_ag','slosh','falphaL','falphaT']
        for item in fixVars: self.model.w.arg(item).setConstant()
        nll = self.model.w.pdf('model').createNLL(self.model.w.data('data'), *self.fitArgs[:-1])
        minu = r.RooMinuit(nll)
        minu.setPrintLevel(-1)
        minu.setNoWarn()
        minu.setStrategy(2)
        minu.migrad()

    @roo.quiet
    def doFit(self):
        w = self.model.w
        contourPoints=None
        while not contourPoints:
            pll = self.minimize()
            p0, contourPoints = self.contourPoints(pll)

        ell = enclosing_ellipse([p[:2] for p in contourPoints],p0[:2])
        self.profVal = p0[:2]
        self.profErr = ell.sigmas2
        self.profPLL = p0[2]

        self.scales = np.array([w.arg(a).getVal() for a in ['Ac_y_ttqq', 'Ac_y_ttqg']])
        self.scalesPhi= [w.arg('Ac_phi_%s'%n).getVal() for n in ['ttqq','ttgg','ttag','ttqg','tt']]
        self.correction = w.arg('Ac_y_ttgg').getVal() * w.arg('f_gg').getVal()
        self.fractionHats = [w.arg('f_%s_hat' % a).getVal() for a in ['gg','qg','qq','ag']]
        self.fitX,self.fitY = [float(i) for i in self.profVal*self.scales]
        self.fitXX = float(self.profErr[0,0] * self.scales[0]**2)
        self.fitXY = float(self.profErr[0,1] * self.scales[0]*self.scales[1])
        self.fitYY = float(self.profErr[1,1] * self.scales[1]**2)
        self.contourPointsX,self.contourPointsY, = zip(*[[float(i) for i in self.scales*p] for p in contourPoints])

        if not self.quiet:
            print>>self.log, self.profVal, self.profPLL
            print>>self.log, self.profErr
            for item in ['d_qq','d_xs_tt','d_xs_wj',
                         'factor_elqcd','factor_muqcd',
                         'f_gg','f_qq','f_qg','f_ag',
                         'R_ag','slosh','alphaL']:
                print>>self.log, '\t', roo.str(w.arg(item))
        self.pll = pll
        return

    @roo.quiet
    def minimize(self):
        w = self.model.w
        nll = w.pdf('model').createNLL(w.data('data'), *self.fitArgs[:-1])
        minu = r.RooMinuit(nll)
        minu.setPrintLevel(-1)
        minu.setNoWarn()
        for j in range(10):
            minu.setStrategy(2)
            for i in range(10):
                self.fitstatus = minu.migrad()
                print>>self.log, i + 10*j,
                self.log.flush()
                if not self.fitstatus: break
            if not self.fitstatus: break
            minu.setStrategy(1)
            minu.migrad()
        print>>self.log
        self.NLL = nll.getVal()
        pll = nll.createProfile(w.argSet(','.join(self.profileVars)))
        print>>self.log, roo.str(nll)
        print>>self.log, roo.str(pll)
        pll.minimizer().setStrategy(2)
        return pll
        

    def contourPoints(self,pll):
        w = self.model.w
        p0 = [w.arg(a).getVal() for a in self.profileVars]
        errMin = 0.12
        oneSigmaNLL = 1.14
        pllPoints = [(p0[0]+errMin*math.cos(t),p0[1]+errMin*math.sin(t))
                     for t in np.arange(0,2*math.pi,math.pi/8)]
        pllCache = {}
        def pllEval(p):
            p = tuple(p)[:2]
            if p not in pllCache:
                for name,val in zip(self.profileVars,p): w.arg(name).setVal(val)
                pllCache[p] = pll.getVal()
            return pllCache[p]

        p0 += pllEval(p0),        
        targetPLL = oneSigmaNLL + p0[2]

        def contourIntersect(guess,epsilon=0.01,xepsilon=0.001):
            guess += pllEval(guess),
            def point(g): return (p0[0] + g * (guess[0]-p0[0]),
                                  p0[1] + g * (guess[1]-p0[1]))

            def bsearch(lo,hi):
                assert pllEval(lo) < targetPLL
                assert pllEval(hi) > targetPLL
                p = tuple(0.5*(lo[i]+hi[i]) for i in [0,1])
                PLL = pllEval(p)
                if math.sqrt(sum((lo[i]-hi[i])**2 for i in [0,1])) < xepsilon: return p
                if PLL < targetPLL-epsilon: return bsearch(p,hi)
                if PLL > targetPLL+epsilon: return bsearch(lo,p)
                return p

            mlo = guess if guess[2] < targetPLL else None
            mhi = guess if targetPLL < guess[2] else None
            iteration = 0
            while True:
                if guess[2] < p0[2] : return None
                p = point( math.sqrt((targetPLL-p0[2]) / (guess[2]-p0[2])) )
                if abs(targetPLL-pllEval(p)) < epsilon: return p
                if pllEval(p) < targetPLL and (not mlo or pllEval(mlo) < pllEval(p)): mlo = p
                if pllEval(p) > targetPLL and (not mhi or pllEval(p) < pllEval(mhi)): mhi = p
                if iteration>5 and mlo and mhi: return bsearch(mlo,mhi)
                guess = p + (pllEval(p),)
                iteration += 1

        points = []
        for p in pllPoints:
            points.append(contourIntersect(p))
            if not points[-1]: return None,None
        return p0,points


    @staticmethod
    def fields(): 
        return ('#label fqq.Ac_y_qq  fqg.Ac_y_qg  XX  XY  YY fhat_gg ' +
                'fhat_qg fhat_qq fhat_ag Ac_y_qq_hat Ac_y_qg_hat f_gg.Ac_y_gg fitstatus NLL Ac_phi_qq_hat Ac_phi_gg_hat Ac_phi_ag_hat Ac_phi_qg_hat')

    def __str__(self):
        return '\t'.join(str(i) for i in [self.label] +
                         list(self.profVal*self.scales) +
                         [self.profErr[0,0]*self.scales[0]**2,
                          self.profErr[0,1]*self.scales[0]*self.scales[1],
                          self.profErr[1,1]*self.scales[1]**2] +
                         self.fractionHats+
                         list(self.scales) +
                         [self.correction,
                          self.fitstatus,
                          self.NLL]+
                         self.scalesPhi[:4]
                         )

    def ttree(self):
        # Note : ROOT and array.array use opposite conventions for upper/lowercase (un)signed
        #         name     array  ROOT  ROOT_typedef
        types = {int    : ("i", "I", "Int_t"),
                 long   : ("l", "L", "Long_t"),
                 float  : ("f", "F", "Float_t"),
                 bool   : ("B", "O", "Bool_t"),
                 str    : ("c", "C", "Char_t")
             }

        tree = r.TTree('fitresult','')

        selfStuff = ['label','fitX','fitY','fitXX','fitXY','fitYY','NLL','fitstatus','contourPointsX','contourPointsY','correction','fixSM']
        modelStuff = [item%xx for item in ['Ac_y_tt%s','Ac_phi_tt%s','f_%s_hat','f_%s'] for xx in ['qq','qg','ag','gg']]
        modelStuff += ['d_xs_%s'%item for item in ['tt','wj','st','dy']]
        modelStuff += ['expect_%s_%s'%(lep,s) for lep in ['el','mu'] for s in ['tt','wj','mj','st','dy']]
        modelStuff += ['d_lumi','lumi_factor','R_ag','slosh','alphaL','alphaT','falphaL','falphaT','factor_elqcd','factor_muqcd']

        selfPairs = [(item,getattr(self,item)) for item in selfStuff]
        modelPairs = [(item,self.model.w.arg(item).getVal()) for item in modelStuff]
        
        address = {}
        for name,value in selfPairs+modelPairs:
            if type(value) not in [list,tuple]:
                ar,ro,t = types[type(value)]
                address[name] = array.array(ar, value if type(value)==str else [value])
                tree.Branch(name, address[name], "%s/%s"%(name,ro))
            else:
                address[name] = r.std.vector('float')()
                for item in value: address[name].push_back(item)
                tree.Branch(name, address[name])
        tree.Fill()
        return tree
