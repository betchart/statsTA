import math,ROOT as r
from numpy.random import poisson
r.RooMsgService.instance().setGlobalKillBelow(r.RooFit.WARNING) # suppress messages

def compare() :
    '''Check that RooFit internal likelihood for binned data treats each bin with Poisson probability.'''

    # simple binned probabilty distribution
    bins = 10
    exh = r.TH1D('exh','',bins,0,bins)
    for i in range(bins) : exh.Fill(i,2*(i+1))
    exh_i = exh.Integral()
    
    # RooFit model
    bn = r.RooRealVar('bn','bn',0,bins)
    rdh = r.RooDataHist('rdh','',r.RooArgList(bn), exh)
    expected = r.RooHistPdf('expected','epdf',r.RooArgSet(bn),rdh)
    n_expect = r.RooRealVar('n_expect','N',exh_i, exh_i/2, exh_i*2)
    expected_e = r.RooExtendPdf('expected_e','#lambda', expected, n_expect)

    # pseudo-data, explicitly Poisson distributed around pdf (not really necessary)
    data = r.TH1D('data','',bins,0,bins)
    for i in range(bins) : data.SetBinContent(i+1,poisson(exh.GetBinContent(i+1)))
    rdata = r.RooDataHist('rdata','',r.RooArgList(bn), data)
    nll = expected_e.createNLL(rdata)

    def RNLL(N) :
        '''RooFit internal negative log likelihood.'''
        n_expect.setVal(N)
        return nll.getVal()

    def NLL(N) :
        '''explicit negative log likelihood.'''
        def ll(n,lam) : return n*math.log(lam) - lam - math.log(math.factorial(n))
        return sum([-ll(data.GetBinContent(i+1), N*exh.GetBinContent(i+1)/exh_i) for i in range(bins)])
    
    nlls = [(N, NLL(N), RNLL(N)) for N in range(int(exh_i/2),int(exh_i*2))]
    differences = set(round(nll[1]-nll[2],8) for nll in nlls)
    return len(differences) == 1

if __name__=='__main__' :
    print 'Pass' if compare() else 'Fail', 'Test'
