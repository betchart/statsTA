import inputs
import ROOT as r
r.gROOT.SetBatch(1)
r.gStyle.SetOptStat(0)

dists = ['fitTopPtOverSumPt_triD','fitTopTanhRapiditySum_triD','fitTopQueuedBin7TridiscriminantWTopQCD']
processes = ['ttgg','ttqg','ttqq','ttag']
colors = [r.kBlack,r.kRed,r.kBlue,r.kViolet]
fileName = 'graphics/fractions.pdf'

can = r.TCanvas()
can.Divide(2,2)
can.Print(fileName+'[')

for dist in dists :
    xname = dist.split('_')[0].replace('fitTopQueuedBin7TridiscriminantWTopQCD','unrolling of X_T by X_L [reco]').replace('fitTopTanhRapiditySum','tanh|t#bar{t}.y| [reco]').replace('fitTopTanhAvgRapidity','tanh(|t.y+#bar{t}.y|/2) [reco]').replace('fitTopPtOverSumPt','t#bar{t}.pt / (t.pt + #bar{t}.pt) [reco]')
    channels = dict((lepton,inputs.channel_data(lepton,'top',signal=dist,getTT=True)) for lepton in ['el','mu'])
    for i,(lep,c) in enumerate(channels.items()) : 
        titlename = {'el':'Electrons+Jets','mu':'Muons+Jets'}[lep]
        can.cd(i+1)
        denom = c.samples['tt'].datasX[0]
        denom.SetMinimum(0)
        denom.SetTitle('%s;%s;%s'%(titlename,xname,'selected t#bar{t} events / bin / 19.5 / fb'))
        denom.GetYaxis().SetTitleOffset(1.3)
        denom.Draw('hist')
        can.cd(i+3)
        hists = [c.samples[p].datasX[0] for p in processes]
        for h in hists :
            h.Divide(denom)
            h.SetMinimum(0)
            h.SetMaximum(1.0)
            h.SetTitle(';%s;%s'%(xname,'production fraction'))
        for j,color,h in zip(range(len(colors)),colors,hists) :
            h.SetLineColor(color)
            h.Draw('histsame' if j else 'hist')
    can.Print(fileName)
    
can.Print(fileName+']')
