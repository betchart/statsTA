import ROOT as r
r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L tdrstyle.C")
r.setTDRStyle()
r.gStyle.SetOptFit(0)

xstat = 0.36828715173
ystat = 0.108619583584
nll = -6121051.71045

names = ['ensemble_atMeasured','ensemble_B_LF100','ensemble_C_LF100','ensemble_atZero']
tfiles = [r.TFile.Open('data/%s.root'%name) for name in names]

c = r.TCanvas()
c.SetRightMargin(0.06)
outName = 'ensembles.pdf'
c.Print(outName + '[')

plots = ['delta_Aqq','delta_Aqg','error_Aqq','error_Aqg','pullqq','pullqg','dNLL']
labels = ['#delta A_{c}^{y(q#bar{q})} (%)',
          '#delta A_{c}^{y(qg)} (%)',
          '#sigma A_{c}^{y(q#bar{q})} (%)',
          '#sigma A_{c}^{y(qg)} (%)',
          '#frac{#delta}{#sigma}  A_{c}^{y(q#bar{q})} ',
          '#frac{#delta}{#sigma}  A_{c}^{y(qg)} ',
          'NLL-#bar{NLL}',
]

a = r.TArrow()
a.SetLineColor(r.kBlue)
a.SetLineWidth(4)

def par(h,name):
    p = h.GetFunction('gaus').GetParameter(name)
    e = h.GetFunction('gaus').GetParError(h.GetFunction('gaus').GetParNumber(name))
    p1,p2 = divmod(abs(p),1)
    sign = '-' if p<0 else ''
    s = 100000
    return "%s%d&%05d(%05d)"%(sign,p1,s*p2,s*e)
        
print names
print
for item,x in zip(plots,labels):
    print item
    h = [tfile.Get(item) for tfile in tfiles]
    m = max([1.2 * h_.GetMaximum() for h_ in h])
    if 'NLL' in item:
        r.TGaxis.SetMaxDigits(3)
        c.SetRightMargin(0.1)
    for h_,n in zip(h,names): 
        h_.UseCurrentStyle()
        h_.GetYaxis().SetTitle('Pseudo-Experiments')
        h_.GetXaxis().SetTitle(x)
        h_.SetMaximum(m)
        if 'NLL' not in item: print par(h_,'Mean'), '&', par(h_,'Sigma')
    print
    h[0].Draw()
    for i,(color,style) in enumerate(zip([r.kRed,r.kRed,r.kBlack],[False,True,True])):
        h[i+1].SetLineColor(color)
        h[i+1].SetMarkerColor(color)
        h[i+1].SetMarkerSize(1.1)
        if style: h[i+1].SetMarkerStyle(4)
        if 'NLL' not in item: 
            h[i+1].GetFunction('gaus').SetLineColor(color)
            if style: h[i+1].GetFunction('gaus').SetLineStyle(2)
        h[i+1].Draw('same')

    x_ = xstat if 'Aqq' in item else ystat if 'Aqg' in item else nll-tfiles[0].Get('meanNLL').GetMean()
    if 'error' in item or 'NLL' in item: a.DrawArrow(x_,3,x_,0.4*h[0].GetMaximum(),0.05,"<")
    c.Update()
    c.Print(outName)

c.Print(outName + ']')