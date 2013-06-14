import ROOT as r
r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L tdrstyle.C")
r.setTDRStyle()
r.gStyle.SetOptFit(0)

xstat = 0.36828715173
ystat = 0.108619583584

names = ['ensemble_atMeasured','ensemble_atZero']
tfiles = [r.TFile.Open('data/%s.root'%name) for name in names]

c = r.TCanvas()
outName = 'ensembles.pdf'
c.Print(outName + '[')

plots = ['delta_Aqq','delta_Aqg','error_Aqq','error_Aqg','pullqq','pullqg']
labels = ['#delta A_{c}^{y(q#bar{q})} (%)',
          '#delta A_{c}^{y(qg)} (%)',
          '#sigma A_{c}^{y(q#bar{q})} (%)',
          '#sigma A_{c}^{y(qg)} (%)',
          '#frac{#delta}{#sigma}  A_{c}^{y(q#bar{q})} ',
          '#frac{#delta}{#sigma}  A_{c}^{y(qg)} ',
]

a = r.TArrow()
a.SetLineColor(r.kBlue)
a.SetLineWidth(4)

def par(h,name):
    p = h.GetFunction('gaus').GetParameter(name)
    e = h.GetFunction('gaus').GetParError(h.GetFunction('gaus').GetParNumber(name))
    p1,p2 = divmod(abs(p),1)
    sign = -1 if p<0 else 1
    s = 100000
    return "%d&%05d(%05d)"%(sign*p1,s*p2,s*e)
        
print names
print
for item,x in zip(plots,labels):
    print item
    h = [tfile.Get(item) for tfile in tfiles]
    m = max([1.2 * h_.GetMaximum() for h_ in h])
    for h_,n in zip(h,names): 
        h_.UseCurrentStyle()
        h_.GetYaxis().SetTitle('Pseudo-Experiments')
        h_.GetXaxis().SetTitle(x)
        h_.SetMaximum(m)
        print par(h_,'Mean'), '&', par(h_,'Sigma')
    print
    h[1].SetLineColor(r.kRed)
    h[1].SetMarkerColor(r.kRed)
    h[1].GetFunction('gaus').SetLineColor(r.kRed)
    h[0].Draw()
    h[1].Draw('same')

    x_ = xstat if 'Aqq' in item else ystat
    if 'error' in item: a.DrawArrow(x_,3,x_,150,0.05,"<")
    c.Update()
    c.Print(outName)

c.Print(outName + ']')
