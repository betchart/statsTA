from inputs import channel_data
import lib
import ROOT as r
import math
r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L lib/tdrstyle.C")
r.setTDRStyle()
r.tdrStyle.SetErrorX(r.TStyle().GetErrorX())
r.tdrStyle.SetPadTopMargin(0.065)
r.TGaxis.SetMaxDigits(3)
r.tdrStyle.SetEndErrorSize(6)
#r.tdrStyle.SetPadRightMargin(0.06)

def unqueue(h):
    return lib.unQueuedBins(h,5,[-1,1],[-1,1])

threeD = True
check = False
extra = True
channels = dict([(lep, channel_data(lep, 'top', signal='fitTopQueuedBin5_TridiscriminantWTopQCD', threeD=threeD, extra=extra)) for lep in ['el', 'mu']])
linetypes = [1, 2]

comps = ['ttgg','ttag','ttqg','ttqq']
colors = [r.kBlack, r.kGreen, r.kBlue, r.kRed]

if check:
    channels_def = dict([(lep, channel_data(lep, 'top', signal='fitTopQueuedBin5_TridiscriminantWTopQCD', threeD=threeD, extra=extra)) for lep in ['el', 'mu']])
    for lep,ch in channels.items():
        for comp in comps:
            c = ch.samples[comp].datas
            d = channels_def[lep].samples[comp].datas
            for a,b in zip(c,d):
                b.Add(a,-1)
                print ' '.join(str(b.GetBinContent(i)) for i in range(2+b.GetNbinsX()))

projections = {}

for lt,(lep,ch) in zip(linetypes,channels.items()):
    for comp,color in zip(comps,colors): 
        symm = unqueue(ch.samples[comp].datas[1])
        anti = unqueue(ch.samples[comp].datas[2])
        quad = unqueue(ch.samples[comp].datas[3 if threeD else 0])
        bias = unqueue(ch.samples[comp].datas[4 if threeD else 0])
        [h.SetLineColor(color) for h in [symm,anti,quad,bias]]
        [h.SetLineStyle(lt) for h in [symm,anti,quad,bias]]
        [h.SetMarkerColor(color) for h in [symm,anti,quad,bias]]
        if lep=='mu': [h.SetMarkerStyle(4) for h in [symm,anti,quad,bias]]
        
        sf = 1./symm.Integral()
        [h.Scale(sf) for h in [symm,anti,quad,bias]]

        projections[(lep,comp)] = [(symm.ProjectionX('symmx'+lep+comp), anti.ProjectionX('antix'+lep+comp), quad.ProjectionX('quadx'+lep+comp), bias.ProjectionX('biasx'+lep+comp)),
                                   (symm.ProjectionY('symmy'+lep+comp), anti.ProjectionY('antiy'+lep+comp), quad.ProjectionY('quady'+lep+comp), bias.ProjectionY('biasy'+lep+comp)),
                                   (symm.ProjectionZ('symmz'+lep+comp), anti.ProjectionZ('antiz'+lep+comp), quad.ProjectionZ('quadz'+lep+comp), bias.ProjectionZ('biasz'+lep+comp))
                               ]

def extrema(A,B, func):
    if type(A) not in [list,tuple]:
        return func(A,B)
    return type(A)([extrema(a,b, func) for a,b in zip(A,B)])

def MAX(A,B) : return extrema(A,B, max)
def MIN(A,B) : return extrema(A,B, min)

maxs = reduce( MAX, [[tuple([i.GetMaximum()+i.GetBinError(4) for i in pair]) for pair in L] for L in projections.values()], 3*[(0,0)])
mins = reduce( MIN, [[tuple([i.GetMinimum()-i.GetBinError(4) for i in pair]) for pair in L] for L in projections.values()], 3*[(0,0)])

c = r.TCanvas()
fn = 'graphics/template.pdf'
c.Print(fn+'[')

for i,label in enumerate(['X_{L}^{rec}','X_{T}^{rec}']):
    for j,sublabel in enumerate(['symmetrized ','antisymmetrized ', 'M^{-}#bf{x}^{-} ','M^{-}#bf{x}^{+} '][:None if threeD else 2]):
        init = False
        for k,comp in enumerate(reversed(comps)):
            for lep in channels:
                h = projections[(lep,comp)][i][j]
                if j:
                    h.SetBinError(3,0)
                    h.SetMinimum(-0.007)
                    h.SetMaximum(0.007)
                else:
                    h.SetMinimum(0)
                    h.SetMaximum(0.265)
                h.GetXaxis().SetNdivisions(5,4,0,False)
                #h.SetMaximum(1.1*maxs[i][j])
                #h.SetMinimum(1.1*mins[i][j])
                h.SetMarkerSize(1.5 - 0.3*math.sqrt(k) - (0 if lep=='mu' else 0.3))
                h.SetLineWidth(4-k)
                h.GetXaxis().SetTitleOffset(0.9)
                h.GetYaxis().SetTitle(sublabel + ' probability')
                h.GetXaxis().SetTitle(label)
                h.Draw('same e1' if init else 'e1')
                init = True
        c.Print(fn)


c.Print(fn+']')

