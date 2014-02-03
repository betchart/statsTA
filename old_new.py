import sys
from inputs import channel_data
import lib
import ROOT as r
import math
from lib.__autoBook__ import autoBook
r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L tdrstyle.C")
r.setTDRStyle()
r.tdrStyle.SetErrorX(r.TStyle().GetErrorX())
r.tdrStyle.SetPadTopMargin(0.065)
r.TGaxis.SetMaxDigits(3)
r.tdrStyle.SetEndErrorSize(6)
#r.tdrStyle.SetPadRightMargin(0.06)

def unqueue(h):
    return lib.unQueuedBins(h,5,[-1,1],[-1,1])

threeD = False

comps = ['ttag','ttqg','ttqq']
colors = [r.kBlack, r.kGreen, r.kBlue, r.kRed]

projections = {}
book = autoBook("stuff")

for extra in [True,False,'only']:
    for template in [None]+range(1000):
        print template
        #sys.stdout.flush()
        channels = dict([(lep, channel_data(lep, 'top', signal='fitTopQueuedBin5_TridiscriminantWTopQCD', threeD=threeD, extra=extra, templateID=template,getTT=True)) for lep in ['el', 'mu']])
        for lep,ch in channels.items():
            for comp,color in zip(comps,colors):
                name = '_'.join(['extra' if extra=='only' else 'total' if extra else 'orig' ,lep, comp])
                if template==None: name += '_actual'
                tot = unqueue(ch.samples[comp].datas[0].ProjectionX())
                v = (100*lib.asymmetry(tot.ProjectionX())[0],
                     100*lib.asymmetry(tot.ProjectionY())[0])
                book.fill( v, name, (100,100), (-3,-3), (3,3))

tfile = r.TFile.Open("jiggled_asymmetries.root","RECREATE")
for key,hist in book.items():
    hist.Write()
tfile.Close()
