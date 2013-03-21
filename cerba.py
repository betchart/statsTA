
import os,inputs,ROOT as r

os.system('mkdir -p cerba')
txtname = 'cerba/cerba_efficiencies.txt'
with open(txtname,'w') as txt:
    print>>txt, '# sample    efficiency-5j/4j      efficiency3dcut'
    for lep in ['el','mu'] :
        chan = inputs.channel_data(lep, 'top', signal='jetMoments2Sum_triD')
        chan4 = inputs.channel_data(lep, 'top', signal='fitTopTanhRapiditySum_triD')

        print>>txt
        print>>txt, '# ', lep
        filename = "cerba/cerba_%s.root"%lep
        f = r.TFile.Open(filename,'RECREATE')
        for item in chan.samples:
            eff5jet = chan.samples[item].datas[0].Integral()/chan4.samples[item].datas[0].Integral()
            h = chan.samples[item].datas[0].ProjectionX(item,3,3)
            eff3d = h.Integral() / chan.samples[item].datas[0].Integral()
            if item!='data': h.Scale(1./h.Integral())
            h.Write()
            print>>txt, '\t'.join([str(i) for i in [item, eff5jet,eff3d]])
        f.Close()
        print "Wrote", filename
print "Wrote", txtname
