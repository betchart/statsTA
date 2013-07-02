import ROOT as r
import sys

def get(d,key) :
    dirs = key.split('/')[:-1]
    key = key.split('/')[-1]
    try:
        for d_ in dirs: d = d.GetKey(d_).ReadObj()
        obj = d.GetKey(key).ReadObj()
    except: obj = d.Get(key)
    return obj

if len(sys.argv)<3 or not all('.root' in v for v in sys.argv[1:]):
    print "Usage: compare_files.py <file1.root> <file2.root>"
    exit()


def compareListsOfKeys(i1,i2):
    KL1 = i1.GetListOfKeys()
    KL2 = i2.GetListOfKeys()
    ks1 = set(k.GetName() for k in KL1)
    ks2 = set(k.GetName() for k in KL2)
    return ( set.intersection(ks1,ks2), 
             (ks1 & ks2) - ks1, 
             (ks1 & ks2) - ks2 )

def sameHistos(h1,h2):
    nbinsX = h1.GetNbinsX()
    nbinsY = h1.GetNbinsY()
    bins = nbinsX == h2.GetNbinsX() and nbinsY == h2.GetNbinsY()
    loX = h1.GetXaxis().GetBinLowEdge(0) == h2.GetXaxis().GetBinLowEdge(0)
    hiX = h1.GetXaxis().GetBinLowEdge(nbinsX+1) == h2.GetXaxis().GetBinLowEdge(nbinsX+1)

    pDiff = max([(abs(v1-v2)/(v1+v2) if v1+v2 else (v1+v2))
                 for v1,v2 in [( h1.GetBinContent(iX,iY),
                                 h2.GetBinContent(iX,iY) ) 
                               for iX in range(nbinsX+2) for iY in range(nbinsY+2)]])
    same = pDiff < 1e-6
    return all([bins,loX,hiX,same])
    
def checkdirs(d1,d2):
    keys,not1,not2 = compareListsOfKeys( d1, d2 )
    if not1: print "Missing from arg1:", not1
    if not2: print "Missing from arg2:", not2

    for item in sorted(keys):
        name = '/'.join([':'.join(d1.GetPath().split(':')[1:]),item])[1:]
        i1 = get(f1,name)
        i2 = get(f2,name)
        if not i1 :
            print "nil i1:", name
        elif not i2:
            print "nil i2:", name
        elif type(i1)==r.TDirectoryFile:
            checkdirs(i1,i2)
        elif not sameHistos(i1,i2):
            print 'diff:', name
        del i1
        del i2
    return

f1 = r.TFile.Open(sys.argv[1])
f2 = r.TFile.Open(sys.argv[2])

checkdirs(f1,f2)
print 'Done'
f1.Close() # hangs
print 'close f1'
f2.Close()
print 'close f2'
