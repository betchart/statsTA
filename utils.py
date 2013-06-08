from multiprocessing import Process,JoinableQueue
import math,traceback,sys,itertools, ROOT as r
try: import numpy as np
except: pass

#####################################
hyphens="-"*115
#####################################
def get(d,key) :
    dirs = key.split('/')[:-1]
    key = key.split('/')[-1]
    try:
        for d_ in dirs: d = d.GetKey(d_).ReadObj()
        obj = d.GetKey(key).ReadObj()
    except: obj = d.Get(key)
    return obj
#####################################
def asymmetry(hist) :
    bins = [(hist.GetBinCenter(i),hist.GetBinContent(i),hist.GetBinError(i)) for i in range(hist.GetNbinsX()+2)]
    bans = [(y2-y1, y2+y1, e1**2+e2**2) for (x1,y1,e1),(x2,y2,e2) in zip(bins,bins[::-1]) if x1<x2 and y2+y1]
    num = sum( n for n,_,_ in bans )
    err2 = sum( e2 for _,_,e2 in bans )
    denom = sum( d for _,d,_ in bans )
    return  num / denom, math.sqrt(err2) / denom
#####################################
def alphaMax(symm,anti):
    tmp = symm.Clone('tmp')
    tmp.Divide(anti)
    return min(filter(None,[abs(tmp.GetBinContent(iX,iY))
                            for iX in range(1,1+tmp.GetNbinsX())
                            for iY in range(1,1+tmp.GetNbinsY())]))
#####################################
def tCanvasPrintPdf(canvas, fileName, verbose=True, option='', title=""):
    illegal = [':','[',']','(',')']
    for ill in illegal : fileName = fileName.replace(ill,"_")
    canvas.Print(fileName+".pdf" + option,"pdf"+((' Title:'+title) if title else ''))
    if verbose : print "Output file: %s.pdf"%fileName
#####################################
def operateOnListUsingQueue(nCores,workerFunc,inList,daemon=True) :
    q = JoinableQueue()
    listOfProcesses=[]
    for i in range(nCores):
        p = Process(target = workerFunc, args = (q,))
        p.daemon = daemon
        p.start()
        listOfProcesses.append(p)
    map(q.put,inList)
    q.join()# block until all tasks are done
    #clean up
    for process in listOfProcesses :
        process.terminate()
#####################################
class qWorker(object) :
    def __init__(self,func = None) : self.func = func
    def __call__(self,q) :
        while True:
            item = q.get()
            try:
                if self.func : self.func(*item)
                else: item()
            except Exception as e:
                traceback.print_tb(sys.exc_info()[2], limit=20, file=sys.stdout)
                print e.__class__.__name__,":", e
            q.task_done()
#####################################        
def dependence(TH2, name="", limit=1, inSigma = False) :
    if not TH2: return None
    if TH2.GetDirectory() : TH2.GetDirectory().cd()
    dep = TH2.Clone(name if name else TH2.GetName()+"_dependence")
    dep.SetZTitle("dependence")
    norm = TH2.Integral()
    projX = TH2.ProjectionX()
    projY = TH2.ProjectionY()
    for iX in range(1,TH2.GetNbinsX()+1) :
        for iY in range(1,TH2.GetNbinsY()+1) :
            X = projX.GetBinContent(iX)
            Y = projY.GetBinContent(iY)
            bin = TH2.GetBin(iX,iY)
            XY = max(0,TH2.GetBinContent(bin))
            dBin = math.log(norm*XY/X/Y) if XY else 0
            eX = projX.GetBinError(iX)
            eY = projX.GetBinError(iY)
            eXY = TH2.GetBinError(bin)
            eBin = math.sqrt((eXY/XY)**2 + (eX/X)**2 + (eY/Y)**2) if XY else 1# faulty assumption of independent errors
            #dep.SetBinContent(bin, min(maximum,max(minimum,math.log(norm*XY/X/Y)/(eBin if eBin and inSigma else 1))) if XY else 0)
            dep.SetBinContent(bin, min(limit,max(-limit,dBin/(eBin if eBin and inSigma else 1))) if XY else 0)
            dep.SetBinError(bin,0)
    dep.SetMinimum(-limit)
    dep.SetMaximum(limit)
    
    return dep
#####################################
def intFromBits(bits) :
    return sum([j[1] * (1<<j[0]) for j in enumerate(reversed(bits))])
#####################################
def splitList(List,item) :
    if item not in List: return [List]
    i = List.index(item)
    return [List[:i+1]] + splitList(List[i+1:],item)
#####################################
def pages(blocks,size) :
    iBreak = next((i-1 for i in range(len(blocks)) if len(sum(blocks[:i],[]))>size), None)
    return [blocks] if not iBreak else [blocks[:iBreak]] + pages(blocks[iBreak:],size)
#####################################
def justNameTitle(tkey) :
    name,title = tkey.GetName(),tkey.GetTitle()
    name = name.replace('-SLASH-','/').replace('-SEMI-',';').replace('-COLON-',':')
    L = len(title)
    return ( (name,"") if name == title else
             (name[:-L],title) if name[-L:] == title else
             (name,title) )
#####################################
def symmAnti(hist) :
    nbins = hist.GetNbinsX()
    symm = hist.Clone(hist.GetName()+'_symm')
    anti = hist.Clone(hist.GetName()+'_anti')
    for j in [0] if hist.GetDimension()==1 else range(2+hist.GetNbinsY()) :
        for i in range(2+nbins) :
            a,b = hist.GetBinContent(i,j), hist.GetBinContent(1+nbins-i,j)
            e = 0.5 * math.sqrt( hist.GetBinError(i,j)**2 + hist.GetBinError(1+nbins-i,j)**2 )
            symm.SetBinContent(i,j,0.5*(a+b)) ; symm.SetBinError(i,j,e)
            anti.SetBinContent(i,j,0.5*(a-b)) ; anti.SetBinError(i,j,e)
    return symm,anti
#####################################
def unQueuedBins(queued,N,limits1, limits2 ):
    hist = ( r.TH2D('un'+queued.GetName(), '', N, limits1[0], limits1[1], N, limits2[0], limits2[1]) if queued.GetDimension() == 1 else
             r.TH3D('un'+queued.GetName(), '', N, limits1[0], limits1[1], N, limits2[0], limits2[1], queued.GetNbinsY(), queued.GetYaxis().GetBinLowEdge(1), queued.GetYaxis().GetBinLowEdge(queued.GetNbinsY()+1))
             )
    for i,bins in enumerate([(sB-j,j) if sB%2 else (j,sB-j) for sB in range(N) for j in range(sB+1)]) :
        if queued.GetDimension() == 1:
            hist.SetBinContent(1+bins[0],1+bins[1],queued.GetBinContent(i+1))
            hist.SetBinContent(N-bins[0],N-bins[1],queued.GetBinContent(N*N-i))

            hist.SetBinError(1+bins[0],1+bins[1],queued.GetBinError(i+1))
            hist.SetBinError(N-bins[0],N-bins[1],queued.GetBinError(N*N-i))
        else:
            for j in range(queued.GetNbinsY()):
                hist.SetBinContent(1+bins[0],1+bins[1],j+1,queued.GetBinContent(i+1,j+1))
                hist.SetBinContent(N-bins[0],N-bins[1],j+1,queued.GetBinContent(N*N-i,j+1))

                hist.SetBinError(1+bins[0],1+bins[1],j+1,queued.GetBinError(i+1,j+1))
                hist.SetBinError(N-bins[0],N-bins[1],j+1,queued.GetBinError(N*N-i,j+1))

    return hist


#####################################
def roundString(val, err, width=None, noSci = False, noErr = False) :
    err_digit = int(math.floor(math.log(abs(err))/math.log(10))) if err else 0
    val_digit = int(math.floor(math.log(abs(val))/math.log(10))) if val else 0
    dsp_digit = max(err_digit,val_digit)
    sci = (val_digit<-1 or err_digit>0) and not noSci

    precision = val_digit-err_digit if sci else -err_digit

    display_val = val/pow(10.,dsp_digit) if sci else val
    display_err = str(int(round(err/pow(10,err_digit))))

    while True:
        display_sci = ("e%+d"%dsp_digit) if sci else ""
        returnVal = "%.*f(%s)%s"%(precision,display_val,display_err,display_sci) if not noErr else "%.*f%s"%(precision,display_val,display_sci)
        if (not width) or len(returnVal) <= width or precision < 1: break
        else:
            display_err = "-"
            if not precision :
                display_val*=10
                dsp_digit-=1
            precision-=1
    return returnVal
#####################################
