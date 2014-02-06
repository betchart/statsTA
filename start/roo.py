import ROOT as r
r.gROOT.SetStyle("Plain")
r.gStyle.SetPalette(1)
r.TH1.SetDefaultSumw2(True)
r.gErrorIgnoreLevel = 2000
r.gROOT.SetBatch(True)

def quiet(function) :
    def wrapped(*args,**kwargs) :
        msg = r.RooMsgService.instance()
        cache = msg.globalKillBelow()            # get current message level
        msg.setGlobalKillBelow([r.RooFit.INFO,
                                r.RooFit.WARNING,
                                r.RooFit.ERROR,
                                r.RooFit.FATAL
                                ][-2]) # suppress messages
        val = function(*args,**kwargs)                 # call function
        msg.setGlobalKillBelow(cache)            # resume prior message level
        return val
    return wrapped


@quiet
def wimport(w, *args, **kwargs) : getattr(w, "import")(*args,**kwargs)


@quiet
def wimport_const(w, name, value) : getattr(w, "import")(r.RooConstVar(*(2*[name]+[value])))


def factory(w, command) :  w.factory(command)

def str(arg):
    ss = r.stringstream()
    arg.printStream(ss,arg.defaultPrintContents(""),arg.defaultPrintStyle(""))
    return ss.str().strip()
