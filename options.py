from optparse import OptionParser
parser = OptionParser("usage: %prog [options]")
def argOrTrue(option, opt, value, parser) :
    peek = next(iter(parser.rargs),None)
    if peek and peek[0]!='-' : del parser.rargs[0]
    setattr(parser.values, option.dest, peek if peek and peek[0]!='-' else True)

parser.add_option("--partitions", dest="partitions", default=None, metavar='p1,p2,...', action="callback", callback=argOrTrue, help='specify list of partitions, or list partitions')
parser.add_option("--systematics", dest="systematics", default=None, metavar='sys1,sys2,...', action="callback", callback=argOrTrue, help='specify list of systematics, or list systematics')
parser.add_option("--ensembles", dest="ensembles", default=None, metavar='lo:hi', help='ensembles by slice notation')
parser.add_option("--templates", dest="templates", default=None, metavar='lo:hi', help='templates by slice notation')
parser.add_option("--visualize", dest='visualize', default=False, action='store_true', help='project the fits')

def opts() :
    options,args = parser.parse_args()
    if options.partitions==None:
        parser.print_help()
        exit()
    return options

def default(options = []) :
    options,args = parser.parse_args(options)
    return options
