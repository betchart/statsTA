import sys
import utils

pat = sys.argv[1]

nums = [int(f.replace(pat,'').replace('.root','')) for f in utils.getCommandOutput("ls %s*.root"%pat)["stdout"].split()]
print ','.join([str(s) for s in sorted(set(range(max(nums))) - set(nums))])
