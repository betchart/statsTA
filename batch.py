import os
import lib
import tempfile

class batch(object):
    def __init__(self,cmds, nCores=4, site=None):
        self.nCores = nCores
        self.cmds = cmds

        if site == 'ic':
            localdir =  '/'.join(__file__.split('/')[:-1])
            prelude = '''
export SCRAM_ARCH=slc5_amd64_gcc472
source /vols/cms/grid/setup.sh
cd /home/hep/bbetchar/work/CMSSW_6_2_0_patch1/src
eval `scram runtime -sh`
cd %s
'''%localdir
            def sub(iJob,cmd):
                with tempfile.NamedTemporaryFile('w',prefix='job%04d_'%iJob,suffix='.sh') as tmp:
                    print>>tmp, prelude
                    print>>tmp, cmd
                    tmp.flush()
                    os.system('./ic_cmsSub.sh %s'%tmp.name)
            lib.operateOnListUsingQueue(self.nCores, lib.qWorker(sub), enumerate(self.cmds))
        else:
            self.default()

    def default(self):
        lib.operateOnListUsingQueue(self.nCores, lib.qWorker(os.system), zip(self.cmds))
