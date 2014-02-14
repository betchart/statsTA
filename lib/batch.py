import os
import __init__ as lib
import tempfile

class batch(object):

    def __init__(self,cmds, nCores=4, site=None):
        self.nCores = nCores
        self.cmds = cmds

        if site == 'ic': self.ic()
        if site == 'fnal': self.fnal()
        else: self.default()


    def default(self):
        lib.operateOnListUsingQueue(self.nCores, lib.qWorker(os.system), zip(self.cmds))


    def ic(self):
        with open('lib/ic_cmsJob.sh') as f: prelude = f.readlines()

        def sub(iJob,cmd):
            with tempfile.NamedTemporaryFile('w',prefix='job%04d_'%iJob,suffix='.sh') as tmp:
                print>>tmp, ''.join(prelude)
                print>>tmp, "cd", os.environ['PWD']
                print>>tmp, cmd
                tmp.flush()
                os.system('./lib/ic_cmsSub.sh %s'%tmp.name)
        lib.operateOnListUsingQueue(self.nCores, lib.qWorker(sub), enumerate(self.cmds))


    def fnal(self):
        with open('lib/fnal_cmsJob.sh') as f: prelude = f.readlines()
        with open('lib/fnal_cmsTemplate.condor') as f: condtemplate = f.readlines()
        delit = False
        def sub(iJob,cmd):
            with tempfile.NamedTemporaryFile(prefix='job%04d_'%iJob,suffix='.sh',dir="%s/condor"%os.environ['PWD'], delete=delit) as scr:
                print>>scr, ''.join(prelude)
                print>>scr, "cd", os.environ['PWD']
                print>>scr, cmd
                scr.flush()
                os.system('chmod +x %s'%scr.name)

                with tempfile.NamedTemporaryFile(prefix='job%04d_'%iJob,suffix='.condor',dir="%s/condor"%os.environ['PWD'], delete=delit) as cond:
                    print>>cond, ''.join(condtemplate).replace('JOBFLAG',scr.name.split('/')[-1]).replace('OUTFLAG','./')
                    cond.flush()

                    subCmd = "; ".join(["lib/fnal_cmsSub.sh %s" % cond.name])
                    os.system(subCmd)

        lib.operateOnListUsingQueue(self.nCores, lib.qWorker(sub), enumerate(self.cmds))
