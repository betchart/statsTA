import os
import utils

class batch(object):
    def __init__(self,cmds, nCores=4, site=None):
        self.nCores = nCores
        self.cmds = cmds

        if site == 'ic':
            pass
        else:
            self.default()

    def default(self):
        utils.operateOnListUsingQueue(self.nCores, utils.qWorker(os.system), zip(self.cmds))
