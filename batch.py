import os
import utils

class batch(object):
    def __init__(self,cmds, nCores=4):
        self.nCores = nCores
        self.cmds = cmds
        self.default()

    def default(self):
        utils.operateOnListUsingQueue(self.nCores, utils.qWorker(os.system), zip(self.cmds))
