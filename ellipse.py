import math
import numpy as np
import itertools

class ellipse(object):

    def __init__(self, matrix=None, mean=None, sigmas2=None):

        assert not bool(mean)^(sigmas2!=None)         # both or neither
        assert bool(matrix)^bool(mean and sigmas2!=None)  # pick one

        for item in ['matrix', 'mean', 'sigmas2']: setattr(self, item, eval(item))

        if not matrix: 
            self.fillMatrix()
        else: 
            self.fillMean()
            self.fillSigmas()
        self.fillParametric()

    def fillMatrix(self):
        E,R = np.linalg.eig(self.sigmas2)
        twobytwo = R.dot(np.diag([1/e for e in E])).dot(R.T)
        shape = np.hstack([np.vstack([twobytwo, [0,0]]), 
                           [[0],[0],[-1]]])
        trans = np.eye(3)
        trans[0,2] = -self.mean[0]
        trans[1,2] = -self.mean[1]
        self.matrix = trans.T.dot(shape).dot(trans)

    def fillMean(self):
        pass

    def fillSigmas(self):
        pass


    @staticmethod
    def getParametric(matrix):
        E,R = np.linalg.eig(matrix)
        for p in itertools.permutations(range(3),3):
            d = np.diag(np.sqrt(np.abs(E[p,:]))).dot(R[:,p].T)
            if np.allclose(matrix, d.T.dot(np.diag([1,1,-1]).dot(d))):
                return np.linalg.inv(d) 

    def fillParametric(self):
        self.parametric = self.getParametric(self.matrix)

    def eval(self,t):
        point = self.parametric.dot([math.cos(t), math.sin(t), 1])
        return point[0]/point[2], point[1]/point[2]
