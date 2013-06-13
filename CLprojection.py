import numpy as np
import math

def oneSigmaCLprojection(sigmas2):
    A = np.linalg.inv(sigmas2) * 2 * 1.14
    return math.sqrt(A[1,1] / np.linalg.det(A))
