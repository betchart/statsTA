import math, numpy as np
import itertools

class paraboloid(object):
    '''Useful properties of a vertical paraboloid defined by six points.'''
    
    def __init__(self, points=[]):
        if len(points)!=6 : raise Exception("NPointsNot6")
        self.ABCDEF = np.linalg.solve( np.array([[x**2, y**2, 2*x*y, 2*x, 2*y, 1] for x,y,_ in points]),
                                       [z for _,_,z in points])
    def z(self, x, y):
        A,B,C,D,E,F = self.ABCDEF
        return A*x*x + B*y*y + 2*C*x*y + 2*D*x + 2*E*y + F

    def gradient(self, x, y):
        A,B,C,D,E,F = self.ABCDEF
        return (2*(A*x + C*y + D),
                2*(B*y + C*x + E))

    @property
    def xymin(self):
        A,B,C,D,E,_ = self.ABCDEF
        x = (C*E - B*D) / (A*B - C**2)
        y = -(C*x + E) / B
        return x,y

    @property
    def zmin(self) : return self.z(*self.xymin)
    
    def dxy(self,dz):
        A,B,C,_,_,_ = self.ABCDEF
        ellipse = np.array([[A,C],[C,B]]) / dz
        return np.linalg.inv(ellipse)

    def ellipse(self,dz):
        A,B,C,D,E,F = self.ABCDEF
        return np.array([[A,C,D],
                         [C,B,E],
                         [D,E,F-dz]])

    def parametricEllipse(self,dz):
        '''Matrix multiplying [cos(t),sin(t),1] to yield points on the ellipse.'''
        ellipse = self.ellipse(dz)
        E,R = np.linalg.eig(ellipse)
        for p in itertools.permutations(range(3),3):
            d = np.diag(np.sqrt(np.abs(E[p,:]))).dot(R[:,p].T)
            if np.allclose(matrix, d.T.dot(np.diag([1,1,-1]).dot(d))):
                return np.linalg.inv(d) 

    def major_angle(self):
        A,B,C,_,_,_ = self.ABCDEF
        ellipse = np.array([[A,C],[C,B]])
        E,R = np.linalg.eig(ellipse)
        return math.atan2(*ellipse[0])
        
