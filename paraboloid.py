import math, numpy as np

class paraboloid(object) :
    '''Useful properties of a paraboloid defined by six points.'''
    
    def __init__(self, points = [] ) :
        if len(points)!=6 : raise "NPointsNot6"
        self.ABCDEF = np.linalg.solve( np.array([[x**2, y**2, x*y, x, y, 1] for x,y,_ in points]),
                                       [z for _,z in points])
    def z(self,x,y) :
        A,B,C,D,E,F = self.ABCDEF
        return A*x*x + B*y*y + C*x*y + D*x + E*y + F

    def gradient(self,x,y) :
        A,B,C,D,E,F = self.ABCDEF
        return (2*A*x + C*y + D,
                2*B*y + C*x + E)

    @property
    def xymin(self) :
        A,B,C,D,E,_ = self.ABCDEF
        x = (C*E - 2*B*D) / (4*A*B - C**2)
        y = -0.5 * (C*x + E) / B
        return x,y

    @property
    def zmin(self) : return self.z(*self.xymin)
    
    def contour_matrix(self,dz) :
        return None

