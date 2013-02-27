import math, numpy as np

class parabola(object) :
    '''Useful properties of a parabola defined by three points.'''
    
    def __init__(self, points = [] ) :
        if len(points)!=3 : raise "NPointsNot3"
        self.ABC = np.linalg.solve( np.array([[x**2, x, 1] for x,_ in points]),
                                    [y for _,y in points])
    def y(self,x) :
        A,B,C = self.ABC
        return A*x*x + B*x + C

    def slope(self,x) :
        A,B,C = self.ABC
        return 2*A*x + B

    @property
    def xmin(self) :
        A,B,C = self.ABC
        return -0.5 * B / A

    @property
    def ymin(self) : return self.y(self.xmin)
    
    def dx(self,dy) :
        A,B,C = self.ABC
        return 0.5 * math.sqrt(B*B - 4*A*(C-dy)) / A

