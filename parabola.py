import math, numpy as np

class parabola(object) :
    '''Useful properties of a parabola defined by three points.'''
    
    def __init__(self, points = [] ) :
        if len(points)!=3 : raise Exception("NPointsNot3")
        self.ABC = np.linalg.solve( np.array([[x**2, 2*x, 1] for x,_ in points]),
                                    [y for _,y in points])
    def y(self,x) :
        A,B,C = self.ABC
        return A*x*x + 2*B*x + C

    def slope(self,x) :
        A,B,C = self.ABC
        return 2*(A*x + B)

    @property
    def xmin(self) :
        A,B,_ = self.ABC
        return - B / A

    @property
    def ymin(self) : return self.y(self.xmin)
    
    def dx(self,dy) :
        A,_,_ = self.ABC
        return math.sqrt(dy/A) if dy/A>0 else None
