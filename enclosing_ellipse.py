# Implementation of algorithms described in slides by Charles F. Van Loan, Cornell University

import numpy as np
import math
from ellipse import ellipse
from goldenSection import goldenSectionSearch,resphi

def s_given_foci(pointset, foci):
    '''String length given foci.'''
    return max(sum(math.sqrt(np.dot(*(2*[f-p]))) for f in foci)
               for p in pointset)


def d_s_given_center_tilt(pointset, c, t):
    '''Foci spacing and string length given center and tilt.'''

    dlo,dhi = (0, 2*max([math.sqrt(np.dot(*(2*[c-p]))) for p in pointset]))
    def s_given_d(d):
        delta = 0.5*d*np.array([math.cos(t),math.sin(t)])
        foci = [c+delta,c-delta]
        return s_given_foci(pointset, foci)
        
    def area(d): 
        s = s_given_d(d)
        a = 0.25*math.pi*s*math.sqrt(s**2-d**2)
        return a
    b = dlo+resphi*(dhi-dlo)
    D = goldenSectionSearch(area, zip((dlo,b,dhi),3*[None]), 0.0001*dhi)
    return D , s_given_d(D)


def tilt_d_s_given_center(pointset, c):
    '''Tilt, foci spacing, and string length given center. '''
    tlo,thi = 0,math.pi
    def area(t):
        d,s = d_s_given_center_tilt(pointset, c, t)
        return 0.25*math.pi*s*math.sqrt(s**2-d**2)
    b = tlo+resphi*(thi-tlo)
    T = goldenSectionSearch(area, zip((tlo,b,thi),3*[None]), 0.05*thi)
    return T, + d_s_given_center_tilt(c,T)


def enclosing_ellipse(pointset, center):
    #t,d,s = tilt_d_s_given_center(pointset,np.array(center))
    pmax = max(pointset, key=lambda p: np.dot(*(2*[np.array(center)-p])))
    t = math.atan2((pmax[1]-center[1]),(pmax[0]-center[0]))
    d,s = d_s_given_center_tilt(pointset, np.array(center), t)
    a2 = 0.25 * s**2
    b2 = 0.25 * (s**2-d**2)
    c = math.cos(t)
    s = math.sin(t)
    R = np.array([[c,-s],[s,c]])
    sigmas2 = np.dot(R,[[a2,0],[0,b2]]).dot(R.T) 
    return ellipse( mean = center, sigmas2 = sigmas2 )
