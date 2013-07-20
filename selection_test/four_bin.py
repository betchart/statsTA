import math
import sys
import numpy as np

# POWHEG 2-bin parameters
a_s = 0.0032
a_t = 0.0056
eff = 0.1
R = eff / (1-eff)
a_sbar = R*(a_t/eff - a_s)

class twoBinSlopes(object):
    def __init__(self,g):
        eps = eff / g
        deps = math.sqrt(1-g)
        da = (a_t/a_s-1) / deps
        a_in = a_s * (1-da)
        a_out = a_s * (1+da)
        eff_in = eps * (1+deps)
        eff_out = eps * (1-deps)
        for item in ['g','eps','deps','da','a_in','a_out','eff_in','eff_out']:
            setattr(self,item,eval(item))

    def As(self,fin,fout):
        ain = fin*self.a_in
        aout = fout*self.a_out
        return 0.5 * (ain+aout)

    def slopeUnfold(self,fin,fout):
        return ( (self.a_tot(fin*self.a_in,fout*self.a_out) - a_t) / 
                 (self.unfold2(self.As(fin,fout)) - a_t) )

    def slopeTemplate(self,fin,fout):
        return ( (self.a_tot(fin*self.a_in,fout*self.a_out) - a_t) / 
                 (self.template2(self.As(fin,fout)) - a_t) )

    def a_tot(self, ain, aout):
        return 100 * (ain/self.eff_in + aout/self.eff_out) / (1/self.eff_in + 1/self.eff_out)
    
    def unfold2(self, As):
        num = As * (1-a_s**2)*R + As*(1-a_s*a_sbar) + a_sbar-a_s
        den = (1-a_s**2)*R + 1-a_s*a_sbar + As*(a_sbar-a_s)
        return 100* num / den
    
    def template2(self, As):
        return 100* As * (a_t / a_s)
        



def line(g):
    slopes = twoBinSlopes(g)
    return '\t'.join([str(i) for i in [
        slopes.eff_out/slopes.eff_in,
        slopes.slopeUnfold(2,1),
        slopes.slopeTemplate(2,1),
        slopes.slopeUnfold(1,2),
        slopes.slopeTemplate(1,2),
        slopes.slopeUnfold(2,2),
        slopes.slopeTemplate(2,2),
        g,
        100*slopes.eff_in,
        100*slopes.eff_out,
        100*slopes.a_in,
        100*slopes.a_out
    ]])

with open('slopes.txt','w') as f:
    for g in np.arange(0.19,0.99,0.005):
        print>>f, line(g)

