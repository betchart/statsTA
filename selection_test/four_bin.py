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
        ain = fin * self.a_in
        aout = fout * self.a_out
        return 0.5 * (ain + aout)

    def slopeUnfold(self,fin,fout):
        num = (self.unfold2(self.As(fin,fout)) - 100*a_t)
        den = (self.a_tot(fin*self.a_in,fout*self.a_out) - 100*a_t)
        return num/den

    def slopeTemplate(self,fin,fout):
        num = (self.template2(self.As(fin,fout)) - 100*a_t)
        den = (self.a_tot(fin*self.a_in,fout*self.a_out) - 100*a_t)
        return num/den

    def a_tot(self, ain, aout):
        return 100 * (ain/self.eff_in + aout/self.eff_out) / (1/self.eff_in + 1/self.eff_out)
    
    def unfold2(self, As):
        num = As * (1-a_s**2)*R + As*(1-a_s*a_sbar) + a_sbar-a_s
        den = (1-a_s**2)*R + 1-a_s*a_sbar + As*(a_sbar-a_s)
        return 100* num / den
    
    def template2(self, As):
        return 100* As * (a_t / a_s)

    def percentMistakes(self,fin,fout):
        true = self.a_tot(fin*self.a_in,fout*self.a_out)
        As = self.As(fin,fout)
        unfold = self.unfold2(As)
        template = self.template2(As)
        return [(unfold-true)/(100*a_t), (template-true)/(100*a_t)]
        


def line(g):
    scale = 2
    slopes = twoBinSlopes(g)
    return ' '.join([str(i).ljust(15) for i in [
        100*slopes.eff_in,
        slopes.slopeUnfold(scale,1),
        slopes.slopeTemplate(scale,1),
        slopes.slopeUnfold(1,scale),
        slopes.slopeTemplate(1,scale),
        slopes.slopeUnfold(scale,scale),
        slopes.slopeTemplate(scale,scale),
        g,
        100*slopes.eff_in,
        100*slopes.eff_out,
        100*slopes.a_in,
        100*slopes.a_out
    ]+(slopes.percentMistakes(-1,1)+
       slopes.percentMistakes(1,-1)+
       slopes.percentMistakes(-1,-1))
   ])

with open('slopes.txt','w') as f:
    print>>f, ' '.join([i.ljust(15) for i in 'eff_in    sl_u_in    sl_t_in    sl_u_out    sl_t_out    sl_u_both    sl_t_both    g      eff_in    eff_out    a_in    a_out'.split()])
    for g in np.arange(0.19,0.99,0.005):
        print>>f, line(g)

