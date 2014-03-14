#!/usr/bin/env python

import ROOT as r
r.gROOT.SetBatch(1)
import numpy as np
import math
from lib.ellipse import ellipse
from lib.__autoBook__ import autoBook

r.gStyle.SetOptStat(111111)

class fitresult(object):
    def __init__(self,fname, getfirstentry = False):
        self.tfile = r.TFile.Open(fname)
        self.tree = self.tfile.Get('fitresult')
        if getfirstentry: self.tree.GetEntry(0)

oneSigmaN2LL = 1.14 * 2


class calibrated(object):

    def __init__(self, fname):
        self.cal = fitresult(fname)
        self.book = autoBook(fname.split('/')[-1])
        for e in self.cal.tree:
            Ac_y = self.Ac_y_tt(e)
            Ac_y_gen = e.gen_Ac_y_ttalt
            Ac_y_sig = self.sigmas_Ac_y_tt(e)
            self.book.fill( Ac_y - Ac_y_gen, 'delta_Ac_y', 100, -0.05, 0.05)
            self.book.fill( Ac_y_sig, 'Ac_y_sig', 100, 0, 0.006)
            self.book.fill( (Ac_y - Ac_y_gen) / Ac_y_sig, 'Ac_y_pull', 100, -5, 5)
        print Ac_y_gen

        self.Draw()

    def Ac_y_tt(self, e):
        return sum([e.Ac_y_ttgg*e.f_gg,
                    e.Ac_y_ttqg*e.falphaT,
                    e.Ac_y_ttqq*e.falphaL,
                    e.Ac_y_ttag*e.f_ag*e.alphaT])

    def Ac_phi_tt(self, e):
        return sum([e.Ac_phi_ttgg*e.f_gg,
                    e.Ac_phi_ttqg*e.falphaT,
                    e.Ac_phi_ttqq*e.falphaL,
                    e.Ac_phi_ttag*e.f_ag*e.alphaT])

    def sigmas_Ac_y_tt(self, e):
        R = np.array([[1,1],[-1,1]]) # rotate pi/4, except also scale by sqrt(2)
        sigmas = np.array([[e.fitXX,e.fitXY],[e.fitXY,e.fitYY]]) / oneSigmaN2LL
        return math.sqrt(R.dot(sigmas).dot(R.T)[0,0])

    def Draw(self):
        outName = self.book.title.rstrip('.root') + '.pdf'
        c = r.TCanvas()
        c.Print(outName + '[')
        for hname, h in self.book.items():
            h.Draw()
            c.Print(outName)
        c.Print(outName + ']')
        print "Output in", outName


if __name__=='__main__':
    import sys
    if len(sys.argv) < 2 or sys.argv[1][-5:]!='.root':
        print "Usage: calibrated.py <calibration_fitresults.root>"
        exit()
    calibrated(sys.argv[1])
