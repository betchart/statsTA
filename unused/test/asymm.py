import ROOT as r

w = r.RooWorkspace()
w.factory("x[0,-1,1]")
w.factory("Gaussian::g(x, 0, 1)")
w.factory("EXPR::p_symm('3./8*(1+x**2)',x)")
w.factory("EXPR::p_both('3./8*(1+x**2)-0.1*x',x)")
w.factory("SUM::p( alpha[-3,-5,5] * p_both, p_symm)")

xframe = w.arg('x').frame()
w.pdf('p').plotOn(xframe)
print w.pdf('p').getNorm(w.argSet('x'))
xframe.Draw()

w.Print()

raw_input()
