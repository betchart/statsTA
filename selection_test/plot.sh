set term postscript enhanced color 22 size 8in, 8in
set output 'plot.ps'
set key bottom

rab = 0.1/0.9
aA = 0.0032
aB = 0.0057

aE(aG,AA,AB) = (aG*(1-AA**2)*rab + aG*(1-AA*AB) + (AB-AA)) / ( (1-AA**2)*rab + (1-AA*aB) + aG*(AB-AA))

L(RAB) = RAB / (1+RAB)

aD(aG,AA,AB) = aG * (L(rab)*AA + (1-L(rab))*AB) / AA

set grid
set xlabel 'A@_c^y (%) in lepton+jets selection '
set ylabel 'extrapolated total A@_c^y (%)'

lim = 1
plot [-lim:lim] 100*aD(x/100, aA, aB) title 'POWHEG 2-bin templates', 100*aE(x/100, aA, aB) title 'POWHEG 2-bin unfolding' lc 3, "<echo '0.32 0.56'" title 'POWHEG' lc 2, x title 'bisector'
