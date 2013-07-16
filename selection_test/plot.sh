set term postscript enhanced color 22 size 8in, 5in
set output 'plot.ps'
set key bottom

rab = 0.052
aa = 0.0032
ab = 0.0057

unf(f) = f * aa * (rab + (1-aa*ab+(ab/aa-1)/f) / (1-aa**2)) / (rab + (1-aa*ab+f*aa*(ab-aa)) / (1-aa**2))

L(RAB) = RAB / (1+RAB)

alpha(f) = f * (L(rab)*aa + (1-L(rab))*ab)

set grid
set xlabel 'F'
set ylabel 'A@_c^y'


plot [-3:3] alpha(x) title 'true asymmetry',unf(x) title 'calculation via unfolding' lc 3, "<echo '1 0.0056'" title 'POWHEG' lc 2
