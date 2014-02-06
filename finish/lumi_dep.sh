set term postscript enhanced 22 size 7in, 7in
set output 'graphics/lumi_dep.ps'

set logscale xy
set xlabel 'luminosity / 19.59fb^{-1}'
set ylabel 'factor on expected stat. uncertainty'
set grid
set key spacing 1.7

plot \
    'output/lumi_dep.txt' \
    u    1:2 pt 7  title 'A@_c^{(qq)}', \
    '' u 1:3       title 'A@_c^{(qg)}', \
    1./sqrt(x)     title 'x^{-1/2}'
