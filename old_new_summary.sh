set term postscript landscape enhanced color
set output 'graphics/old_new.ps'

set xlabel 'Reconstruction A@_c^y   (%)'
set ylabel 'Reconstruction A@_c^{/Symbol f}   (%)'

set grid

set size square
set key out left box

set title 'POWHEG Templates: q@^{\261}q --> t@^{\261}t'
plot [-2:2] [-2:2] \
  'ttqq_points.txt' u 1:2 w lines lt 1 lc 1 title 'original e',\
  '' u 3:4 w lines lt 2 lc 1 title 'additional e',\
  '' u 5:6 w lines lt 1 lw 3 lc 1 title 'total e',\
  '' u 7:8 w lines lt 1 lc 3 title 'original {/Symbol m}',\
  '' u 9:10 w lines lt 2 lc 3 title 'additional {/Symbol m}',\
  '' u 11:12 w lines lt 1 lw 3 lc 3 title 'total {/Symbol m}',\
  'ttqg_points.txt' u 1:2 w lines lt 3 lc 9 title 'other',\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle,\
    'ttag_points.txt' u 1:2 w lines lt 3 lc 9 notitle,\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle


set title 'POWHEG Templates: qg --> t@^{\261}t'
plot [-2:2] [-2:2] \
  'ttqg_points.txt' u 1:2 w lines lt 1 lc 1 title 'original e',\
  '' u 3:4 w lines lt 2 lc 1 title 'additional e',\
  '' u 5:6 w lines lt 1 lw 3 lc 1 title 'total e',\
  '' u 7:8 w lines lt 1 lc 3 title 'original {/Symbol m}',\
  '' u 9:10 w lines lt 2 lc 3 title 'additional {/Symbol m}',\
  '' u 11:12 w lines lt 1 lw 3 lc 3 title 'total {/Symbol m}',\
    'ttag_points.txt' u 1:2 w lines lt 3 lc 9 title 'other',\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle,\
  'ttqq_points.txt' u 1:2 w lines lt 3 lc 9 notitle,\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle

set title 'POWHEG Templates: @^{\261}qg --> t@^{\261}t'
plot [-2:2] [-2:2] \
  'ttag_points.txt' u 1:2 w lines lt 1 lc 1 title 'original e',\
  '' u 3:4 w lines lt 2 lc 1 title 'additional e',\
  '' u 5:6 w lines lt 1 lw 3 lc 1 title 'total e',\
  '' u 7:8 w lines lt 1 lc 3 title 'original {/Symbol m}',\
  '' u 9:10 w lines lt 2 lc 3 title 'additional {/Symbol m}',\
  '' u 11:12 w lines lt 1 lw 3 lc 3 title 'total {/Symbol m}',\
  'ttqg_points.txt' u 1:2 w lines lt 3 lc 9 title 'other',\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle,\
  'ttqq_points.txt' u 1:2 w lines lt 3 lc 9 notitle,\
  '' u 3:4 w lines lt 2 lc 9 notitle,\
  '' u 5:6 w lines lt 3 lw 3 lc 9 notitle,\
  '' u 7:8 w lines lt 3 lc 9 notitle,\
  '' u 9:10 w lines lt 2 lc 9 notitle,\
  '' u 11:12 w lines lt 3 lw 3 lc 9 notitle
