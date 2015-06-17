set title "Eigenmode"

set xlabel "x"

set ylabel "y"

set zlabel "z"

set hidden3d

set grid

set terminal pdf

set output 'op-eigenmode.pdf'

splot 'eigenmode.txt' title 'k(x,y)^2' pt 6 ps 0.1

set title "Gaussian-ksq"

set xlabel "x"

set ylabel "y"

set zlabel "z"

set hidden3d

set grid

set terminal pdf

set output 'op-ksq.pdf'

splot 'ksq.txt' title 'k(x,y)^2' pt 6 ps 0.1
