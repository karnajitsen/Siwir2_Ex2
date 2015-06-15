set title "Multigrid solution of Laplace equation in 2 dimensions"
set xlabel "x"
set ylabel "y"
set zlabel "z"
set hidden3d
set grid
set terminal pdf
set output 'output.pdf'
splot 'output' title 'k(x,y)^2' pt 6 ps 0.1
