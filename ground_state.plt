# gnuplot plot file: ground_state.plt
reset
set terminal x11
set xlabel 'x'
set ylabel 'psi(x)'
set key top left
set timestamp
plot "ground_state.dat" using ($1):($2) title 'psi_0(x)'
set out "ground_state.ps"
set terminal postscript color
replot
