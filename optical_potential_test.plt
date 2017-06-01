# gnuplot plot file: optical_potential_test.plt
reset
set terminal x11
set xlabel 'r'
set ylabel 'V(r)'
set key top left
set timestamp
plot "optical_potential_test.dat" using ($1):($2) title 'Real(V(r))',\
     "optical_potential_test.dat" using ($1):($3) title 'Imag(V(r))'
set out "optical_potential_test.ps"
set terminal postscript color
replot
