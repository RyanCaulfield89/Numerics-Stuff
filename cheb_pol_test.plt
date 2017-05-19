# gnuplot plot file: cheb_pol_test.plt
reset
set terminal x11
set xlabel 'x'
set key top left
set timestamp
plot "cheb_pol_test.dat" using ($1):($2) title 'n = 0',\
     "cheb_pol_test.dat" using ($1):($3) title 'n = 1',\
     "cheb_pol_test.dat" using ($1):($4) title 'n = 2',\
     "cheb_pol_test.dat" using ($1):($5) title 'n = 3',\
     "cheb_pol_test.dat" using ($1):($6) title 'n = 4'
set out "cheb_pol_test.ps"
set terminal postscript color
replot
