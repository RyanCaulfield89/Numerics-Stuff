# gnuplot plot file: interpolated_rational.plt
reset
set terminal x11
set xlabel 'x'
set yrange[-1:1]
set key top left
set timestamp
plot "interpolated_rational.dat" using ($1):($2) with lines title 'exact',\
     "interpolated_rational.dat" using ($1):($3) with lines title 'chebyshev nodes',\
     "interpolated_rational.dat" using ($1):($4) with lines title 'even spacing'
set out "interpolated_rational.ps"
set terminal postscript color
replot
