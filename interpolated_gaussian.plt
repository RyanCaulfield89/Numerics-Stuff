# gnuplot plot file: interpolated_gaussian.plt
reset
set terminal x11
set xlabel 'x'
set key top left
set timestamp
plot "interpolated_gaussian.dat" using ($1):($2) title 'actual gaussian',\
     "interpolated_gaussian.dat" using ($1):($3) title 'interpolated gaussian'
set out "interpolated_gaussian.ps"
set terminal postscript color
replot
