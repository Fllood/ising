set xlabel "Monte Carlo time t"
set ylabel "Magnetization per spin m"
set title "Magnetization versus MC time"

plot  "ising_1.dat" using 1:2 w l lt 1

pause -1

# replot on pdf-file
set term push

set terminal pdf
set output 'output_ising_1_mag.pdf'

replot
set term pop

