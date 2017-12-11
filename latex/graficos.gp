# script gnuplot
set terminal epslatex font ",10" 

#set terminal pngcairo size 1350,1262 enhanced font 'Verdana,10'
#set output 'introduction.png'

set xrange[0:12]
set yrange[-3.2:3.2]

set style line 1 lt 1 lc rgb "red" lw 2 pt 7 ps 1.2
set style line 2 lt 1 lc rgb "blue" lw 2 pt 7 ps 1.2
set style line 8 lt 2 lc rgb "black" lw 8

set border 3
set xtics nomirror
set ytics nomirror

set xlabel "$t$"
set ylabel "$\\psi(t)$"

A = 1.0
phi = 0.0
omg = 1.4

f(t) = A * cos(omg * t + phi)

set key top left
set key spacing 1.2

set output "graficoA.tex"

plot "rk2.dat" title "Runge-Kutta 2" ls 1, \
     "euler.dat" title "Euler" ls 2, \
     f(x) title "Sol. analítica" ls 8

set output "graficoB.tex"

plot "rk2B.dat" title "Runge-Kutta 2" ls 1, \
     "eulerB.dat" title "Euler" ls 2, \
     f(x) title "Sol. analítica" ls 8


