set xlabel "Size of inner matrices products"
set title "Comparison CPU performance"
#set yrange[0:0.00001]
#set key

plot 'results_p3_Ofast.dat' u 2:4 t "Fortran function using Ofast" w lp ls 3 lt 3 pt 7 linecolor 7, 'results_p3.dat' u 2:4 t "Fortran function" w l ls 3 lt 3  linecolor 7

#set title "Comparison CPU performance"

set terminal postscript 
set output "pr3_fort.ps"
rep 
set terminal wxt
